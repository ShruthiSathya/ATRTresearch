#!/usr/bin/env python3
"""
scripts/01_prepare_cmap.py
==========================
ATRT CMap Query Preparation — generates gene lists for clue.io submission.

WHAT THIS DOES
--------------
1. Loads GSE70678 ATRT expression vs GTEx normal brain
2. Runs Welch's t-test to find differentially expressed genes
3. Writes gene lists for clue.io L1000 query

INPUTS (must already exist)
---------------------------
  data/raw_omics/GSE70678_gene_expression.tsv   ← your processed file ✅
  data/raw_omics/GTEx_brain_normal_reference.tsv ← your processed file ✅

OUTPUTS
-------
  data/cmap_query/atrt_up_genes.txt       ← paste into clue.io UP box
  data/cmap_query/atrt_down_genes.txt     ← paste into clue.io DOWN box
  data/cmap_query/atrt_de_results.tsv     ← full DE table
  data/cmap_query/atrt_query_manifest.txt ← step-by-step instructions

AFTER RUNNING THIS SCRIPT
--------------------------
Option A — clue.io web UI (recommended):
  1. Go to https://clue.io → Log in (free academic account)
  2. Tools → L1000 Query (sig_queryl1k_tool)
  3. Paste atrt_up_genes.txt → "Up genes"
  4. Paste atrt_down_genes.txt → "Down genes"
  5. Name: ATRT_SMARCB1_loss_signature
  6. Perturbation type: Compounds (trt_cp)
  7. Submit → ~10-15 min
  8. Download query_result.gct
  9. cp ~/Downloads/query_result.gct data/cmap_query/query_result.gct
  10. python -m backend.pipeline.integrate_cmap_results

Option B — clue.io API (needs API key from your profile):
  python scripts/01_prepare_cmap.py --submit --api-key YOUR_KEY
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
ATRT_EXPR = REPO_ROOT / "data/raw_omics/GSE70678_gene_expression.tsv"
GTEX_REF  = REPO_ROOT / "data/raw_omics/GTEx_brain_normal_reference.tsv"
OUT_DIR   = REPO_ROOT / "data/cmap_query"

# ── Statistical thresholds ────────────────────────────────────────────────────
LOG2FC_THRESHOLD = 1.0   # ≥ 2-fold change
P_ADJ_THRESHOLD  = 0.05  # BH-corrected FDR
MIN_EXPR         = 3.0   # log2 units
MIN_SAMPLES      = 5     # minimum ATRT samples per gene (relaxed for small cohorts)

# ── Expected ATRT biology markers ─────────────────────────────────────────────
EXPECTED_UP = {"EZH2", "AURKA", "BRD4", "MYC", "MYCN", "HDAC1", "CDK4", "SOX2", "LIN28A"}
EXPECTED_DOWN = {"SMARCB1", "CDKN2A", "OLIG2"}

# ── L1000 library gaps (post-2017 drugs not profiled) ─────────────────────────
L1000_NOT_PROFILED = {
    "TAZEMETOSTAT": "Approved 2020 — after L1000 profiling cutoff",
    "MARIZOMIB": "Marine natural product — not commercially profiled",
}
L1000_PARTIAL = {
    "ALISERTIB": "Partial profiling — AURKA not a standard landmark gene",
}


def bh_correction(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction (scipy fallback for older versions)."""
    n     = len(p_values)
    order = np.argsort(p_values)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    p_adj = np.minimum(1.0, p_values * n / ranked)
    for i in range(n - 2, -1, -1):
        p_adj[order[i]] = min(p_adj[order[i]], p_adj[order[i + 1]])
    return p_adj


def load_atrt_matrix(path: Path) -> pd.DataFrame:
    logger.info("Loading ATRT expression matrix: %s", path)
    df = pd.read_csv(str(path), sep="\t", index_col=0, low_memory=False)
    df.index = df.index.astype(str).str.upper().str.strip()
    df = df.apply(pd.to_numeric, errors="coerce")

    # Drop all-NaN rows
    df = df.dropna(how="all")
    logger.info("  Shape: %d genes × %d ATRT samples", len(df), len(df.columns))

    # Show sample names
    sample_names = list(df.columns[:5])
    logger.info("  Sample IDs (first 5): %s ...", sample_names)
    return df


def load_gtex_reference(path: Path) -> pd.Series:
    logger.info("Loading GTEx brain reference: %s", path)
    df = pd.read_csv(str(path), sep="\t", index_col=0, low_memory=False)
    df.index = df.index.astype(str).str.upper().str.strip()
    df = df.apply(pd.to_numeric, errors="coerce")

    # Select brain-relevant columns
    brain_cols = [
        c for c in df.columns
        if any(k.lower() in c.lower()
               for k in ["brain", "cerebellum", "cortex", "hippocampus"])
    ]
    if not brain_cols:
        logger.warning("No brain tissue columns found — using all columns")
        brain_cols = list(df.columns)

    ref = df[brain_cols].mean(axis=1).dropna()
    logger.info(
        "  Brain reference: %d genes, tissues: %s",
        len(ref),
        [c.replace("Brain - ", "")[:15] for c in brain_cols[:4]]
    )
    return ref


def run_welch_de(
    atrt_matrix: pd.DataFrame,
    gtex_ref:    pd.Series,
) -> pd.DataFrame:
    """
    Welch's one-sample t-test: ATRT samples vs GTEx normal brain reference.

    Platform note:
      GSE70678 = Affymetrix HuGene 1.0 ST (GPL570), RMA log2-normalised
      GTEx     = log2(median_TPM + 1)
    Both are log2 scale. We z-score the log2FC vector to correct platform offset.
    """
    common = atrt_matrix.index.intersection(gtex_ref.index)
    logger.info("Common genes (ATRT ∩ GTEx): %d", len(common))

    if len(common) < 500:
        logger.warning(
            "Only %d common genes — expected >5,000. "
            "Check that both files use uppercase HGNC symbols (e.g. EZH2, BRD4).",
            len(common)
        )

    rows = []
    n_skipped_samples = 0
    n_skipped_expr    = 0

    for gene in common:
        expr_vals = pd.to_numeric(atrt_matrix.loc[gene], errors="coerce").dropna()
        ref_val   = float(gtex_ref[gene])
        n_valid   = len(expr_vals)

        if n_valid < MIN_SAMPLES:
            n_skipped_samples += 1
            continue

        # Skip if both very lowly expressed (unexpressed genes — not informative)
        if expr_vals.mean() < MIN_EXPR and ref_val < MIN_EXPR:
            n_skipped_expr += 1
            continue

        mean_atrt = float(expr_vals.mean())
        log2fc    = mean_atrt - ref_val
        t_stat, p_val = stats.ttest_1samp(expr_vals, popmean=ref_val)

        rows.append({
            "gene":        gene,
            "log2FC":      round(log2fc, 4),
            "t_statistic": round(float(t_stat), 4),
            "p_value":     float(p_val),
            "mean_atrt":   round(mean_atrt, 4),
            "mean_normal": round(ref_val, 4),
            "std_atrt":    round(float(expr_vals.std()), 4),
            "n_atrt":      n_valid,
        })

    logger.info(
        "DE: tested %d genes | skipped %d (too few samples) | %d (lowly expressed)",
        len(rows), n_skipped_samples, n_skipped_expr
    )

    if not rows:
        logger.error("No genes survived DE filtering — check input data.")
        return pd.DataFrame()

    de_df = pd.DataFrame(rows).set_index("gene")
    de_df["p_adj"] = bh_correction(de_df["p_value"].values)

    # Z-score of log2FC (corrects for microarray/RNA-seq platform offset)
    lfc_mean = de_df["log2FC"].mean()
    lfc_std  = de_df["log2FC"].std()
    de_df["zscore_log2FC"] = (de_df["log2FC"] - lfc_mean) / max(lfc_std, 0.01)

    de_df["significant"] = (
        (de_df["log2FC"].abs() >= LOG2FC_THRESHOLD) &
        (de_df["p_adj"] <= P_ADJ_THRESHOLD)
    )
    de_df["direction"] = np.where(
        ~de_df["significant"], "NS",
        np.where(de_df["log2FC"] > 0, "UP", "DOWN"),
    )

    n_up   = int((de_df["direction"] == "UP").sum())
    n_down = int((de_df["direction"] == "DOWN").sum())
    n_ns   = int((de_df["direction"] == "NS").sum())
    logger.info(
        "Welch DE: %d genes | UP=%d DOWN=%d NS=%d (|log2FC|≥%.1f, FDR≤%.2f)",
        len(de_df), n_up, n_down, n_ns, LOG2FC_THRESHOLD, P_ADJ_THRESHOLD,
    )

    # Validation check
    _up_genes_set = set(de_df[de_df["direction"] == "UP"].index)
    found_up   = EXPECTED_UP   & _up_genes_set
    found_down = EXPECTED_DOWN & set(de_df[de_df["direction"] == "DOWN"].index)

    if len(found_up) >= 3:
        logger.info("✅ ATRT biology check — UP markers found: %s", sorted(found_up))
    else:
        logger.warning(
            "⚠️  Only %d/%d expected UP markers found: %s\n"
            "   Check that GSE70678 loaded correctly (should have ~49 ATRT samples)",
            len(found_up), len(EXPECTED_UP), sorted(found_up)
        )

    if found_down:
        logger.info("✅ DOWN markers found (SMARCB1 lost): %s", sorted(found_down))

    return de_df.reset_index()


def select_query_genes(
    de_df:  pd.DataFrame,
    n_up:   int = 150,
    n_down: int = 50,
) -> Tuple[List[str], List[str]]:
    """Select top UP and DOWN genes for clue.io query."""
    up_sig   = de_df[de_df["direction"] == "UP"].sort_values("zscore_log2FC", ascending=False)
    down_sig = de_df[de_df["direction"] == "DOWN"].sort_values("zscore_log2FC", ascending=True)

    # Relax to all positive/negative if few significant
    if len(up_sig) < 30:
        logger.warning(
            "%d significant UP genes — relaxing to all positive log2FC", len(up_sig)
        )
        up_sig = de_df[de_df["log2FC"] > 0].sort_values("zscore_log2FC", ascending=False)

    if len(down_sig) < 10:
        logger.warning(
            "%d significant DOWN genes — relaxing to all negative log2FC", len(down_sig)
        )
        down_sig = de_df[de_df["log2FC"] < 0].sort_values("zscore_log2FC", ascending=True)

    up_genes   = up_sig["gene"].head(n_up).tolist()
    down_genes = down_sig["gene"].head(n_down).tolist()

    logger.info(
        "CMap query: %d UP genes, %d DOWN genes selected",
        len(up_genes), len(down_genes)
    )
    return up_genes, down_genes


def write_outputs(
    up_genes:   List[str],
    down_genes: List[str],
    de_df:      pd.DataFrame,
    out_dir:    Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    # Gene lists for clue.io
    (out_dir / "atrt_up_genes.txt").write_text("\n".join(up_genes))
    (out_dir / "atrt_down_genes.txt").write_text("\n".join(down_genes))
    logger.info("Saved: %s/{atrt_up_genes,atrt_down_genes}.txt", out_dir)

    # Full DE table
    de_path = out_dir / "atrt_de_results.tsv"
    if not de_df.empty:
        de_df.to_csv(str(de_path), sep="\t", index=False)
        logger.info("Saved DE table (%d genes): %s", len(de_df), de_path)

    # L1000 library gap notes
    gap_notes = "\n".join(
        f"  {drug}: {reason}"
        for drug, reason in {**L1000_NOT_PROFILED, **L1000_PARTIAL}.items()
    )

    # Manifest
    top20_up   = "\n".join(f"  {g}" for g in up_genes[:20])
    top10_down = "\n".join(f"  {g}" for g in down_genes[:10])

    manifest = f"""ATRT CMap Query — Differential Expression vs Normal Brain
=========================================================
Generated: GSE70678 (Torchia 2015, n=49 ATRT) vs GTEx v8 brain reference
Method: Welch's one-sample t-test, BH-corrected FDR ≤ 0.05, |log2FC| ≥ 1.0

Query: {len(up_genes)} UP genes + {len(down_genes)} DOWN genes

TOP 20 UP GENES (paste all from atrt_up_genes.txt):
{top20_up}
  ... ({len(up_genes)} total)

TOP 10 DOWN GENES (paste all from atrt_down_genes.txt):
{top10_down}
  ... ({len(down_genes)} total)

══════════════════════════════════════════════════════
STEP-BY-STEP: CLUE.IO WEB SUBMISSION (recommended)
══════════════════════════════════════════════════════
1. Go to: https://clue.io
2. Sign up / log in (free academic account)
3. Click: Tools → Query → L1000 Query
4. "Up-regulated genes" box: paste contents of atrt_up_genes.txt
5. "Down-regulated genes" box: paste contents of atrt_down_genes.txt
6. Query name: ATRT_SMARCB1_loss_signature
7. Perturbation type: Compounds (trt_cp)
8. Cell lines: all
9. Click SUBMIT — results take ~10-15 minutes
10. When done: Download → query_result.gct
11. Move file:
    cp ~/Downloads/query_result.gct data/cmap_query/query_result.gct
    OR rename and place at: data/cmap_query/query_result.gct
12. Run integration:
    python -m backend.pipeline.integrate_cmap_results

══════════════════════════════════════════════════════
STEP-BY-STEP: API SUBMISSION
══════════════════════════════════════════════════════
python scripts/01_prepare_cmap.py --submit --api-key YOUR_CLUE_API_KEY

Get your API key from: https://clue.io → Account → API key

══════════════════════════════════════════════════════
L1000 LIBRARY GAPS — these drugs return neutral prior (0.50)
══════════════════════════════════════════════════════
{gap_notes}

══════════════════════════════════════════════════════
EXPECTED STRONG REVERSERS (norm_cs < -0.9 on clue.io)
══════════════════════════════════════════════════════
These drugs should reverse the ATRT SMARCB1-loss signature:
  PANOBINOSTAT  (HDAC1/2) — reverses HDAC upregulation
  BIRABRESIB    (BRD4)    — reverses BRD4/MYC super-enhancers
  ABEMACICLIB   (CDK4/6)  — reverses CDK4/CCND1 upregulation
  ONC201        (DRD2)    — profiled as TIC-10 in L1000
  VORINOSTAT    (HDAC)    — pan-HDAC reversal

While waiting for CMap results, pipeline runs with neutral prior (0.50)
for all drugs. Scores will improve significantly with real CMap data.
"""
    (out_dir / "atrt_query_manifest.txt").write_text(manifest)
    logger.info("Saved manifest: %s/atrt_query_manifest.txt", out_dir)


def submit_via_api(
    up_genes:   List[str],
    down_genes: List[str],
    api_key:    str,
    query_name: str = "ATRT_SMARCB1_loss_signature",
) -> Optional[str]:
    """Submit to clue.io API and poll for result."""
    import requests
    import time

    headers = {
        "x-user-key":   api_key,
        "Content-Type": "application/json",
        "Accept":       "application/json",
    }
    payload = {
        "name":    query_name,
        "tool_id": "sig_queryl1k_tool",
        "parameters": {
            "pert_type":   "trt_cp",
            "cell_id":     "all",
            "up_genes":    " ".join(up_genes[:500]),
            "dn_genes":    " ".join(down_genes[:500]),
            "ignore_dose": "1",
        },
    }

    logger.info("Submitting to clue.io API ...")
    try:
        resp = requests.post(
            "https://api.clue.io/api/jobs",
            headers=headers,
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()
    except requests.HTTPError as e:
        if e.response.status_code == 401:
            logger.error("❌ 401 Unauthorized — check your API key")
        else:
            logger.error("❌ HTTP %d: %s", e.response.status_code, e.response.text[:200])
        return None
    except Exception as e:
        logger.error("❌ API error: %s", e)
        return None

    job_id = resp.json().get("id")
    logger.info("✅ Job submitted: %s", job_id)
    logger.info("Polling for completion (max 20 min) ...")

    for _ in range(40):  # 40 × 30s = 20 min
        time.sleep(30)
        try:
            status_resp = requests.get(
                f"https://api.clue.io/api/jobs/{job_id}",
                headers=headers,
                timeout=15,
            )
            status = status_resp.json().get("status", "unknown")
            logger.info("  Status: %s", status)
            if status == "completed":
                # Try to download result
                result_url = status_resp.json().get("result_url", "")
                if result_url:
                    dl = requests.get(result_url, headers=headers, stream=True, timeout=60)
                    out_path = OUT_DIR / "query_result.gct"
                    out_path.parent.mkdir(parents=True, exist_ok=True)
                    with open(str(out_path), "wb") as f:
                        for chunk in dl.iter_content(8192):
                            f.write(chunk)
                    logger.info("✅ Result downloaded: %s", out_path)
                    logger.info("Run: python -m backend.pipeline.integrate_cmap_results")
                    return str(out_path)
                else:
                    logger.info("Job complete — download manually from clue.io dashboard")
                    return None
            elif status in ("failed", "error"):
                logger.error("❌ Job failed: %s", status_resp.json())
                return None
        except Exception as e:
            logger.warning("Poll error: %s", e)

    logger.warning("⏱️  Timed out — check clue.io dashboard manually")
    return None


def main():
    parser = argparse.ArgumentParser(description="ATRT CMap Query Preparation")
    parser.add_argument("--n-up",    default=150, type=int, help="Number of UP genes")
    parser.add_argument("--n-down",  default=50,  type=int, help="Number of DOWN genes")
    parser.add_argument("--submit",  action="store_true",   help="Submit to clue.io API")
    parser.add_argument("--api-key", default=None,          help="clue.io API key")
    args = parser.parse_args()

    logger.info("=" * 65)
    logger.info("ATRT CMap Query Preparation")
    logger.info("Computing ATRT vs normal brain differential expression")
    logger.info("=" * 65)

    # Validate inputs
    missing = []
    for path, label in [(ATRT_EXPR, "GSE70678"), (GTEX_REF, "GTEx brain")]:
        if not path.exists():
            missing.append((path, label))

    if missing:
        for path, label in missing:
            logger.error("Missing %s: %s", label, path)
        logger.error(
            "\nYour data files are at:\n"
            "  data/raw_omics/GSE70678_gene_expression.tsv\n"
            "  data/raw_omics/GTEx_brain_normal_reference.tsv\n"
            "These should already exist based on the image — check paths."
        )
        sys.exit(1)

    # Run DE
    atrt_matrix = load_atrt_matrix(ATRT_EXPR)
    gtex_ref    = load_gtex_reference(GTEX_REF)
    de_df       = run_welch_de(atrt_matrix, gtex_ref)

    if de_df.empty:
        logger.error("DE computation failed — check input data.")
        sys.exit(1)

    # Select genes and write outputs
    up_genes, down_genes = select_query_genes(de_df, n_up=args.n_up, n_down=args.n_down)
    write_outputs(up_genes, down_genes, de_df, OUT_DIR)

    # Print summary
    print("\n" + "=" * 65)
    print("✅ DONE — CMap query files generated")
    print("=" * 65)
    print(f"\nTop 15 UP genes in ATRT vs normal brain:")
    top15 = de_df[de_df["direction"] == "UP"].sort_values("zscore_log2FC", ascending=False).head(15)
    for _, row in top15.iterrows():
        print(f"  {row['gene']:<12}  log2FC={row['log2FC']:+.2f}  z={row['zscore_log2FC']:+.2f}  FDR={row['p_adj']:.3f}")

    print(f"\nTop 10 DOWN genes (expected: SMARCB1, CDKN2A, etc.):")
    top10d = de_df[de_df["direction"] == "DOWN"].sort_values("zscore_log2FC", ascending=True).head(10)
    for _, row in top10d.iterrows():
        print(f"  {row['gene']:<12}  log2FC={row['log2FC']:+.2f}  z={row['zscore_log2FC']:+.2f}  FDR={row['p_adj']:.3f}")

    print(f"\nFiles written to: {OUT_DIR}/")
    print(f"  atrt_up_genes.txt   ({len(up_genes)} genes)")
    print(f"  atrt_down_genes.txt ({len(down_genes)} genes)")
    print(f"  atrt_de_results.tsv (full table)")
    print(f"  atrt_query_manifest.txt (instructions)")
    print("\nNext step: see atrt_query_manifest.txt for clue.io submission")
    print("=" * 65)

    # API submission if requested
    if args.submit:
        if not args.api_key:
            logger.error("--api-key required for --submit. Get it from clue.io → Account.")
            sys.exit(1)
        submit_via_api(up_genes, down_genes, args.api_key)


if __name__ == "__main__":
    main()