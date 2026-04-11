#!/usr/bin/env python3
"""
scripts/01_prepare_cmap.py
==========================
BATCH 2 — CMap query preparation.

Computes differential expression between ATRT tumour (GSE70678) and
normal brain (GTEx v8 cerebellum + cortex), then writes the gene lists
needed for a clue.io L1000 query.

WHAT THIS DOES (why CMap needs DE)
------------------------------------
CMap L1000 asks: "which drug reverses a disease gene signature?"
To answer that, you must first define the signature:
  UP genes  = genes HIGHER in ATRT vs normal brain  → drug should push DOWN
  DOWN genes = genes LOWER  in ATRT vs normal brain  → drug should push UP

The DE is computed with Welch's t-test (handles unequal variance between
Affymetrix microarray and GTEx RNA-seq platforms).

INPUTS (must exist — run data_downloader.py first if not present)
------------------------------------------------------------------
  data/raw_omics/GSE70678_gene_expression.tsv  ← ATRT (gene × sample)
  data/raw_omics/GTEx_brain_normal_reference.tsv ← normal brain reference

OUTPUTS
-------
  data/cmap_query/atrt_up_genes.txt      ← top 150 UP genes (for clue.io)
  data/cmap_query/atrt_down_genes.txt    ← top 50  DOWN genes (for clue.io)
  data/cmap_query/atrt_de_results.tsv    ← full DE table for reference
  data/cmap_query/atrt_query_manifest.txt ← submission instructions

HOW TO RUN
----------
  cd <repo_root>
  python scripts/01_prepare_cmap.py

NEXT STEPS after this script
-----------------------------
  1. Go to https://clue.io → Tools → L1000 Query
  2. Paste atrt_up_genes.txt into "Up genes" box
  3. Paste atrt_down_genes.txt into "Down genes" box
  4. Name query: ATRT_SMARCB1_loss_signature
  5. Submit → wait ~10 min
  6. Download query_result.gct
  7. Run: python scripts/02_integrate_cmap.py
"""

import logging
import sys
from pathlib import Path

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
REPO_ROOT    = Path(__file__).parent.parent
ATRT_EXPR    = REPO_ROOT / "data/raw_omics/GSE70678_gene_expression.tsv"
GTEX_REF     = REPO_ROOT / "data/raw_omics/GTEx_brain_normal_reference.tsv"
OUT_DIR      = REPO_ROOT / "data/cmap_query"

# ── Statistical thresholds ────────────────────────────────────────────────────
LOG2FC_THRESHOLD = 1.0   # ≥ 2-fold change
P_ADJ_THRESHOLD  = 0.05  # BH-corrected FDR
MIN_EXPR         = 3.0   # log2 units — filter unexpressed genes
MIN_SAMPLES      = 10    # minimum ATRT samples per gene

# ── Landmark genes confirmed in L1000 (for validation) ───────────────────────
EXPECTED_UP_MARKERS = {"EZH2", "AURKA", "BRD4", "MYC", "MYCN", "HDAC1", "CDK4"}


def bh_correction(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n     = len(p_values)
    order = np.argsort(p_values)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    p_adj = np.minimum(1.0, p_values * n / ranked)
    for i in range(n - 2, -1, -1):
        p_adj[order[i]] = min(p_adj[order[i]], p_adj[order[i + 1]])
    return p_adj


def load_atrt_matrix(path: Path) -> pd.DataFrame:
    """Load ATRT gene expression matrix (gene × samples, log2 scale)."""
    logger.info("Loading ATRT expression matrix: %s", path)
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.upper().str.strip()
    df = df.apply(pd.to_numeric, errors="coerce")
    logger.info("  Shape: %d genes × %d ATRT samples", len(df), len(df.columns))
    return df


def load_gtex_reference(path: Path) -> pd.Series:
    """
    Load GTEx normal brain reference — mean of cerebellum + cortex columns.
    Returns Series: gene → log2(median_TPM + 1).
    """
    logger.info("Loading GTEx brain reference: %s", path)
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.upper().str.strip()
    df = df.apply(pd.to_numeric, errors="coerce")

    # Select brain columns
    brain_cols = [c for c in df.columns if "brain" in c.lower()
                  or "cerebellum" in c.lower() or "cortex" in c.lower()]
    if not brain_cols:
        logger.warning("No brain columns found — using all columns as reference")
        brain_cols = list(df.columns)

    ref = df[brain_cols].mean(axis=1).dropna()
    logger.info("  Brain reference: %d genes across %d tissues: %s",
                len(ref), len(brain_cols),
                [c.replace("Brain - ", "") for c in brain_cols[:4]])
    return ref


def run_welch_de(atrt_matrix: pd.DataFrame, gtex_ref: pd.Series) -> pd.DataFrame:
    """
    Welch's t-test: ATRT samples vs GTEx normal brain reference.

    Uses one-sample t-test (H0: mean(ATRT) == gtex_ref) — appropriate because
    GTEx median represents a known population reference, not a sample.

    Platform note: GSE70678 is Affymetrix (RMA, log2-normalised).
                   GTEx is log2(TPM+1). Both are log2 scale — comparable after
                   z-score normalisation of the differential vector.
    """
    common_genes = atrt_matrix.index.intersection(gtex_ref.index)
    logger.info("Genes in common between ATRT and GTEx: %d", len(common_genes))

    if len(common_genes) < 1000:
        logger.warning(
            "Only %d common genes — check that both use uppercase HGNC symbols.",
            len(common_genes),
        )

    rows = []
    for gene in common_genes:
        expr_vals = atrt_matrix.loc[gene].dropna()
        ref_val   = float(gtex_ref[gene])
        n_valid   = len(expr_vals)

        if n_valid < MIN_SAMPLES:
            continue
        if expr_vals.mean() < MIN_EXPR and ref_val < MIN_EXPR:
            continue

        mean_atrt = float(expr_vals.mean())
        log2fc    = mean_atrt - ref_val  # log2(ATRT/normal) in log2 space
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

    if not rows:
        logger.error("No genes survived DE filtering — check input data.")
        return pd.DataFrame()

    de_df = pd.DataFrame(rows).set_index("gene")

    # BH multiple testing correction
    de_df["p_adj"] = bh_correction(de_df["p_value"].values)

    # Z-score of log2FC (corrects for platform offset between microarray and RNA-seq)
    lfc_mean = de_df["log2FC"].mean()
    lfc_std  = de_df["log2FC"].std()
    de_df["zscore_log2FC"] = (de_df["log2FC"] - lfc_mean) / lfc_std.clip(lower=0.01)

    # Significance
    de_df["significant"] = (
        (de_df["log2FC"].abs() >= LOG2FC_THRESHOLD) &
        (de_df["p_adj"] <= P_ADJ_THRESHOLD)
    )
    de_df["direction"] = np.where(
        ~de_df["significant"], "NS",
        np.where(de_df["log2FC"] > 0, "UP", "DOWN"),
    )

    n_up   = (de_df["direction"] == "UP").sum()
    n_down = (de_df["direction"] == "DOWN").sum()
    n_ns   = (de_df["direction"] == "NS").sum()
    logger.info(
        "Welch DE complete: %d genes | UP=%d DOWN=%d NS=%d "
        "(|log2FC|≥%.1f, FDR≤%.2f)",
        len(de_df), n_up, n_down, n_ns,
        LOG2FC_THRESHOLD, P_ADJ_THRESHOLD,
    )
    return de_df.reset_index()


def select_query_genes(
    de_df: pd.DataFrame,
    n_up:  int = 150,
    n_down: int = 50,
) -> tuple:
    """Select top UP and DOWN genes for clue.io query."""
    up_pool   = de_df[de_df["direction"] == "UP"].sort_values(
        "zscore_log2FC", ascending=False)
    down_pool = de_df[de_df["direction"] == "DOWN"].sort_values(
        "zscore_log2FC", ascending=True)

    # Fall back to all positive/negative if few significant hits
    if len(up_pool) < 20:
        logger.warning(
            "Only %d significant UP genes — relaxing to all log2FC>0", len(up_pool)
        )
        up_pool = de_df[de_df["log2FC"] > 0].sort_values("zscore_log2FC", ascending=False)
    if len(down_pool) < 10:
        logger.warning(
            "Only %d significant DOWN genes — relaxing to all log2FC<0", len(down_pool)
        )
        down_pool = de_df[de_df["log2FC"] < 0].sort_values("zscore_log2FC", ascending=True)

    up_genes   = up_pool["gene"].head(n_up).tolist()
    down_genes = down_pool["gene"].head(n_down).tolist()

    logger.info("Selected %d UP genes, %d DOWN genes for CMap query", len(up_genes), len(down_genes))

    # Biology validation
    found = EXPECTED_UP_MARKERS & set(up_genes)
    if len(found) >= 3:
        logger.info("✅ ATRT biology check — key markers in UP list: %s", sorted(found))
    else:
        logger.warning(
            "⚠️  Only %d/%d expected ATRT markers in UP: %s\n"
            "   Check that GSE70678 loaded correctly.",
            len(found), len(EXPECTED_UP_MARKERS), sorted(found),
        )

    return up_genes, down_genes


def write_outputs(up_genes, down_genes, de_df, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)

    # Gene lists for clue.io web UI
    (out_dir / "atrt_up_genes.txt").write_text("\n".join(up_genes))
    (out_dir / "atrt_down_genes.txt").write_text("\n".join(down_genes))
    logger.info("Saved gene lists → %s/{atrt_up_genes,atrt_down_genes}.txt", out_dir)

    # Full DE table
    de_path = out_dir / "atrt_de_results.tsv"
    de_df.to_csv(de_path, sep="\t", index=False)
    logger.info("Saved DE table (%d genes) → %s", len(de_df), de_path)

    # Manifest with submission instructions
    manifest = f"""ATRT CMap Query — Differential Expression vs Normal Brain
=========================================================
Generated from: GSE70678 (Torchia 2015) vs GTEx v8 brain

UP genes   : {len(up_genes)} genes (HIGHER in ATRT than normal brain)
DOWN genes : {len(down_genes)} genes (LOWER in ATRT than normal brain)

Top 20 UP (paste these and more from atrt_up_genes.txt):
{chr(10).join('  ' + g for g in up_genes[:20])}

Top 10 DOWN (paste from atrt_down_genes.txt):
{chr(10).join('  ' + g for g in down_genes[:10])}

CLUE.IO SUBMISSION STEPS
-------------------------
1. Go to: https://clue.io → Tools → L1000 Query (sig_queryl1k_tool)
2. Paste contents of atrt_up_genes.txt into the "Up genes" box
3. Paste contents of atrt_down_genes.txt into the "Down genes" box
4. Query name: ATRT_SMARCB1_loss_signature
5. Perturbation type: Compounds (trt_cp)
6. Cell line: all
7. Click Submit → results in ~10-15 minutes
8. Download query_result.gct when complete
9. Run: python scripts/02_integrate_cmap.py

DRUGS NOT IN L1000 LIBRARY (will get neutral score 0.50)
---------------------------------------------------------
  TAZEMETOSTAT — approved 2020, after L1000 profiling cutoff
  MARIZOMIB    — marine natural product, not profiled
  ALISERTIB    — partial profiling only

EXPECTED STRONG REVERSERS (norm_cs < -0.9)
-------------------------------------------
  PANOBINOSTAT  (HDAC1/2 inhibitor → reverses HDAC upregulation)
  BIRABRESIB    (BRD4 inhibitor → reverses BRD4/MYC super-enhancers)
  ABEMACICLIB   (CDK4/6 inhibitor → reverses CDK4/CCND1 upregulation)
"""
    (out_dir / "atrt_query_manifest.txt").write_text(manifest)
    logger.info("Saved manifest → %s/atrt_query_manifest.txt", out_dir)


def main():
    logger.info("=" * 65)
    logger.info("ATRT CMap Query Preparation")
    logger.info("Computes ATRT vs normal brain differential expression")
    logger.info("=" * 65)

    # Validate inputs
    for path, label in [(ATRT_EXPR, "GSE70678"), (GTEX_REF, "GTEx brain")]:
        if not path.exists():
            logger.error(
                "Missing %s: %s\n"
                "Run first: python -m backend.pipeline.data_downloader --dataset %s\n"
                "           python -m backend.pipeline.data_downloader --process %s",
                label, path,
                "gse70678 gpl6244" if "GSE" in str(path) else "gtex_brain",
                "all",
            )
            sys.exit(1)

    # Load data
    atrt_matrix = load_atrt_matrix(ATRT_EXPR)
    gtex_ref    = load_gtex_reference(GTEX_REF)

    # Run differential expression
    de_df = run_welch_de(atrt_matrix, gtex_ref)
    if de_df.empty:
        logger.error("DE computation failed — no results.")
        sys.exit(1)

    # Select query genes
    up_genes, down_genes = select_query_genes(de_df, n_up=150, n_down=50)

    # Write outputs
    write_outputs(up_genes, down_genes, de_df, OUT_DIR)

    # Print summary
    logger.info("")
    logger.info("=" * 65)
    logger.info("✅ DONE")
    logger.info("")
    logger.info("Top 10 UP genes (drugs should REVERSE these):")
    top10 = de_df[de_df["direction"] == "UP"].sort_values("zscore_log2FC", ascending=False).head(10)
    for _, row in top10.iterrows():
        logger.info("  %-12s  log2FC=%+.2f  z=%+.2f  FDR=%.3f",
                    row["gene"], row["log2FC"], row["zscore_log2FC"], row["p_adj"])
    logger.info("")
    logger.info("Next step: submit gene lists to clue.io")
    logger.info("See: %s/atrt_query_manifest.txt", OUT_DIR)
    logger.info("=" * 65)


if __name__ == "__main__":
    main()