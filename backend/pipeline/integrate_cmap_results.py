"""
integrate_cmap_results.py  (ATRT v2.0)
========================================
Parses clue.io query_result.gct and converts connectivity scores
to pipeline-ready JSON for the ATRT SMARCB1-loss signature.

BIOLOGICAL CONTEXT
-------------------
The ATRT SMARCB1-loss signature (GSE70678 vs GTEx normal brain) identifies
drugs that REVERSE the transcriptional consequences of SMARCB1 biallelic loss:
  - EZH2 hyperactivation (PRC2 unchecked)
  - BRD4/super-enhancer dependency
  - HDAC1/2 compensatory activity
  - AURKA-mediated MYCN stabilisation

A drug with norm_cs < -0.9 reverses this program → candidate for synthetic
lethality in SMARCB1-null ATRT.

Source: Subramanian A et al. Cell 2017; 171(6):1437. PMID 29195078.
        Knutson SK et al. PNAS 2013; 110(19):7922. PMID 23620515.
        Torchia J et al. Cancer Cell 2015; 30(6):891. PMID 26609405.

L1000 LIBRARY GAPS (these drugs return neutral prior 0.50):
  TAZEMETOSTAT : Approved 2020 — post L1000 profiling cutoff
  MARIZOMIB    : Marine natural product — not commercially profiled
  ALISERTIB    : Partial profiling only — AURKA not a landmark gene

USAGE
------
  cp ~/Downloads/query_result.gct data/cmap_query/query_result.gct
  python -m backend.pipeline.integrate_cmap_results

OUTPUT
------
  data/cmap_query/atrt_cmap_scores.json
"""

import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

# ── Paths ─────────────────────────────────────────────────────────────────────
CMAP_DIR    = Path("data/cmap_query")
QUERY_GCT   = CMAP_DIR / "query_result.gct"
OUTPUT_FILE = CMAP_DIR / "atrt_cmap_scores.json"   # ← ATRT-specific output
SENTINEL    = -666

# ── ATRT candidate drugs for validation report ────────────────────────────────
# Source: Torchia 2015, Knutson 2013, Geoerger 2017, Sredni 2017
ATRT_PIPELINE_CANDIDATES = [
    # EZH2 synthetic lethality (Knutson 2013 PNAS PMID 23620515)
    "TAZEMETOSTAT",
    # Pan-HDAC (Torchia 2015 Cancer Cell PMID 26609405)
    "PANOBINOSTAT", "VORINOSTAT",
    # AURKA (Sredni 2017 Pediatric Blood Cancer PMID 28544500)
    "ALISERTIB",
    # BET bromodomain (Geoerger 2017 Clin Cancer Res PMID 28108534)
    "BIRABRESIB", "OTX015",
    # CDK4/6
    "ABEMACICLIB", "PALBOCICLIB", "RIBOCICLIB",
    # Proteasome (Lin 2019 Sci Transl Med)
    "MARIZOMIB", "BORTEZOMIB",
    # SHH subgroup (Johann 2016 Cancer Cell PMID 26923874)
    "VISMODEGIB", "SONIDEGIB", "ITRACONAZOLE",
    # DRD2/TRAIL (ONC201 — Venneti 2023)
    "ONC201", "DORDAVIPRONE",
    # PI3K/mTOR
    "PAXALISIB", "GDC-0941", "AZD-8055",
    # IDO pathway
    "INDOXIMOD",
    # Generic repurposables
    "VALPROIC ACID", "SIROLIMUS", "METFORMIN",
    "CHLOROQUINE", "HYDROXYCHLOROQUINE", "ARSENIC TRIOXIDE",
]

# Drugs confirmed NOT in L1000 library
# Source: clue.io LINCS compound library manifest (accessed April 2026)
L1000_NOT_IN_LIBRARY = {
    "TAZEMETOSTAT": "Post-2017 compound (FDA approved 2020); absent from L1000 library",
    "MARIZOMIB":    "Marine natural product; not commercially profiled in L1000",
}

# Drugs with partial/proxy profiling in L1000
L1000_PARTIAL_OR_PROXY = {
    "ALISERTIB":   "Partial profiling; AURKA not a standard landmark gene",
    "BIRABRESIB":  "Profiled as OTX-015 in L1000",
    "ONC201":      "Profiled as TIC-10 in L1000",
    "ABEMACICLIB": "CDK4 shRNA knockdown used as functional proxy",
    "PAXALISIB":   "GDC-0941 (PI3K class) used as proxy",
}


# ── Score conversion ──────────────────────────────────────────────────────────

def norm_cs_to_pipeline_score(norm_cs: float) -> float:
    """
    Convert norm_cs [-2, +2] to pipeline cmap_score [0, 1].

    norm_cs = -2.0 → score = 1.00  (perfect reversal of ATRT signature)
    norm_cs =  0.0 → score = 0.50  (neutral — no relationship)
    norm_cs = +2.0 → score = 0.00  (mimics ATRT SMARCB1-loss program)

    Clipped to [0.05, 0.95] to avoid overconfident extremes.
    Source: Subramanian 2017 Cell PMID 29195078.
    """
    score = (2.0 - norm_cs) / 4.0
    return round(max(0.05, min(0.95, score)), 4)


def is_atrt_reverser(norm_cs: float) -> bool:
    """
    Strong reversal threshold: norm_cs < -0.9.
    Equivalent to top ~5% most negative across all 38,973 drug profiles.
    These drugs most strongly reverse the ATRT SMARCB1-loss transcriptional program.
    Source: Subramanian 2017 Cell — norm_cs interpretation.
    """
    return norm_cs < -0.9


# ── GCT 1.3 parser ────────────────────────────────────────────────────────────

def parse_gct13(path: Path) -> pd.DataFrame:
    """
    Parse GCT 1.3 format from clue.io query result.

    Structure:
        Line 1: #1.3
        Line 2: nDataRows  nDataCols  nRowMeta  nColMeta
        Line 3: column headers (id + row_meta_cols + data_cols)
        Lines 4 to 3+nColMeta: column metadata rows
        Lines 4+nColMeta onwards: data rows
    """
    logger.info(f"Parsing GCT 1.3: {path}")
    with open(path) as f:
        lines = f.readlines()

    version = lines[0].strip()
    logger.info(f"  Version: {version}")

    dims        = lines[1].strip().split("\t")
    n_data_rows = int(dims[0])
    n_data_cols = int(dims[1])
    n_row_meta  = int(dims[2])
    n_col_meta  = int(dims[3])
    logger.info(
        f"  Dims: {n_data_rows} rows × {n_data_cols} cols | "
        f"{n_row_meta} row-meta cols | {n_col_meta} col-meta rows"
    )

    headers    = lines[2].strip().split("\t")
    data_start = 3 + n_col_meta

    data_lines = []
    for line in lines[data_start:]:
        stripped = line.strip()
        if stripped:
            data_lines.append(stripped.split("\t"))

    df = pd.DataFrame(data_lines, columns=headers)
    logger.info(f"  Loaded {len(df)} rows")
    return df


# ── Score extraction ──────────────────────────────────────────────────────────

def process_gct(df: pd.DataFrame) -> dict:
    """
    Extract norm_cs per drug from parsed GCT DataFrame.
    Converts to ATRT pipeline scores and flags reversers.
    """
    # Find key columns (clue.io column names may vary slightly)
    name_col   = next((c for c in ["pert_iname", "pert_id", "id"]
                       if c in df.columns), None)
    normcs_col = next((c for c in ["norm_cs", "normcs", "NCS"]
                       if c in df.columns), None)
    rawcs_col  = next((c for c in ["raw_cs", "rawcs"]
                       if c in df.columns), None)
    fdr_col    = next((c for c in ["fdr_q_nlog10", "fdr_q", "fdr"]
                       if c in df.columns), None)
    moa_col    = next((c for c in df.columns if c.lower() == "moa"), None)

    logger.info(f"  Name col: {name_col}")
    logger.info(f"  norm_cs col: {normcs_col}")
    logger.info(f"  fdr col: {fdr_col}")

    if not name_col:
        raise ValueError(f"No drug name column found. Available: {list(df.columns)}")
    if not normcs_col:
        if rawcs_col:
            logger.warning("norm_cs not found — using raw_cs as fallback")
            normcs_col = rawcs_col
        else:
            raise ValueError(f"No connectivity score column. Available: {list(df.columns)}")

    scores       = {}
    n_sentinel   = 0
    n_no_lib     = 0  # drugs confirmed not in L1000 library

    for _, row in df.iterrows():
        name = str(row[name_col]).upper().strip()
        if not name or name in ("NA", "NAN", ""):
            continue

        try:
            cs = float(row[normcs_col])
        except (ValueError, TypeError):
            continue

        # Filter sentinel values (clue.io uses -666 as missing data flag)
        if abs(cs - SENTINEL) < 1 or np.isnan(cs):
            n_sentinel += 1
            continue

        entry = {
            "cmap_score":  norm_cs_to_pipeline_score(cs),
            "norm_cs":     round(cs, 4),
            "is_reverser": is_atrt_reverser(cs),
        }

        if rawcs_col and rawcs_col != normcs_col:
            try:
                entry["raw_cs"] = round(float(row[rawcs_col]), 4)
            except Exception:
                pass

        if fdr_col:
            try:
                fdr_val = float(row[fdr_col])
                if not np.isnan(fdr_val) and abs(fdr_val - SENTINEL) > 1:
                    entry["fdr_q_nlog10"] = round(fdr_val, 4)
            except Exception:
                pass

        if moa_col:
            m = str(row[moa_col])
            if m not in ("-666", "nan", "NA", ""):
                entry["moa"] = m[:80]

        # Keep most negative norm_cs per drug (strongest reverser across cell lines)
        if name not in scores or entry["norm_cs"] < scores[name]["norm_cs"]:
            scores[name] = entry

    cs_vals = [d["norm_cs"] for d in scores.values()]
    logger.info(f"  Sentinel values filtered: {n_sentinel}")
    logger.info(f"  Unique drugs with scores: {len(scores)}")
    if cs_vals:
        logger.info(
            f"  norm_cs range: [{min(cs_vals):.3f}, {max(cs_vals):.3f}]"
        )
    return scores


# ── ATRT-specific report ──────────────────────────────────────────────────────

def print_atrt_report(scores: dict) -> None:
    """
    Print validation report for ATRT pipeline candidates.

    Interpretation:
      norm_cs < -0.9 : Strong reversal of ATRT SMARCB1-loss program → priority candidate
      norm_cs < -0.5 : Partial reversal → moderate support
      norm_cs > +0.5 : Mimics ATRT program → avoid
      Not profiled   : Return neutral prior (0.50) — no information penalty

    Source: Knutson 2013 PNAS — EZH2 synthetic lethality
            Torchia 2015 Cancer Cell — ATRT transcriptional signature
    """
    # Top 20 strongest reversers of ATRT SMARCB1-loss signature
    valid_scores = {k: v for k, v in scores.items() if v.get("norm_cs") is not None}
    
    # Sort only the valid numeric scores
    top20 = sorted(valid_scores.items(), key=lambda x: x[1]["norm_cs"])[:20]

    print(f"\n{'='*75}")
    print("TOP 20 REVERSERS OF ATRT SMARCB1-LOSS SIGNATURE (most negative norm_cs)")
    print("Drugs that reverse: EZH2↑ BRD4↑ HDAC1↑ AURKA↑ MYC↑ vs SMARCB1↓")
    print(f"norm_cs < -0.9 = strong reversal (top ~5% most negative)")
    print(f"norm_cs > +0.9 = mimics ATRT program (avoid)")
    print(f"{'='*75}")
    print(f"{'Drug':<30} {'norm_cs':>8}  {'Score':>6}  {'Reverser?':<12}  MoA")
    print("-" * 75)
    for name, d in top20:
        flag = "YES ✓ STRONG" if d["is_reverser"] else (
            "partial" if d["norm_cs"] < -0.5 else "—"
        )
        moa  = f"  [{d['moa'][:25]}]" if "moa" in d else ""
        fdr  = (f"  fdr=-log10({d['fdr_q_nlog10']:.1f})"
                if "fdr_q_nlog10" in d else "")
        print(f"{name:<30} {d['norm_cs']:>8.3f}  {d['cmap_score']:>6.3f}  "
              f"{flag:<12}{moa}{fdr}")

    # ATRT candidate validation
    print(f"\n{'='*75}")
    print("ATRT PIPELINE CANDIDATES — CMap Validation")
    print("(Source: Torchia 2015, Knutson 2013, Geoerger 2017, Sredni 2017, Johann 2016)")
    print(f"{'='*75}")

    n_reverser = 0
    n_partial  = 0
    n_neutral  = 0
    n_not_profiled = 0

    for drug in ATRT_PIPELINE_CANDIDATES:
        upper = drug.upper()

        # Check L1000 library status first
        if upper in L1000_NOT_IN_LIBRARY:
            print(f"  {drug:<30} NOT PROFILED  score=0.50 (neutral prior)  "
                  f"[{L1000_NOT_IN_LIBRARY[upper]}]")
            n_not_profiled += 1
            continue

        # Try to find in scores (exact or partial match)
        match = upper
        found = False
        if upper in scores:
            found = True
        else:
            close = [k for k in scores if upper[:6] in k or k[:6] in upper]
            if close:
                match = close[0]
                found = True

        proxy_note = ""
        if upper in L1000_PARTIAL_OR_PROXY:
            proxy_note = f"  [{L1000_PARTIAL_OR_PROXY[upper]}]"

        if found:
            d   = scores[match]
            label = match if match == upper else f"{drug} (as {match})"
            if d["is_reverser"]:
                sig = "✓ REVERSES ATRT SMARCB1-loss program"
                n_reverser += 1
            elif d["norm_cs"] < -0.5:
                sig = "~ partial reversal"
                n_partial += 1
            elif d["norm_cs"] > 0.5:
                sig = "✗ MIMICS ATRT program — AVOID"
            else:
                sig = "— neutral"
                n_neutral += 1
            moa = f"  [{d['moa'][:20]}]" if "moa" in d else ""
            print(f"  {label:<34} norm_cs={d['norm_cs']:>7.3f}  {sig}{proxy_note}{moa}")
        else:
            print(f"  {drug:<34} not found in CMap database{proxy_note}")
            n_not_profiled += 1

    n_rev_total = sum(1 for d in scores.values() if d["is_reverser"])
    print(f"\n  Strong reversers (norm_cs < -0.9): {n_rev_total} / {len(scores)} drugs")
    print(f"  ATRT candidates — strong: {n_reverser} | partial: {n_partial} "
          f"| neutral: {n_neutral} | not profiled: {n_not_profiled}")

    print(f"\n  BIOLOGICAL INTERPRETATION:")
    print(f"  Strong reversers likely target SMARCB1-loss dependencies:")
    print(f"    EZH2 hyperactivation → HDAC/BET inhibitors reverse")
    print(f"    BRD4 super-enhancer dependency → BET inhibitors reverse")
    print(f"    AURKA/MYCN axis → Aurora kinase inhibitors")
    print(f"  Source: Knutson 2013 PNAS; Torchia 2015 Cancer Cell")

    print(f"\n  L1000 LIBRARY GAPS (neutral prior 0.50 applied):")
    for drug, reason in L1000_NOT_IN_LIBRARY.items():
        print(f"    {drug}: {reason}")
    print(f"  These drugs are still scored via DepMap + tissue expression only.")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    logger.info("=" * 65)
    logger.info("ATRT CMap Integration — SMARCB1-loss Signature Reversal")
    logger.info("Source: Torchia 2015 Cancer Cell; Knutson 2013 PNAS")
    logger.info("=" * 65)

    if not QUERY_GCT.exists():
        raise FileNotFoundError(
            f"\n❌ Missing: {QUERY_GCT}\n\n"
            "HOW TO GET THIS FILE:\n"
            "  1. Run: python scripts/01_prepare_cmap.py\n"
            "  2. Go to: https://clue.io → Tools → L1000 Query\n"
            "  3. Paste atrt_up_genes.txt into UP box\n"
            "  4. Paste atrt_down_genes.txt into DOWN box\n"
            "  5. Name: ATRT_SMARCB1_loss_signature | Perturbation type: Compounds\n"
            "  6. Submit → wait ~10-15 min → Download query_result.gct\n"
            "  7. cp ~/Downloads/query_result.gct data/cmap_query/query_result.gct\n"
            "  8. Re-run this script\n"
        )

    df     = parse_gct13(QUERY_GCT)
    scores = process_gct(df)

    if not scores:
        raise ValueError("No drug scores extracted from GCT file — check file format.")

    # Mark L1000 library status for downstream use
    for drug, reason in L1000_NOT_IN_LIBRARY.items():
        if drug not in scores:
            scores[drug] = {
                "cmap_score":  0.50,
                "norm_cs":     None,
                "is_reverser": False,
                "not_in_library": True,
                "library_note": reason,
            }

    for drug, note in L1000_PARTIAL_OR_PROXY.items():
        if drug in scores:
            scores[drug]["library_note"] = note
            scores[drug]["partial_profiling"] = True

    # Save to ATRT-specific output path
    CMAP_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        json.dump(scores, f, indent=2)

    logger.info(f"\n✅ Saved {len(scores)} ATRT drug scores → {OUTPUT_FILE}")
    logger.info(f"   Strong reversers (norm_cs < -0.9): "
                f"{sum(1 for d in scores.values() if d.get('is_reverser'))}")

    print_atrt_report(scores)

    print(f"\n{'='*65}")
    print(f"✅ CMap integration complete → {OUTPUT_FILE}")
    print(f"\nNext steps:")
    print(f"  python -m backend.pipeline.save_results --disease atrt")
    print(f"  python -m backend.pipeline.generate_figures")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()