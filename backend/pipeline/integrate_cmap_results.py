import json, logging, numpy as np, pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

CMAP_DIR    = Path("data/cmap_query")
QUERY_GCT   = CMAP_DIR / "query_result.gct"
OUTPUT_FILE = CMAP_DIR / "cmap_scores_pipeline.json"
SENTINEL    = -666


def norm_cs_to_pipeline_score(norm_cs: float) -> float:
    """
    Convert norm_cs [-2, +2] to pipeline cmap_score [0, 1].
    norm_cs = -2.0 → score = 1.00  (perfect reversal)
    norm_cs =  0.0 → score = 0.50  (neutral)
    norm_cs = +2.0 → score = 0.00  (mimics disease)
    Clipped to [0.05, 0.95].
    """
    score = (2.0 - norm_cs) / 4.0
    return round(max(0.05, min(0.95, score)), 4)


def is_reverser(norm_cs: float) -> bool:
    """
    Strong reversal threshold for norm_cs scale.
    Equivalent to tau < -90 in the old scale.
    norm_cs < -0.9 = strong reversal (top ~5% most negative).
    """
    return norm_cs < -0.9


def parse_gct13(path: Path) -> pd.DataFrame:
    """
    Parse GCT 1.3 format.

    Structure:
        Line 1: #1.3
        Line 2: nDataRows  nDataCols  nRowMeta  nColMeta
        Line 3: column headers (id + row_meta_cols + data_cols)
        Line 4: col metadata descriptor row
        Lines 5 to 4+nColMeta: column metadata rows
        Lines 5+nColMeta onwards: data rows
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
    logger.info(f"  Dims: {n_data_rows} rows × {n_data_cols} cols | "
                f"{n_row_meta} row-meta cols | {n_col_meta} col-meta rows")

    # Line 3: all column headers
    headers = lines[2].strip().split("\t")
    logger.info(f"  Headers: {headers}")

    # Skip n_col_meta metadata rows (lines 4 to 3+n_col_meta)
    data_start = 3 + n_col_meta

    # Read data rows
    data_lines = []
    for line in lines[data_start:]:
        stripped = line.strip()
        if stripped:
            data_lines.append(stripped.split("\t"))

    df = pd.DataFrame(data_lines, columns=headers)
    logger.info(f"  Loaded {len(df)} rows")
    return df


def process(df: pd.DataFrame) -> dict:
    """Extract norm_cs per drug and convert to pipeline scores."""

    # Find key columns
    name_col   = next((c for c in ["pert_iname","pert_id","id"] if c in df.columns), None)
    normcs_col = next((c for c in ["norm_cs","normcs","NCS"] if c in df.columns), None)
    rawcs_col  = next((c for c in ["raw_cs","rawcs"] if c in df.columns), None)
    fdr_col    = next((c for c in ["fdr_q_nlog10","fdr_q","fdr"] if c in df.columns), None)
    moa_col    = next((c for c in df.columns if c.lower() == "moa"), None)

    logger.info(f"  Name col: {name_col}")
    logger.info(f"  norm_cs col: {normcs_col}")
    logger.info(f"  raw_cs col: {rawcs_col}")
    logger.info(f"  fdr col: {fdr_col}")

    if not name_col:
        raise ValueError(f"No drug name column. Available: {list(df.columns)}")
    if not normcs_col:
        # Fall back to raw_cs
        if rawcs_col:
            logger.warning("norm_cs not found — using raw_cs instead")
            normcs_col = rawcs_col
        else:
            raise ValueError(f"No connectivity score column. Available: {list(df.columns)}")

    scores = {}
    n_sentinel = 0

    for _, row in df.iterrows():
        name = str(row[name_col]).upper().strip()
        if not name or name in ("NA", "NAN", ""):
            continue

        try:
            cs = float(row[normcs_col])
        except (ValueError, TypeError):
            continue

        # Filter sentinel values
        if abs(cs - SENTINEL) < 1 or np.isnan(cs):
            n_sentinel += 1
            continue

        entry = {
            "cmap_score":  norm_cs_to_pipeline_score(cs),
            "norm_cs":     round(cs, 4),
            "is_reverser": is_reverser(cs),
        }

        # Add raw_cs if available
        if rawcs_col:
            try:
                entry["raw_cs"] = round(float(row[rawcs_col]), 4)
            except:
                pass

        # Add FDR if available
        if fdr_col:
            try:
                fdr_val = float(row[fdr_col])
                if not np.isnan(fdr_val) and abs(fdr_val - SENTINEL) > 1:
                    entry["fdr_q_nlog10"] = round(fdr_val, 4)
            except:
                pass

        # Add MoA
        if moa_col:
            m = str(row[moa_col])
            if m not in ("-666", "nan", "NA", ""):
                entry["moa"] = m[:80]

        # Keep most negative norm_cs per drug name
        if name not in scores or entry["norm_cs"] < scores[name]["norm_cs"]:
            scores[name] = entry

    cs_vals = [d["norm_cs"] for d in scores.values()]
    logger.info(f"  Sentinel filtered: {n_sentinel}")
    logger.info(f"  Unique drugs: {len(scores)}")
    if cs_vals:
        logger.info(f"  norm_cs range: [{min(cs_vals):.3f}, {max(cs_vals):.3f}]")
    return scores


def print_report(scores: dict):
    top20 = sorted(scores.items(), key=lambda x: x[1]["norm_cs"])[:20]
    print(f"\n{'='*70}")
    print("Top 20 H3K27M REVERSERS (most negative norm_cs)")
    print("norm_cs < -0.9 = strong reversal  |  norm_cs > +0.9 = mimics disease")
    print(f"{'='*70}")
    print(f"{'Drug':<28} {'norm_cs':>8}  {'Score':>6}  {'Reverser?'}")
    print("-"*65)
    for name, d in top20:
        flag = "YES ✓" if d["is_reverser"] else "—"
        moa  = f"  [{d['moa'][:28]}]" if "moa" in d else ""
        fdr  = f"  fdr=-log10({d['fdr_q_nlog10']:.2f})" if "fdr_q_nlog10" in d else ""
        print(f"{name:<28} {d['norm_cs']:>8.3f}  {d['cmap_score']:>6.3f}  {flag}{moa}{fdr}")

    candidates = [
        "BIRABRESIB","OTX015","PANOBINOSTAT","MARIZOMIB",
        "ABEMACICLIB","ONC201","ONATASERTIB","AZD-8055",
        "PAXALISIB","VORINOSTAT","TAZEMETOSTAT","INDOXIMOD",
    ]
    print(f"\n{'='*70}")
    print("Your pipeline candidates:")
    print(f"{'='*70}")
    for drug in candidates:
        match = drug
        if drug not in scores:
            # Try partial
            close = [k for k in scores if drug[:6] in k or k[:6] in drug]
            if close:
                match = close[0]
        if match in scores:
            d   = scores[match]
            sig = ("✓ REVERSES H3K27M" if d["is_reverser"] else
                   ("~ partial"         if d["norm_cs"] < -0.5 else
                    ("✗ mimics disease" if d["norm_cs"] > 0.5  else "— neutral")))
            moa = f"  [{d['moa'][:25]}]" if "moa" in d else ""
            label = f"{drug}" + (f" (as {match})" if match != drug else "")
            print(f"  {label:<30} norm_cs={d['norm_cs']:>7.3f}  {sig}{moa}")
        else:
            print(f"  {drug:<30} not in CMap database")

    n_rev = sum(1 for d in scores.values() if d["is_reverser"])
    print(f"\n  Strong reversers (norm_cs < -0.9): {n_rev} / {len(scores)} drugs")
    print(f"  Interpretation: these drugs reverse the H3K27M DIPG gene expression program")


def main():
    logger.info("="*55)
    logger.info("Integrating clue.io query_result.gct (v6)")
    logger.info("="*55)

    if not QUERY_GCT.exists():
        raise FileNotFoundError(
            f"Missing: {QUERY_GCT}\n"
            "Run:\n"
            "  cp ~/Downloads/my_analysis.sig_queryl1k_tool.*"
            "/arfs/TAG/query_result.gct data/cmap_query/query_result.gct"
        )

    df     = parse_gct13(QUERY_GCT)
    scores = process(df)

    CMAP_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        json.dump(scores, f, indent=2)
    logger.info(f"\n✅ Saved {len(scores)} drug scores → {OUTPUT_FILE}")

    print_report(scores)
    print(f"\n✅ Done. Re-run pipeline:")
    print(f"   python -m backend.pipeline.save_results")

if __name__ == "__main__":
    main()