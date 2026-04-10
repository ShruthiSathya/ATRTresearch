"""
prepare_cmap_query.py
=====================
Prepare ATRT gene signature for clue.io CMap L1000 query.

WORKFLOW
---------
1. Load GSE70678 ATRT bulk RNA-seq (Torchia 2015)
2. Compute differential expression: ATRT samples vs normal brain
3. Select top 50–150 up/down genes meeting threshold
4. Convert gene symbols to Entrez IDs via MyGene.info API
5. Save gene lists for clue.io query submission
6. Print instructions for clue.io submission

AFTER THIS SCRIPT
-----------------
1. Go to https://clue.io and create a free account
2. Click 'Query' → 'L1000 Query'
3. Paste contents of data/cmap_query/atrt_up_genes.txt into UP box
4. Paste contents of data/cmap_query/atrt_down_genes.txt into DOWN box
5. Set name: 'ATRT_SMARCB1_signature'
6. Click Submit (results ready in ~10 minutes)
7. Download results → save as data/cmap_query/atrt_query_result.gct
8. Run integrate_cmap_results.py --output atrt_cmap_scores.json

BIOLOGICAL RATIONALE FOR ATRT SIGNATURE
-----------------------------------------
The ATRT CMap signature captures genes upregulated due to SMARCB1 loss:
  - PRC2/EZH2 targets (H3K27me3 regulated genes)
  - MYC transcriptional program
  - BET bromodomain targets
  - AURKA/MYCN axis

Drugs with strong negative connectivity score (norm_cs < -0.9) reverse this
signature and are prioritised as candidates.

Source: Torchia 2015 Cancer Cell; Subramanian 2017 Cell (CMap L1000 methods)
"""

import logging
import requests
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

try:
    from .pipeline_config import PATHS, GENOMICS, CMAP
except ImportError:
    from pipeline_config import PATHS, GENOMICS, CMAP

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

RNA_PATH   = Path(PATHS["scrna"])    # GSE70678
OUTPUT_DIR = Path("data/cmap_query")

# Differential expression thresholds
UP_THRESHOLD   =  1.0    # z-score / log2FC
DOWN_THRESHOLD = -1.0

# clue.io recommends 50–150 genes per direction for best signal
MAX_GENES_PER_DIRECTION = 150
MIN_GENES_PER_DIRECTION = 20

# ATRT vs normal column indicators (matches GSE70678 naming)
ATRT_INDICATORS   = GENOMICS["rna_h3k27m_col_indicators"]    # ATRT, MRT, tumor...
NORMAL_INDICATORS = GENOMICS["rna_normal_col_indicators"]    # normal, brain, cortex...


def load_rna_data() -> pd.DataFrame:
    if not RNA_PATH.exists():
        raise FileNotFoundError(
            f"GSE70678 not found at {RNA_PATH}\n"
            "Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70678\n"
            "Or use GEOparse:\n"
            "  pip install GEOparse\n"
            "  python -c \"\n"
            "  import GEOparse, pandas as pd\n"
            "  gse = GEOparse.get_GEO('GSE70678', destdir='data/raw_omics/')\n"
            "  \""
        )
    logger.info(f"Loading {RNA_PATH}...")
    df = pd.read_csv(RNA_PATH, sep="\t", comment="!", index_col=0)
    df.index = df.index.astype(str).str.upper().str.strip()
    # Collapse duplicate gene symbols — take max expression
    df = df.groupby(level=0).max()
    logger.info(f"  Shape: {df.shape[0]} genes × {df.shape[1]} samples")
    return df


def identify_columns(df: pd.DataFrame):
    """Separate ATRT and normal columns using GSE70678 naming conventions."""
    atrt_cols = [
        c for c in df.columns
        if any(ind.lower() in c.lower() for ind in ATRT_INDICATORS)
        and not any(n.lower() in c.lower() for n in NORMAL_INDICATORS)
    ]
    normal_cols = [
        c for c in df.columns
        if any(n.lower() in c.lower() for n in NORMAL_INDICATORS)
    ]

    if not atrt_cols:
        # Fallback: use all columns that are not normal
        atrt_cols = [c for c in df.columns if c not in normal_cols]
        logger.warning(
            "No ATRT indicator columns found by name indicators %s\n"
            "Using %d non-normal columns as ATRT samples.\n"
            "Column sample: %s",
            ATRT_INDICATORS, len(atrt_cols), list(df.columns[:10])
        )

    logger.info(f"  ATRT samples ({len(atrt_cols)}): {atrt_cols[:5]}{'...' if len(atrt_cols) > 5 else ''}")
    logger.info(f"  Normal samples ({len(normal_cols)}): {normal_cols}")
    return atrt_cols, normal_cols


def compute_signature(
    df: pd.DataFrame, atrt_cols: list, normal_cols: list
) -> tuple:
    """Compute differential expression: ATRT mean vs normal mean."""
    df_num    = df.apply(pd.to_numeric, errors="coerce")
    atrt_mean = df_num[atrt_cols].mean(axis=1)

    if normal_cols:
        normal_mean = df_num[normal_cols].mean(axis=1)
        diff        = atrt_mean - normal_mean
        method      = f"mean(ATRT, n={len(atrt_cols)}) − mean(normal, n={len(normal_cols)})"
        logger.info(f"  Computing differential: {method}")
    else:
        diff   = atrt_mean
        method = f"mean(ATRT, n={len(atrt_cols)}) — no normal reference"
        logger.warning(
            "No normal samples found — using raw ATRT mean expression.\n"
            "This is less precise. Consider downloading GSE70678 with normal samples."
        )

    diff = diff.dropna()
    logger.info(f"  Genes with valid scores: {len(diff)}")
    logger.info(f"  Score range: [{diff.min():.2f}, {diff.max():.2f}]")
    return diff, method


def select_genes(diff: pd.Series) -> tuple:
    """Select top up- and down-regulated genes for clue.io query."""
    up_all   = diff[diff > UP_THRESHOLD].sort_values(ascending=False)
    down_all = diff[diff < DOWN_THRESHOLD].sort_values(ascending=True)

    logger.info(f"  Genes above {UP_THRESHOLD}: {len(up_all)}")
    logger.info(f"  Genes below {DOWN_THRESHOLD}: {len(down_all)}")

    up_genes   = up_all.head(MAX_GENES_PER_DIRECTION)
    down_genes = down_all.head(MAX_GENES_PER_DIRECTION)

    if len(up_genes) < MIN_GENES_PER_DIRECTION:
        logger.warning(
            f"Only {len(up_genes)} upregulated genes — query may be weak.\n"
            f"Consider lowering UP_THRESHOLD from {UP_THRESHOLD}."
        )
    if len(down_genes) < MIN_GENES_PER_DIRECTION:
        logger.warning(
            f"Only {len(down_genes)} downregulated genes. "
            "clue.io can run with UP genes only — this is still valid."
        )

    logger.info(
        f"\n  Top 10 upregulated ATRT genes:\n"
        + "\n".join(f"    {g}: {s:.2f}" for g, s in up_genes.head(10).items())
        + f"\n\n  Top 10 downregulated ATRT genes:\n"
        + "\n".join(f"    {g}: {s:.2f}" for g, s in down_genes.head(10).items())
    )

    # Validate key ATRT biology genes appear in signature
    expected_up = {"MYC", "MYCN", "EZH2", "BRD4", "AURKA", "CDK4", "HDAC1"}
    found_up = expected_up & set(up_genes.index)
    if found_up:
        logger.info(f"\n  ✅ Key ATRT targets in upregulated signature: {found_up}")
    else:
        logger.warning(
            "\n  ⚠️ None of expected ATRT targets {%s} in top upregulated genes.\n"
            "  The signature may not reflect SMARCB1-loss biology — check GSE70678 loading.",
            ", ".join(expected_up)
        )

    return up_genes, down_genes


def convert_to_entrez(gene_symbols: list) -> dict:
    """Convert gene symbols to Entrez IDs using MyGene.info API (free, no key)."""
    logger.info(f"\n  Converting {len(gene_symbols)} symbols to Entrez IDs (MyGene.info)...")
    try:
        payload  = {
            "q":       ",".join(gene_symbols),
            "scopes":  "symbol",
            "fields":  "entrezgene",
            "species": "human",
        }
        response = requests.post(
            "https://mygene.info/v3/querymany",
            data=payload,
            timeout=30,
        )
        response.raise_for_status()
        results = response.json()

        symbol_to_entrez = {}
        for hit in results:
            if "notfound" not in hit and "entrezgene" in hit:
                symbol_to_entrez[hit["query"]] = str(hit["entrezgene"])

        logger.info(f"  Converted {len(symbol_to_entrez)}/{len(gene_symbols)} symbols")
        return symbol_to_entrez

    except Exception as e:
        logger.warning(f"  MyGene.info conversion failed: {e}")
        logger.warning("  Gene symbol files will still work for the clue.io web UI")
        return {}


def save_outputs(
    up_genes: pd.Series, down_genes: pd.Series, method: str
) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    up_symbols   = list(up_genes.index)
    down_symbols = list(down_genes.index)

    # Save symbol files (for clue.io web interface)
    up_path   = OUTPUT_DIR / "atrt_up_genes.txt"
    down_path = OUTPUT_DIR / "atrt_down_genes.txt"

    with open(up_path, "w") as f:
        f.write("\n".join(up_symbols))
    with open(down_path, "w") as f:
        f.write("\n".join(down_symbols))

    logger.info(f"\n  Saved {len(up_symbols)} UP genes → {up_path}")
    logger.info(f"  Saved {len(down_symbols)} DOWN genes → {down_path}")

    # Convert and save Entrez ID files (for clue.io API)
    all_symbols  = up_symbols + down_symbols
    symbol_map   = convert_to_entrez(all_symbols)

    if symbol_map:
        up_entrez   = [symbol_map[s] for s in up_symbols if s in symbol_map]
        down_entrez = [symbol_map[s] for s in down_symbols if s in symbol_map]

        with open(OUTPUT_DIR / "atrt_up_entrez.txt", "w") as f:
            f.write("\n".join(up_entrez))
        with open(OUTPUT_DIR / "atrt_down_entrez.txt", "w") as f:
            f.write("\n".join(down_entrez))

        logger.info(
            f"  Saved Entrez ID files ({len(up_entrez)} UP, {len(down_entrez)} DOWN)"
        )

    # Save human-readable summary
    instructions = [
        "ATRT SMARCB1-loss Gene Signature for clue.io Query",
        "=" * 55,
        "",
        f"Data source: GSE70678 (Torchia 2015, Cancer Cell)",
        f"Method: {method}",
        f"Upregulated genes: {len(up_symbols)}  (threshold: > {UP_THRESHOLD})",
        f"Downregulated genes: {len(down_symbols)}  (threshold: < {DOWN_THRESHOLD})",
        "",
        "Top 20 upregulated genes (highest differential expression in ATRT vs normal):",
    ]
    for gene, score in up_genes.head(20).items():
        instructions.append(f"  {gene:<15}  score={score:.2f}")
    instructions += [
        "",
        "Top 20 downregulated genes (lowest expression in ATRT):",
    ]
    for gene, score in down_genes.head(20).items():
        instructions.append(f"  {gene:<15}  score={score:.2f}")
    instructions += [
        "",
        "=" * 55,
        "HOW TO SUBMIT TO clue.io:",
        "=" * 55,
        "1. Go to https://clue.io — create free account",
        "2. Click 'Query' → 'L1000 Query'",
        "3. UP box: paste contents of atrt_up_genes.txt",
        "4. DOWN box: paste contents of atrt_down_genes.txt",
        "5. Name: 'ATRT_SMARCB1_signature'",
        "6. Cell line: Use all available (broadest coverage)",
        "7. Click Submit — results ready in ~10 minutes",
        "8. Download: look for query_result.gct in your results",
        "9. Move to: data/cmap_query/atrt_query_result.gct",
        "10. Run: python -m backend.pipeline.integrate_cmap_results --output atrt_cmap_scores.json",
        "",
        "INTERPRETING RESULTS:",
        "  norm_cs < -0.9 = strong reversal of ATRT signature (prioritise these drugs)",
        "  norm_cs > +0.9 = mimics ATRT signature (deprioritise these drugs)",
        "  |norm_cs| < 0.5 = neutral (no significant connectivity)",
        "",
        "EXPECTED TOP REVERSERS:",
        "  Tazemetostat (EZH2-i), Panobinostat (HDAC-i), Birabresib (BET-i),",
        "  Alisertib (AURKA-i) should appear as strong reversers of SMARCB1-loss signature.",
    ]

    summary_path = OUTPUT_DIR / "atrt_query_summary.txt"
    with open(summary_path, "w") as f:
        f.write("\n".join(instructions))

    print("\n" + "\n".join(instructions))


def main():
    logger.info("=" * 55)
    logger.info("ATRT CMap Query Preparation (v1.0)")
    logger.info("Source: GSE70678 (Torchia 2015 Cancer Cell)")
    logger.info("=" * 55)

    df                      = load_rna_data()
    atrt_cols, normal_cols  = identify_columns(df)
    diff, method            = compute_signature(df, atrt_cols, normal_cols)
    up_genes, down_genes    = select_genes(diff)
    save_outputs(up_genes, down_genes, method)

    logger.info("\n✅ Done. Files saved to data/cmap_query/")
    logger.info("Next: follow instructions in atrt_query_summary.txt")


if __name__ == "__main__":
    main()