
import pandas as pd
import numpy as np
import json
import logging
import requests
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION — adjust these if needed
# ─────────────────────────────────────────────────────────────────────────────

RNA_PATH     = Path("data/validation/cbtn_genomics/rna_zscores.txt")
OUTPUT_DIR   = Path("data/cmap_query")

# Thresholds (matching your pipeline_config.py GENOMICS settings)
UP_THRESHOLD   =  1.0   # z-score > 1.0 = upregulated in H3K27M vs normal
DOWN_THRESHOLD = -1.0   # z-score < -1.0 = downregulated

# clue.io works best with 50-150 genes per direction.
# Too many genes = noisy query. Too few = no signal.
MAX_GENES_PER_DIRECTION = 150
MIN_GENES_PER_DIRECTION = 10

# H3K27M tumour column indicators (matching your pipeline)
H3K27M_INDICATORS = ["Pons", "K27M", "DIPG", "H3K27M", "pontine"]
NORMAL_INDICATORS = ["normal", "control", "cortex", "healthy"]


def load_rna_data() -> pd.DataFrame:
    if not RNA_PATH.exists():
        raise FileNotFoundError(
            f"rna_zscores.txt not found at {RNA_PATH}\n"
            "Make sure you're running from the project root (gbmresearch/)"
        )
    logger.info(f"Loading {RNA_PATH}...")
    df = pd.read_csv(RNA_PATH, sep="\t", index_col=0)
    logger.info(f"  Shape: {df.shape[0]} genes × {df.shape[1]} samples")
    return df


def identify_columns(df: pd.DataFrame):
    """Separate H3K27M tumour columns from normal columns."""
    h3k27m_cols = [
        c for c in df.columns
        if any(ind.lower() in c.lower() for ind in H3K27M_INDICATORS)
        and not any(n.lower() in c.lower() for n in NORMAL_INDICATORS)
    ]
    normal_cols = [
        c for c in df.columns
        if any(n.lower() in c.lower() for n in NORMAL_INDICATORS)
    ]

    if not h3k27m_cols:
        # Fallback: use all non-normal columns
        h3k27m_cols = [c for c in df.columns if c not in normal_cols]
        logger.warning(
            "No H3K27M indicator columns found by name — "
            f"using all {len(h3k27m_cols)} non-normal columns"
        )

    logger.info(f"  H3K27M tumour columns ({len(h3k27m_cols)}): {h3k27m_cols}")
    logger.info(f"  Normal columns ({len(normal_cols)}): {normal_cols}")
    return h3k27m_cols, normal_cols


def compute_signature(df: pd.DataFrame, h3k27m_cols, normal_cols):
    """
    Compute differential expression signature.

    If both tumour and normal columns exist: use mean(tumour) - mean(normal).
    If only tumour columns: use mean z-score directly (already normalised).
    """
    # Convert to numeric
    df_numeric = df.apply(pd.to_numeric, errors="coerce")

    h3k27m_mean = df_numeric[h3k27m_cols].mean(axis=1)

    if normal_cols:
        normal_mean  = df_numeric[normal_cols].mean(axis=1)
        diff         = h3k27m_mean - normal_mean
        method       = "mean(H3K27M) - mean(normal)"
        logger.info(f"  Computing differential: {method}")
    else:
        diff   = h3k27m_mean
        method = "mean(H3K27M) z-score (no normal reference)"
        logger.warning(
            f"No normal columns found — using raw H3K27M mean z-scores. "
            f"This is less precise than a differential."
        )

    diff = diff.dropna()
    logger.info(f"  Genes with valid scores: {len(diff)}")
    return diff, method


def select_genes(diff: pd.Series):
    """
    Select top up and down genes for clue.io query.

    clue.io recommendation: 50-150 genes per direction.
    Uses the most extreme scores for best signal.
    """
    up_all   = diff[diff > UP_THRESHOLD].sort_values(ascending=False)
    down_all = diff[diff < DOWN_THRESHOLD].sort_values(ascending=True)

    logger.info(f"  Genes above {UP_THRESHOLD}: {len(up_all)}")
    logger.info(f"  Genes below {DOWN_THRESHOLD}: {len(down_all)}")

    # Take top N most extreme
    up_genes   = up_all.head(MAX_GENES_PER_DIRECTION)
    down_genes = down_all.head(MAX_GENES_PER_DIRECTION)

    if len(up_genes) < MIN_GENES_PER_DIRECTION:
        logger.warning(
            f"Only {len(up_genes)} upregulated genes — query may be weak. "
            f"Consider lowering UP_THRESHOLD from {UP_THRESHOLD}."
        )
    if len(down_genes) < MIN_GENES_PER_DIRECTION:
        logger.warning(
            f"Only {len(down_genes)} downregulated genes found. "
            f"clue.io can run with UP genes only — this is still valid."
        )

    return up_genes, down_all.head(MAX_GENES_PER_DIRECTION)


def convert_to_entrez(gene_symbols: list) -> dict:
    """
    Convert gene symbols to Entrez IDs using MyGene.info API.
    Free, no API key required. Returns {symbol: entrez_id}.
    """
    logger.info(f"  Converting {len(gene_symbols)} symbols to Entrez IDs via MyGene.info...")
    try:
        url      = "https://mygene.info/v3/querymany"
        payload  = {
            "q":       ",".join(gene_symbols),
            "scopes":  "symbol",
            "fields":  "entrezgene",
            "species": "human",
        }
        response = requests.post(url, data=payload, timeout=30)
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
        logger.warning("  Symbol files will still work for the clue.io web interface")
        return {}


def save_outputs(up_genes: pd.Series, down_genes: pd.Series, method: str):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    up_symbols   = list(up_genes.index)
    down_symbols = list(down_genes.index)

    # ── Save gene symbol files ──────────────────────────────────────────────
    up_path   = OUTPUT_DIR / "h3k27m_up_genes.txt"
    down_path = OUTPUT_DIR / "h3k27m_down_genes.txt"

    with open(up_path, "w") as f:
        f.write("\n".join(up_symbols))
    with open(down_path, "w") as f:
        f.write("\n".join(down_symbols))

    logger.info(f"  Saved {len(up_symbols)} UP genes → {up_path}")
    logger.info(f"  Saved {len(down_symbols)} DOWN genes → {down_path}")

    # ── Convert and save Entrez ID files ────────────────────────────────────
    all_symbols = up_symbols + down_symbols
    symbol_map  = convert_to_entrez(all_symbols)

    if symbol_map:
        up_entrez   = [symbol_map[s] for s in up_symbols if s in symbol_map]
        down_entrez = [symbol_map[s] for s in down_symbols if s in symbol_map]

        with open(OUTPUT_DIR / "h3k27m_up_entrez.txt", "w") as f:
            f.write("\n".join(up_entrez))
        with open(OUTPUT_DIR / "h3k27m_down_entrez.txt", "w") as f:
            f.write("\n".join(down_entrez))

        logger.info(f"  Saved Entrez ID files (needed for API method)")

    # ── Save human-readable summary ──────────────────────────────────────────
    summary = [
        "H3K27M DIPG Gene Signature for clue.io Query",
        "=" * 50,
        f"Method: {method}",
        f"Upregulated genes: {len(up_symbols)}",
        f"Downregulated genes: {len(down_symbols)}",
        "",
        "Top 20 upregulated genes (highest z-score):",
    ]
    for gene, score in up_genes.head(20).items():
        summary.append(f"  {gene:15s}  z={score:.2f}")
    summary += [
        "",
        "Top 20 downregulated genes (lowest z-score):",
    ]
    for gene, score in down_genes.head(20).items():
        summary.append(f"  {gene:15s}  z={score:.2f}")
    summary += [
        "",
        "HOW TO USE ON clue.io:",
        "1. Go to https://clue.io and create a free account",
        "2. Click 'Query' in the top menu",
        "3. In the UP box: paste contents of h3k27m_up_genes.txt",
        "4. In the DOWN box: paste contents of h3k27m_down_genes.txt",
        "5. Set name: 'H3K27M_DIPG_signature'",
        "6. Click Submit — results ready in ~10 minutes",
        "7. Download results → save as data/cmap_query/cmap_results.json",
        "8. Run: python integrate_cmap_results.py",
    ]

    with open(OUTPUT_DIR / "query_summary.txt", "w") as f:
        f.write("\n".join(summary))

    print("\n" + "\n".join(summary))
    return up_symbols, down_symbols


def main():
    logger.info("=" * 55)
    logger.info("Preparing H3K27M DIPG signature for clue.io")
    logger.info("=" * 55)

    df                   = load_rna_data()
    h3k27m_cols, normal_cols = identify_columns(df)
    diff, method         = compute_signature(df, h3k27m_cols, normal_cols)
    up_genes, down_genes = select_genes(diff)
    up_symbols, _        = save_outputs(up_genes, down_genes, method)

    logger.info("\n✅ Done. Files saved to data/cmap_query/")
    logger.info("Next: follow the instructions in query_summary.txt")


if __name__ == "__main__":
    main()