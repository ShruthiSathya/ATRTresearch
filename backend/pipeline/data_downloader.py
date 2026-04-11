"""
data_downloader.py  (v3.3 - ALL BUGS FIXED)
============================================
FIXES FROM v3.2:
  1. All paths are now ABSOLUTE (rooted at repo root) — eliminates CWD sensitivity.
  2. GTEx URL corrected: adult-gtex bucket → gtex_analysis_v8 bucket (old URL 404'd).
  3. process_gse70678() now uses parse_geo_series_matrix() from probe_mapper.py
     instead of comment='!' hack (which silently dropped the ID_REF header row
     when the file had a blank line before the data table begin marker).
  4. os.makedirs replaced with Path.mkdir(parents=True, exist_ok=True) throughout.
  5. data/ prefix added to all output_path / processed entries in DATASETS dict.
  6. --check now checks both raw and processed file existence.

HOW TO USE
----------
From the REPO ROOT directory:

  # Step 1 — download raw files (~2 GB total)
  python -m backend.pipeline.data_downloader --dataset all

  # Step 2 — process raw → gene-symbol TSVs
  python -m backend.pipeline.data_downloader --process all

  # Verify all files exist
  python -m backend.pipeline.data_downloader --check

  # Individual steps
  python -m backend.pipeline.data_downloader --dataset gpl6244 gse70678
  python -m backend.pipeline.data_downloader --process gpl6244 gse70678

FILES PRODUCED
--------------
  data/raw_omics/GPL6244.annot.gz                    (raw)
  data/raw_omics/GPL6244_probe_map.tsv               (processed — probe → gene symbol)
  data/raw_omics/GSE70678_series_matrix.txt.gz       (raw)
  data/raw_omics/GSE70678_gene_expression.tsv        (processed — gene × sample, log2)
  data/raw_omics/GTEx_v8_gene_median_tpm.gct.gz      (raw)
  data/raw_omics/GTEx_brain_normal_reference.tsv     (processed — normal brain log2(TPM+1))

DATA SOURCES
------------
GPL6244   : NCBI GEO FTP — Affymetrix HuGene 1.0 ST annotation
GSE70678  : NCBI GEO FTP — Torchia 2015 ATRT bulk RNA-seq (49 samples)
GTEx v8   : GTEx portal — gene median TPM across 54 tissues

REFERENCES
----------
Torchia J et al. Cancer Cell 2015; 30(6):891. PMID 26609405.
GTEx Consortium. Science 2020; 369(6509):1318. PMID 32913098.
Irizarry RA et al. Biostatistics 2003 — RMA normalisation for Affymetrix.
"""

import argparse
import gzip
import logging
import sys
import urllib.request
import io as _io
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Repo-root-relative paths  (FIX: was CWD-relative 'raw_omics/...')
# ─────────────────────────────────────────────────────────────────────────────

REPO_ROOT = Path(__file__).resolve().parent.parent.parent   # backend/pipeline/ → repo root
DATA_DIR  = REPO_ROOT / "data" / "raw_omics"

# ─────────────────────────────────────────────────────────────────────────────
# Dataset catalogue
# ─────────────────────────────────────────────────────────────────────────────

DATASETS: Dict[str, Dict] = {
    "gpl6244": {
        "name":        "GPL6244 — Affymetrix HuGene-1_0-st-v1 probe annotation",
        # NCBI GEO FTP — stable URL, confirmed working April 2026
        "url":         "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz",
        "output_path": str(DATA_DIR / "GPL6244.annot.gz"),
        "processed":   str(DATA_DIR / "GPL6244_probe_map.tsv"),
        "required":    True,
    },
    "gse70678": {
        "name":        "GSE70678 — Torchia 2015 ATRT bulk RNA-seq (49 ATRT + normals)",
        # NCBI GEO FTP — stable URL
        "url":         "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70678/matrix/GSE70678_series_matrix.txt.gz",
        "output_path": str(DATA_DIR / "GSE70678_series_matrix.txt.gz"),
        "processed":   str(DATA_DIR / "GSE70678_gene_expression.tsv"),
        "required":    True,
    },
    "gtex_brain": {
        "name":        "GTEx v8 — Gene median TPM across 54 tissues",
        # FIX: old URL (adult-gtex bucket) returns 404.
        # Correct v8 URL confirmed at gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
        "url":         "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
        "output_path": str(DATA_DIR / "GTEx_v8_gene_median_tpm.gct.gz"),
        "processed":   str(DATA_DIR / "GTEx_brain_normal_reference.tsv"),
        "required":    True,
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Post-processing helpers
# ─────────────────────────────────────────────────────────────────────────────

def process_gpl6244() -> bool:
    """
    Build probe → gene symbol mapping from GPL6244 annotation.
    Skips NCBI comment header lines (start with #) then reads the
    platform table. Multi-gene probes use the first symbol.
    """
    path = Path(DATASETS["gpl6244"]["output_path"])
    out  = Path(DATASETS["gpl6244"]["processed"])

    if not path.exists():
        logger.error("GPL6244 raw file not found: %s  —  run --dataset gpl6244 first", path)
        return False

    logger.info("Parsing %s ...", path)
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        # Find the data table start — skip all comment / metadata lines
        # GPL annotation format:
        #   Lines 1-N : metadata starting with #
        #   First non-# line : column headers (ID\tGene symbol\t...)
        #   Subsequent lines : data rows
        start_line = 0
        for i, line in enumerate(lines):
            stripped = line.strip()
            if not stripped.startswith("#") and stripped:
                start_line = i
                break

        if start_line == 0 and lines[0].startswith("#"):
            logger.error("Could not find data table start in GPL6244 annotation")
            return False

        df = pd.read_csv(
            _io.StringIO("".join(lines[start_line:])),
            sep="\t",
            dtype=str,
            low_memory=False,
        )
        df.columns = df.columns.str.strip()

        # Locate probe ID and gene symbol columns (handles minor column name variants)
        id_col  = next((c for c in df.columns if c.lower() in ["id", "probe set id"]), None)
        sym_col = next((c for c in df.columns if "symbol" in c.lower()), None)
        if not id_col or not sym_col:
            logger.error("Could not find ID/symbol columns. Available: %s", list(df.columns[:10]))
            return False

        probe_map = df[[id_col, sym_col]].dropna()
        probe_map.columns = ["probe_id", "gene_symbol"]
        # Multi-gene probes: "BRCA1 /// BRCA2" → "BRCA1"
        probe_map["gene_symbol"] = probe_map["gene_symbol"].apply(
            lambda x: str(x).split("///")[0].strip().upper()
        )
        # Drop unmapped / control probes
        probe_map = probe_map[
            probe_map["gene_symbol"].str.len() > 0
        ]
        probe_map = probe_map[~probe_map["probe_id"].str.startswith("AFFX-", na=False)]

        out.parent.mkdir(parents=True, exist_ok=True)
        probe_map.to_csv(out, sep="\t", index=False)
        logger.info("✅ GPL6244 probe map: %d probes → %s", len(probe_map), out)
        return True

    except Exception as e:
        logger.exception("process_gpl6244 failed: %s", e)
        return False


def process_gse70678() -> bool:
    """
    Map Affymetrix probe IDs → HGNC gene symbols for GSE70678.

    FIX vs v3.2:
      - Uses parse_geo_series_matrix() from probe_mapper.py to correctly
        extract the expression table from between the GEO matrix table markers.
      - The old approach (comment='!') failed because pd.read_csv drops ALL
        rows starting with '!', including the blank line that sometimes
        precedes the ID_REF header, causing a silent parse failure.
    """
    matrix_path = Path(DATASETS["gse70678"]["output_path"])
    map_path    = Path(DATASETS["gpl6244"]["processed"])
    out_path    = Path(DATASETS["gse70678"]["processed"])

    if not matrix_path.exists():
        logger.error("GSE70678 matrix not found: %s  —  run --dataset gse70678 first", matrix_path)
        return False
    if not map_path.exists():
        logger.error("Probe map not found: %s  —  run --process gpl6244 first", map_path)
        return False

    logger.info("Processing %s ...", matrix_path)
    try:
        # Import here to avoid circular deps when run standalone
        try:
            from .probe_mapper import parse_geo_series_matrix, ProbeMapper
        except ImportError:
            # Running as a standalone script from the pipeline directory
            sys.path.insert(0, str(Path(__file__).parent))
            from probe_mapper import parse_geo_series_matrix, ProbeMapper

        # Step 1: Parse series matrix → probe DataFrame
        probe_df, meta = parse_geo_series_matrix(matrix_path)
        logger.info("  Series matrix: %d probes × %d samples", len(probe_df), len(probe_df.columns))

        # Step 2: Load probe → gene map
        mapper = ProbeMapper(cache_path=map_path)
        mapper._load_from_cache()   # already processed — skip download

        # Step 3: Map and aggregate (median per gene)
        gene_df = mapper.map_probes_to_genes(probe_df, aggregation="median")
        logger.info("  After mapping: %d genes × %d samples", len(gene_df), len(gene_df.columns))

        # Step 4: Log some QC stats
        n_atrt   = sum(1 for c in gene_df.columns if any(
            k in c.upper() for k in ["ATRT", "MRT", "TUMOR", "T"]
        ))
        n_normal = sum(1 for c in gene_df.columns if any(
            k in c.upper() for k in ["NORMAL", "BRAIN", "CTRL"]
        ))
        logger.info("  Approx sample types — ATRT-like: %d | Normal-like: %d", n_atrt, n_normal)

        # Step 5: Save
        out_path.parent.mkdir(parents=True, exist_ok=True)
        gene_df.to_csv(out_path, sep="\t")
        logger.info("✅ GSE70678 gene expression saved: %s", out_path)
        return True

    except Exception as e:
        logger.exception("process_gse70678 failed: %s", e)
        return False


def process_gtex() -> bool:
    """
    Extract brain-specific median TPM from GTEx v8 GCT file.

    GTEx GCT format (v1.2):
      Line 1 : #1.2
      Line 2 : nDataRows  nDataCols
      Line 3 : Name  Description  <sample1>  <sample2> ...
      Line 4+ : data rows
    skiprows=2 is correct — reads line 3 as the column header.

    Output: log2(median_TPM + 1) for brain tissues only.
    """
    path = Path(DATASETS["gtex_brain"]["output_path"])
    out  = Path(DATASETS["gtex_brain"]["processed"])

    if not path.exists():
        logger.error("GTEx raw file not found: %s  —  run --dataset gtex_brain first", path)
        return False

    logger.info("Processing GTEx v8 gene median TPM (%s) ...", path)
    try:
        # GCT v1.2: skip lines 0 (#1.2) and 1 (dimensions) — line 2 is header
        df = pd.read_csv(path, sep="\t", skiprows=2, compression="gzip", low_memory=False)
        logger.info("  GTEx full matrix: %d rows × %d columns", len(df), len(df.columns))

        # Brain columns for ATRT normal-brain reference
        # Prioritise cerebellum (infratentorial) + cortex (supratentorial)
        brain_keywords = ["brain", "cerebellum", "cortex", "frontal"]
        brain_cols = [c for c in df.columns
                      if any(k in c.lower() for k in brain_keywords)]

        if not brain_cols:
            logger.warning("No brain columns found — using all tissues as reference")
            brain_cols = [c for c in df.columns if c not in ("Name", "Description")]

        logger.info("  Brain columns selected (%d): %s", len(brain_cols),
                    [c.replace("Brain - ", "") for c in brain_cols[:6]])

        # Use 'Description' as gene symbol index (= HGNC symbol in GTEx v8)
        if "Description" not in df.columns:
            logger.error("'Description' column not found. Available: %s", list(df.columns[:8]))
            return False

        res = df[["Description"] + brain_cols].copy()
        res.columns = ["gene_symbol"] + brain_cols

        # Numeric conversion and groupby gene symbol (some symbols duplicated)
        for col in brain_cols:
            res[col] = pd.to_numeric(res[col], errors="coerce")

        res = res.groupby("gene_symbol").median(numeric_only=True)

        # log2(TPM + 1) transformation — same scale as GSE70678 RMA values
        res = np.log2(res + 1)

        out.parent.mkdir(parents=True, exist_ok=True)
        res.to_csv(out, sep="\t")
        logger.info("✅ GTEx brain reference: %d genes × %d tissues → %s",
                    len(res), len(brain_cols), out)
        return True

    except Exception as e:
        logger.exception("process_gtex failed: %s", e)
        return False


def download_dataset(key: str) -> bool:
    """Download a single dataset by key. Returns True on success."""
    if key not in DATASETS:
        logger.error("Unknown dataset key: %s. Valid keys: %s", key, list(DATASETS.keys()))
        return False

    ds       = DATASETS[key]
    out_path = Path(ds["output_path"])
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists():
        size_mb = out_path.stat().st_size / 1e6
        logger.info("Already downloaded: %s (%.1f MB) — skipping", out_path.name, size_mb)
        return True

    logger.info("Downloading %s ...", ds["name"])
    logger.info("  URL: %s", ds["url"])
    try:
        urllib.request.urlretrieve(ds["url"], str(out_path))
        size_mb = out_path.stat().st_size / 1e6
        logger.info("✅ %s (%.1f MB) → %s", key, size_mb, out_path)
        return True
    except Exception as e:
        logger.error("Download failed for %s: %s", key, e)
        if out_path.exists():
            out_path.unlink()   # clean up partial download
        return False


def check_files() -> None:
    """Print status of all raw and processed files."""
    print("\n" + "=" * 65)
    print("DATA FILE STATUS")
    print("=" * 65)
    all_ok = True
    for key, ds in DATASETS.items():
        raw_path  = Path(ds["output_path"])
        proc_path = Path(ds["processed"]) if ds.get("processed") else None

        raw_ok  = raw_path.exists()
        proc_ok = proc_path.exists() if proc_path else None

        raw_str  = f"✅ ({raw_path.stat().st_size/1e6:.1f} MB)" if raw_ok else "❌ MISSING"
        proc_str = (f"✅ ({proc_path.stat().st_size/1e3:.0f} KB)" if proc_ok
                    else "❌ MISSING" if proc_ok is False
                    else "—")

        if not raw_ok or proc_ok is False:
            all_ok = False

        print(f"\n  [{key}] {ds['name'][:55]}")
        print(f"    Raw      : {raw_str}")
        if proc_path:
            print(f"    Processed: {proc_str}")

    print("\n" + "=" * 65)
    if all_ok:
        print("✅ All files present — pipeline is ready to run.")
    else:
        print("⚠️  Missing files. Run:")
        print("   python -m backend.pipeline.data_downloader --dataset all")
        print("   python -m backend.pipeline.data_downloader --process all")
    print("=" * 65 + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="ATRT Pipeline Data Downloader v3.3",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full setup (download + process all datasets)
  python -m backend.pipeline.data_downloader --dataset all --process all

  # Individual steps
  python -m backend.pipeline.data_downloader --dataset gpl6244 gse70678 gtex_brain
  python -m backend.pipeline.data_downloader --process gpl6244 gse70678 gtex

  # Check file status
  python -m backend.pipeline.data_downloader --check
        """,
    )
    parser.add_argument(
        "--dataset",
        nargs="+",
        choices=list(DATASETS.keys()) + ["all"],
        help="Dataset(s) to download",
    )
    parser.add_argument(
        "--process",
        nargs="+",
        choices=["gpl6244", "gse70678", "gtex", "all"],
        help="Dataset(s) to process (raw → gene TSV)",
    )
    parser.add_argument("--check", action="store_true", help="Show file status and exit")
    args = parser.parse_args()

    if args.check:
        check_files()
        return

    if not args.dataset and not args.process:
        parser.print_help()
        return

    # Download
    if args.dataset:
        keys = list(DATASETS.keys()) if "all" in args.dataset else args.dataset
        for k in keys:
            download_dataset(k)

    # Process
    if args.process:
        steps = args.process
        run_gpl6244 = "all" in steps or "gpl6244" in steps
        run_gse70678 = "all" in steps or "gse70678" in steps
        run_gtex    = "all" in steps or "gtex" in steps

        # Order matters: gpl6244 must run before gse70678
        if run_gpl6244:
            process_gpl6244()
        if run_gse70678:
            process_gse70678()
        if run_gtex:
            process_gtex()

    # Final status
    check_files()


if __name__ == "__main__":
    main()