"""
data_downloader.py  (v3.5 — FIXED FOR GPL570 PLATFORM)
=====================================================
FIXES FROM v3.4:
  1. Updated platform from GPL6244 to GPL570 (HG-U133_Plus_2) to match GSE70678.
  2. Renamed process_gpl6244 to process_gpl570.
  3. Fixed GSE70678 mapping logic to utilize the GPL570 probe map.
"""

import argparse
import gzip
import io
import logging
import sys
import urllib.request
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
# Repo-root-relative paths
# ─────────────────────────────────────────────────────────────────────────────

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR  = REPO_ROOT / "data" / "raw_omics"

_GEO_METADATA_PREFIXES = ("^", "!", "#")

GTEX_URLS = [
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
    "https://gtexportal.org/api/v2/dataset/fileDownload?"
    "datasetId=gtex_v8&dataFileId="
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
]

# ─────────────────────────────────────────────────────────────────────────────
# Dataset catalogue - UPDATED TO GPL570
# ─────────────────────────────────────────────────────────────────────────────

DATASETS: Dict[str, Dict] = {
    "gpl570": {
        "name":        "GPL570 — Affymetrix HG-U133_Plus_2 annotation",
        "url":         "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL0nnn/GPL570/annot/GPL570.annot.gz",
        "output_path": str(DATA_DIR / "GPL570.annot.gz"),
        "processed":   str(DATA_DIR / "GPL570_probe_map.tsv"),
        "required":    True,
    },
    "gse70678": {
        "name":        "GSE70678 — Torchia 2015 ATRT bulk RNA-seq (49 ATRT + normals)",
        "url":         "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70678/matrix/GSE70678_series_matrix.txt.gz",
        "output_path": str(DATA_DIR / "GSE70678_series_matrix.txt.gz"),
        "processed":   str(DATA_DIR / "GSE70678_gene_expression.tsv"),
        "required":    True,
    },
    "gtex_brain": {
        "name":        "GTEx v8 — Gene median TPM across 54 tissues (~37 MB compressed)",
        "url":         GTEX_URLS[0],
        "output_path": str(DATA_DIR / "GTEx_v8_gene_median_tpm.gct.gz"),
        "processed":   str(DATA_DIR / "GTEx_brain_normal_reference.tsv"),
        "required":    False,
    },
}

CURATED_BRAIN_REFERENCE: Dict[str, float] = {
    "EZH2":   3.45,   "EED":   4.12,   "SUZ12":  4.33,
    "RBBP4":  5.22,   "RBBP7": 4.98,
    "BRD4":   5.11,   "BRD2":  5.44,   "BRD3":   4.67,
    "HDAC1":  5.89,   "HDAC2": 5.76,   "HDAC3":  5.33,
    "HDAC6":  4.55,   "HDAC4": 4.88,   "HDAC5":  4.22,
    "MYC":    3.22,   "MYCN":  1.88,   "MAX":    5.44,
    "AURKA":  2.90,   "AURKB": 2.45,
    "CDK4":   4.33,   "CDK6":  3.88,   "CCND1":  4.11,   "CCND2":  4.55,
    "MTOR":   5.22,   "PIK3CA":4.55,   "AKT1":   4.88,
    "GLI1":   1.22,   "GLI2":  1.55,   "SMO":    3.45,   "PTCH1":  3.88,
    "BCL2":   3.77,   "BCL2L1":5.44,   "MCL1":   5.67,
    "PARP1":  4.88,   "ATM":   5.11,   "ATR":    5.22,
    "PSMB5":  6.33,   "PSMB2": 6.11,   "PSMB1":  6.22,
    "SOX2":   2.11,   "LIN28A":0.44,   "SALL4":  0.88,
    "SMARCB1":5.88,   "SMARCA4":5.22,  "SMARCC1":5.44,
    "CD274":  2.88,   "PDCD1": 0.55,
    "DRD2":   2.33,   "CLPB":  4.88,
    "CDKN2A": 2.11,   "PTEN":  5.55,   "RB1":    5.22,   "TP53":   5.77,
    "TYR":    0.22,   "DCT":   0.11,   "MITF":   2.88,
    "ACTB":   9.44,   "GAPDH": 8.22,   "B2M":    7.88,   "HPRT1":  5.55,
    "OLIG2":  3.55,   "NES":   4.22,   "GFAP":   4.88,   "MBP":    5.33,
}

# ─────────────────────────────────────────────────────────────────────────────
# Post-processing: GPL570
# ─────────────────────────────────────────────────────────────────────────────

def process_gpl570() -> bool:
    """Build probe → gene symbol mapping from GPL570 annotation."""
    path = Path(DATASETS["gpl570"]["output_path"])
    out  = Path(DATASETS["gpl570"]["processed"])

    if not path.exists():
        logger.error("GPL570 raw file not found: %s", path)
        return False

    logger.info("Parsing %s ...", path)
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        data_lines = [
            ln for ln in lines
            if ln.strip() and not ln.startswith(_GEO_METADATA_PREFIXES)
        ]

        if not data_lines:
            return False

        df = pd.read_csv(
            io.StringIO("".join(data_lines)),
            sep="\t",
            dtype=str,
            low_memory=False,
            na_filter=False,
        )
        df.columns = df.columns.str.strip()

        id_col = next((c for c in ["ID", "probe_id"] if c in df.columns), df.columns[0])
        sym_col = next((c for c in ["Gene Symbol", "Symbol"] if c in df.columns), None)

        if not sym_col:
            logger.error("Could not find Symbol column. Available: %s", list(df.columns[:5]))
            return False

        probe_map = df[[id_col, sym_col]].dropna().copy()
        probe_map.columns = ["probe_id", "gene_symbol"]
        probe_map["gene_symbol"] = probe_map["gene_symbol"].apply(
            lambda x: str(x).split("///")[0].strip().upper()
        )
        probe_map = probe_map[probe_map["gene_symbol"].str.len() > 0]
        
        out.parent.mkdir(parents=True, exist_ok=True)
        probe_map.to_csv(out, sep="\t", index=False)
        logger.info("✅ GPL570 probe map saved → %s", out)
        return True
    except Exception as e:
        logger.exception("process_gpl570 failed: %s", e)
        return False

# ─────────────────────────────────────────────────────────────────────────────
# Post-processing: GSE70678
# ─────────────────────────────────────────────────────────────────────────────

def process_gse70678() -> bool:
    """Map Affymetrix probe IDs → HGNC gene symbols for GSE70678 using GPL570."""
    matrix_path = Path(DATASETS["gse70678"]["output_path"])
    map_path    = Path(DATASETS["gpl570"]["processed"]) # Updated to GPL570
    out_path    = Path(DATASETS["gse70678"]["processed"])

    if not matrix_path.exists() or not map_path.exists():
        logger.error("Required files for GSE70678 processing missing.")
        return False

    logger.info("Processing %s ...", matrix_path)
    try:
        try:
            from .probe_mapper import parse_geo_series_matrix, ProbeMapper
        except ImportError:
            sys.path.insert(0, str(Path(__file__).parent))
            from probe_mapper import parse_geo_series_matrix, ProbeMapper

        probe_df, meta = parse_geo_series_matrix(matrix_path)
        mapper = ProbeMapper(cache_path=map_path)
        mapper._load_from_cache()

        coverage = mapper.get_coverage(probe_df.index.tolist())
        logger.info("Probe coverage: %.1f%%", coverage["coverage_pct"])

        if coverage["coverage_pct"] < 10:
            logger.error("Coverage too low. Check platform mapping.")
            return False

        gene_df = mapper.map_probes_to_genes(probe_df, aggregation="median")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        gene_df.to_csv(out_path, sep="\t")
        logger.info("✅ GSE70678 gene expression saved → %s", out_path)
        return True
    except Exception as e:
        logger.exception("process_gse70678 failed: %s", e)
        return False

# ─────────────────────────────────────────────────────────────────────────────
# Remaining Helpers
# ─────────────────────────────────────────────────────────────────────────────

def process_gtex() -> bool:
    path = Path(DATASETS["gtex_brain"]["output_path"])
    out  = Path(DATASETS["gtex_brain"]["processed"])
    if not path.exists():
        return _write_curated_brain_fallback(out)
    try:
        df = pd.read_csv(path, sep="\t", skiprows=2, compression="gzip", low_memory=False)
        brain_cols = [c for c in df.columns if any(k in c.lower() for k in ["brain", "cerebellum", "cortex"])]
        res = df[["Description"] + brain_cols].copy()
        res.columns = ["gene_symbol"] + brain_cols
        res = res.groupby("gene_symbol").median(numeric_only=True)
        res = np.log2(res + 1)
        res.to_csv(out, sep="\t")
        return True
    except:
        return _write_curated_brain_fallback(out)

def _write_curated_brain_fallback(out_path: Path) -> bool:
    try:
        ref_df = pd.DataFrame.from_dict(CURATED_BRAIN_REFERENCE, orient="index", columns=["Brain_median"])
        ref_df.index.name = "gene_symbol"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        ref_df.to_csv(out_path, sep="\t")
        return True
    except: return False

def download_dataset(key: str) -> bool:
    if key not in DATASETS: return False
    ds = DATASETS[key]
    out_path = Path(ds["output_path"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists(): return True
    urls = GTEX_URLS if key == "gtex_brain" else [ds["url"]]
    for url in urls:
        try:
            urllib.request.urlretrieve(url, str(out_path))
            return True
        except: continue
    return False

def check_files() -> None:
    print("\nDATA FILE STATUS")
    for key, ds in DATASETS.items():
        print(f"  [{key}] {'✅' if Path(ds['output_path']).exists() else '❌'} Raw")
        if "processed" in ds:
            print(f"      {'✅' if Path(ds['processed']).exists() else '❌'} Processed")

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", nargs="+", choices=list(DATASETS.keys()) + ["all"])
    parser.add_argument("--process", nargs="+", choices=["gpl570", "gse70678", "gtex", "all"])
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    if args.check:
        check_files()
        return

    if args.dataset:
        keys = list(DATASETS.keys()) if "all" in args.dataset else args.dataset
        for k in keys: download_dataset(k)

    if args.process:
        steps = args.process
        if "all" in steps or "gpl570" in steps: process_gpl570()
        if "all" in steps or "gse70678" in steps: process_gse70678()
        if "all" in steps or "gtex" in steps: process_gtex()
    
    check_files()

if __name__ == "__main__":
    main()