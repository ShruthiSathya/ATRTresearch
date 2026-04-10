"""
data_downloader.py  (v3.2 - FINAL FIX)
=====================================
Automated data acquisition for the ATRT Drug Repurposing Pipeline.
Corrects the pd.read_csv typo and ensures path consistency.
"""

import argparse
import gzip
import logging
import os
import shutil
import time
import urllib.request
import io as _io
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Dataset catalogue
# ─────────────────────────────────────────────────────────────────────────────

DATASETS: Dict[str, Dict] = {
    "gpl6244": {
        "name":        "GPL6244 — Affymetrix HuGene-1_0-st-v1 probe annotation",
        "url":         "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz",
        "output_path": "raw_omics/GPL6244.annot.gz",
        "processed":   "raw_omics/GPL6244_probe_map.tsv",
        "required":    True,
    },
    "gse70678": {
        "name":        "GSE70678 — Torchia 2015 ATRT bulk RNA-seq",
        "url":         "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70678/matrix/GSE70678_series_matrix.txt.gz",
        "output_path": "raw_omics/GSE70678_series_matrix.txt.gz",
        "processed":   "raw_omics/GSE70678_gene_expression.tsv",
        "required":    True,
    },
    "gtex_brain": {
        "name":        "GTEx v8 — Gene median TPM",
        "url":         "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
        "output_path": "raw_omics/GTEx_v8_gene_median_tpm.gct.gz",
        "processed":   "raw_omics/GTEx_brain_normal_reference.tsv",
        "required":    True,
    }
}

# ─────────────────────────────────────────────────────────────────────────────
# Post-processing helpers
# ─────────────────────────────────────────────────────────────────────────────

def process_gpl6244():
    """Builds probe map while skipping NCBI headers dynamically."""
    path = Path("raw_omics/GPL6244.annot.gz")
    out = Path("raw_omics/GPL6244_probe_map.tsv")
    
    logger.info(f"Parsing {path}...")
    try:
        with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
        
        start_line = 0
        for i, line in enumerate(lines):
            if line.startswith("!platform_table_begin") or "ID\t" in line:
                start_line = i if "ID\t" in line else i + 1
                break
        
        df = pd.read_csv(_io.StringIO("".join(lines[start_line:])), 
                         sep="\t", dtype=str, low_memory=False)
        
        df.columns = df.columns.str.strip()
        id_col = next(c for c in df.columns if c.lower() in ['id', 'probe set id'])
        sym_col = next(c for c in df.columns if 'symbol' in c.lower())

        probe_map = df[[id_col, sym_col]].dropna()
        probe_map.columns = ['probe_id', 'gene_symbol']
        probe_map['gene_symbol'] = probe_map['gene_symbol'].apply(
            lambda x: str(x).split('///')[0].strip().upper()
        )
        
        probe_map.to_csv(out, sep="\t", index=False)
        logger.info(f"✅ Created probe map: {out}")
        return True
    except Exception as e:
        logger.error(f"❌ process_gpl6244 failed: {e}")
        return False

def process_gse70678():
    """Maps Probe IDs to Gene Symbols for the ATRT dataset."""
    matrix_path = Path("raw_omics/GSE70678_series_matrix.txt.gz")
    map_path = Path("raw_omics/GPL6244_probe_map.tsv")
    out_path = Path("raw_omics/GSE70678_gene_expression.tsv")

    if not map_path.exists():
        logger.error("Probe map not found. Run --process gpl6244 first.")
        return False

    logger.info(f"Processing {matrix_path}...")
    try:
        # 1. Load the probe map (FIXED: removed index_index typo)
        probe_map = pd.read_csv(map_path, sep="\t")
        mapping_dict = dict(zip(probe_map.probe_id.astype(str), probe_map.gene_symbol))

        # 2. Read the series matrix (skipping the metadata header)
        df = pd.read_csv(matrix_path, sep="\t", compression='gzip', comment='!', index_col=0)

        # 3. Map probes to symbols and aggregate (mean)
        df.index = df.index.astype(str)
        df['symbol'] = df.index.map(mapping_dict)
        
        # Filter rows that have a valid gene symbol and group by symbol
        gene_df = df.dropna(subset=['symbol']).groupby('symbol').mean()

        gene_df.to_csv(out_path, sep="\t")
        logger.info(f"✅ Created gene expression file: {out_path} ({len(gene_df)} genes)")
        return True
    except Exception as e:
        logger.error(f"❌ process_gse70678 failed: {e}")
        return False

def process_gtex():
    path = Path("raw_omics/GTEx_v8_gene_median_tpm.gct.gz")
    out = Path("raw_omics/GTEx_brain_normal_reference.tsv")
    
    logger.info("Processing GTEx Brain data...")
    try:
        # Skip GCT header (2 lines)
        df = pd.read_csv(path, sep='\t', skiprows=2)
        brain_cols = [c for c in df.columns if "Brain" in c]
        res = df[['Description'] + brain_cols].copy()
        res.columns = ['gene_symbol'] + brain_cols
        res = res.groupby('gene_symbol').median()
        
        # log2(TPM + 1) transformation
        res = np.log2(res + 1)
        res.to_csv(out, sep="\t")
        logger.info(f"✅ Saved GTEx reference to {out}")
        return True
    except Exception as e:
        logger.error(f"❌ process_gtex failed: {e}")
        return False

# ─────────────────────────────────────────────────────────────────────────────
# CLI Execution
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=["gpl6244", "gse70678", "gtex_brain", "all"])
    parser.add_argument("--process", choices=["gpl6244", "gse70678", "gtex", "all"])
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    # Ensure local directory exists
    os.makedirs("raw_omics", exist_ok=True)

    if args.dataset:
        keys = ["gpl6244", "gse70678", "gtex_brain"] if args.dataset == "all" else [args.dataset]
        for k in keys:
            out = Path(DATASETS[k]["output_path"])
            logger.info(f"Downloading {DATASETS[k]['url']}...")
            urllib.request.urlretrieve(DATASETS[k]["url"], out)
            logger.info(f"✅ Downloaded {k}")

    if args.process:
        if args.process in ["gpl6244", "all"]: 
            process_gpl6244()
        if args.process in ["gse70678", "all"]: 
            process_gse70678() 
        if args.process in ["gtex", "all"]: 
            process_gtex()

    if args.check:
        print("\nReady Status:")
        for k, v in DATASETS.items():
            status = "✅" if Path(v["output_path"]).exists() else "❌"
            proc_status = "✅" if v.get("processed") and Path(v["processed"]).exists() else "❌"
            print(f"{status} {k} (Raw) | {proc_status} {k} (Processed)")

if __name__ == "__main__":
    main()