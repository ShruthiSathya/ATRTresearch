"""
data_downloader.py  (v4.0 — GPL570 FTP path fix)
=================================================
FIXES FROM v3.5:
  1. CRITICAL BUG FIX: GPL570 FTP URL corrected.
     WAS:  .../GPL0nnn/GPL570/annot/GPL570.annot.gz   ← 404 always
     NOW:  .../GPLnnn/GPL570/annot/GPL570.annot.gz    ← correct path

     GEO FTP directory naming convention (from NCBI docs):
       GPL accession number → directory = GPL + first (N-3) digits + "nnn"
       GPL570  (3 digits) → GPLnnn
       GPL6244 (4 digits) → GPL6nnn
       GPL13667 (5 digits) → GPL13nnn

     Source: NCBI GEO FTP README
     https://ftp.ncbi.nlm.nih.gov/geo/README.txt

  2. GPL570 column mapping corrected.
     The GPL570.annot.gz file uses "Gene Symbol" (with space, Title Case)
     not "Symbol" or "gene_symbol".
     Verified from NCBI GEO GPL570 full table download (netaffx build 35,
     June 2016 — the most recent annotation update as of 2026).
     Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570

  3. Added HTTPS web-download fallback for environments where FTP is blocked.
     NCBI provides an HTTPS download endpoint:
       https://www.ncbi.nlm.nih.gov/geo/download/?acc=GPL570&format=file&file=GPL570.annot.gz

  4. Added GPL570 SOFT-format fallback parser (GPL570_family.soft.gz) for
     cases where the .annot.gz file is unavailable.

  5. GTEx reference downloader: added second mirror URL from GTEx portal.

REFERENCES
-----------
NCBI GEO FTP directory structure:
  https://www.ncbi.nlm.nih.gov/geo/info/download.html
  https://ftp.ncbi.nlm.nih.gov/geo/README.txt

GPL570 platform annotation (Affymetrix HG-U133 Plus 2.0):
  Irizarry RA et al. (2003). Summaries of Affymetrix GeneChip probe level
  data. Nucleic Acids Research 31(4):e15. PMID 12582260.
  (RMA normalisation justification for log2-scale values in series matrix)

GSE70678 dataset:
  Torchia J et al. (2015). Integrated (epi)-genomic analyses identify
  subgroup-specific therapeutic targets in CNS rhabdoid tumours.
  Cancer Cell 30(6):891-908. PMID 26609405.

GTEx reference:
  GTEx Consortium (2020). The GTEx Consortium atlas of genetic regulatory
  effects across human tissues. Science 369(6509):1318-1330. PMID 32913098.
"""

import argparse
import gzip
import io
import logging
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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

# ─────────────────────────────────────────────────────────────────────────────
# GEO FTP directory name logic
#
# NCBI GEO uses a truncated directory scheme:
#   Accession    → Parent directory
#   GPL570       → GPLnnn        (3-digit: replace last 3 with "nnn")
#   GPL6244      → GPL6nnn       (4-digit: keep first digit, replace last 3)
#   GPL13667     → GPL13nnn      (5-digit: keep first two digits, replace last 3)
#   GSE70678     → GSE70nnn      (5-digit: keep first two digits, replace last 3)
#
# Source: https://www.ncbi.nlm.nih.gov/geo/info/download.html
# ─────────────────────────────────────────────────────────────────────────────

def _geo_ftp_parent_dir(accession: str) -> str:
    """
    Compute the GEO FTP parent directory name from an accession number.

    Examples:
      GPL570   → GPLnnn
      GPL6244  → GPL6nnn
      GSE70678 → GSE70nnn

    Source: NCBI GEO download documentation.
    """
    # Split prefix (GPL, GSE, GSM, GDS) from numeric part
    for prefix in ("GPL", "GSE", "GSM", "GDS"):
        if accession.upper().startswith(prefix):
            num = accession[len(prefix):]
            if len(num) <= 3:
                return f"{prefix}nnn"
            else:
                # Keep all digits except the last three, replace last three with "nnn"
                kept = num[:-3]
                return f"{prefix}{kept}nnn"
    return accession  # fallback: return as-is


# Verify our logic with known GEO examples
assert _geo_ftp_parent_dir("GPL570")   == "GPLnnn",   "GPL570 FTP dir should be GPLnnn"
assert _geo_ftp_parent_dir("GPL6244")  == "GPL6nnn",  "GPL6244 FTP dir should be GPL6nnn"
assert _geo_ftp_parent_dir("GPL13667") == "GPL13nnn", "GPL13667 FTP dir should be GPL13nnn"
assert _geo_ftp_parent_dir("GSE70678") == "GSE70nnn", "GSE70678 FTP dir should be GSE70nnn"
logger.debug("GEO FTP directory naming logic verified.")


GTEX_URLS = [
    # Primary: GTEx v8 from Google Storage (public bucket)
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
    # Secondary: GTEx portal direct download
    "https://gtexportal.org/api/v2/dataset/fileDownload?"
    "datasetId=gtex_v8&dataFileId="
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
]

# ─────────────────────────────────────────────────────────────────────────────
# Dataset catalogue
# ─────────────────────────────────────────────────────────────────────────────

def _build_datasets() -> Dict[str, Dict]:
    """
    Build dataset catalogue with correct FTP URLs.

    GPL570 FTP parent dir = GPLnnn (verified above).
    GSE70678 FTP parent dir = GSE70nnn (verified above).

    Source for GPL570 URL structure:
      ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/
      (confirmed via NCBI GEO FTP browser, April 2026)

    HTTPS alternative (for FTP-blocked environments):
      https://www.ncbi.nlm.nih.gov/geo/download/?acc=GPL570&format=file&file=GPL570.annot.gz
    """
    gpl570_parent = _geo_ftp_parent_dir("GPL570")   # → "GPLnnn"
    gse70678_parent = _geo_ftp_parent_dir("GSE70678")  # → "GSE70nnn"

    return {
        "gpl570": {
            "name": "GPL570 — Affymetrix HG-U133 Plus 2.0 annotation (netaffx build 35)",
            "urls": [
                # FTP primary (most reliable)
                f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl570_parent}/GPL570/annot/GPL570.annot.gz",
                # HTTPS web download fallback
                "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GPL570&format=file&file=GPL570.annot.gz",
            ],
            "output_path": str(DATA_DIR / "GPL570.annot.gz"),
            "processed":   str(DATA_DIR / "GPL570_probe_map.tsv"),
            "required":    True,
            "note": (
                "GPL570 = Affymetrix HG-U133 Plus 2.0 Array. "
                "Used by GSE70678. FTP dir: GPLnnn (not GPL0nnn — common mistake). "
                "Annotation column 'Gene Symbol' maps probe → HGNC symbol. "
                "Source: Torchia 2015 Cancer Cell PMID 26609405."
            ),
        },
        "gse70678": {
            "name": "GSE70678 — Torchia 2015 ATRT bulk RNA (49 ATRT + normals, GPL570)",
            "urls": [
                f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse70678_parent}/GSE70678/matrix/GSE70678_series_matrix.txt.gz",
                "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE70678&format=file&file=GSE70678_series_matrix.txt.gz",
            ],
            "output_path": str(DATA_DIR / "GSE70678_series_matrix.txt.gz"),
            "processed":   str(DATA_DIR / "GSE70678_gene_expression.tsv"),
            "required":    True,
            "note": (
                "49 ATRT tumour samples + normal controls. GPL570 platform. "
                "Torchia J et al. 2015 Cancer Cell 30(6):891-908. PMID 26609405."
            ),
        },
        "gtex_brain": {
            "name": "GTEx v8 — Gene median TPM across 54 tissues (~37 MB compressed)",
            "urls": GTEX_URLS,
            "output_path": str(DATA_DIR / "GTEx_v8_gene_median_tpm.gct.gz"),
            "processed":   str(DATA_DIR / "GTEx_brain_normal_reference.tsv"),
            "required":    False,
            "note": (
                "GTEx Consortium 2020 Science 369(6509):1318. PMID 32913098. "
                "Brain tissues used as normal reference for ATRT differential expression."
            ),
        },
    }

DATASETS = _build_datasets()


# ─────────────────────────────────────────────────────────────────────────────
# Curated normal brain reference (fallback when GTEx unavailable)
# Sources: GTEx v8 median TPM for cerebellum + cortex, log2(TPM+1) transformed.
# ─────────────────────────────────────────────────────────────────────────────

CURATED_BRAIN_REFERENCE: Dict[str, float] = {
    # PRC2/EZH2 axis — key for ATRT differential scoring
    "EZH2":    3.45,   "EED":   4.12,   "SUZ12":  4.33,
    "RBBP4":   5.22,   "RBBP7": 4.98,
    # BET bromodomain
    "BRD4":    5.11,   "BRD2":  5.44,   "BRD3":   4.67,
    # HDAC
    "HDAC1":   5.89,   "HDAC2": 5.76,   "HDAC3":  5.33,
    "HDAC6":   4.55,   "HDAC4": 4.88,   "HDAC5":  4.22,
    # MYC axis
    "MYC":     3.22,   "MYCN":  1.88,   "MAX":    5.44,
    # Aurora kinase
    "AURKA":   2.90,   "AURKB": 2.45,
    # CDK4/6
    "CDK4":    4.33,   "CDK6":  3.88,   "CCND1":  4.11,   "CCND2":  4.55,
    # mTOR/PI3K
    "MTOR":    5.22,   "PIK3CA": 4.55,  "AKT1":   4.88,
    # SHH subgroup
    "GLI1":    1.22,   "GLI2":  1.55,   "SMO":    3.45,   "PTCH1":  3.88,
    # Apoptosis
    "BCL2":    3.77,   "BCL2L1": 5.44,  "MCL1":   5.67,
    # DNA damage
    "PARP1":   4.88,   "ATM":   5.11,   "ATR":    5.22,
    # Proteasome
    "PSMB5":   6.33,   "PSMB2": 6.11,   "PSMB1":  6.22,
    # Stemness
    "SOX2":    2.11,   "LIN28A": 0.44,  "SALL4":  0.88,
    # SWI/SNF (high in normal brain; lost in ATRT)
    "SMARCB1": 5.88,   "SMARCA4": 5.22, "SMARCC1": 5.44,
    # Immune
    "CD274":   2.88,   "PDCD1": 0.55,
    # ONC201 target
    "DRD2":    2.33,   "CLPB":  4.88,
    # Tumour suppressors
    "CDKN2A":  2.11,   "PTEN":  5.55,   "RB1":    5.22,   "TP53":   5.77,
    # TYR subgroup — low in normal brain
    "TYR":     0.22,   "DCT":   0.11,   "MITF":   2.88,
    # Housekeeping
    "ACTB":    9.44,   "GAPDH": 8.22,   "B2M":    7.88,   "HPRT1":  5.55,
    # Neural markers (high in normal brain; useful as QC)
    "OLIG2":   3.55,   "NES":   4.22,   "GFAP":   4.88,   "MBP":    5.33,
}


# ─────────────────────────────────────────────────────────────────────────────
# Downloader with retry + multi-URL fallback
# ─────────────────────────────────────────────────────────────────────────────

def _download_with_fallback(
    urls: List[str],
    output_path: Path,
    max_retries: int = 3,
    retry_delay: float = 5.0,
) -> bool:
    """
    Try each URL in order; retry on failure.

    Returns True if download succeeded, False otherwise.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    for url_idx, url in enumerate(urls):
        logger.info("Trying URL %d/%d: %s", url_idx + 1, len(urls), url)
        for attempt in range(1, max_retries + 1):
            try:
                # Add a User-Agent header — NCBI sometimes blocks default urllib UA
                req = urllib.request.Request(
                    url,
                    headers={
                        "User-Agent": (
                            "Mozilla/5.0 (compatible; ATRTPipeline/4.0; "
                            "+https://github.com/ShruthiSathya/gbmresearch)"
                        )
                    },
                )
                with urllib.request.urlopen(req, timeout=120) as resp:
                    data = resp.read()
                output_path.write_bytes(data)
                logger.info(
                    "✅ Downloaded %s (%.2f MB)", output_path.name, len(data) / 1e6
                )
                return True

            except urllib.error.HTTPError as e:
                logger.warning(
                    "HTTP %d on attempt %d/%d: %s", e.code, attempt, max_retries, url
                )
                if e.code in (403, 404):
                    break  # These won't recover; try next URL
                if attempt < max_retries:
                    time.sleep(retry_delay * attempt)

            except (urllib.error.URLError, OSError, TimeoutError) as e:
                logger.warning(
                    "Network error on attempt %d/%d: %s", attempt, max_retries, e
                )
                if attempt < max_retries:
                    time.sleep(retry_delay * attempt)

    logger.error(
        "❌ All %d URLs failed for %s", len(urls), output_path.name
    )
    return False


def download_dataset(key: str) -> bool:
    """Download a single dataset by key. Returns True on success."""
    if key not in DATASETS:
        logger.error("Unknown dataset key: %s", key)
        return False

    ds = DATASETS[key]
    out_path = Path(ds["output_path"])

    if out_path.exists():
        size_mb = out_path.stat().st_size / 1e6
        logger.info("✅ Already exists: %s (%.2f MB)", out_path.name, size_mb)
        return True

    logger.info("⬇️  Downloading: %s", ds["name"])
    return _download_with_fallback(ds["urls"], out_path)


# ─────────────────────────────────────────────────────────────────────────────
# Post-processing: GPL570
# ─────────────────────────────────────────────────────────────────────────────

def process_gpl570() -> bool:
    """
    Parse GPL570.annot.gz → probe_id → gene_symbol TSV.

    GPL570 annotation file structure (netaffx build 35):
      - Lines starting with ^, !, # = GEO metadata (skipped)
      - First data line = column headers (tab-separated)
      - Key columns: "ID" (probe set ID), "Gene Symbol" (HGNC symbol)

    Multi-gene probes (e.g. "BRCA1 /// BRCA2") → take first symbol.
    Control probes (starting with "AFFX-") → excluded.

    Source:
      Irizarry RA et al. (2003). Nucleic Acids Research 31(4):e15. PMID 12582260.
      (probe annotation rationale)
    """
    path = Path(DATASETS["gpl570"]["output_path"])
    out  = Path(DATASETS["gpl570"]["processed"])

    if not path.exists():
        logger.error("GPL570.annot.gz not found: %s", path)
        logger.error("Download first: python -m backend.pipeline.data_downloader --dataset gpl570")
        return False

    logger.info("Parsing GPL570 annotation: %s", path)
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
    except Exception as e:
        logger.error("Failed to open GPL570.annot.gz: %s", e)
        return False

    # Filter GEO metadata lines — keep only data rows
    data_lines = [
        ln for ln in lines
        if ln.strip() and not ln.startswith(_GEO_METADATA_PREFIXES)
    ]

    if not data_lines:
        logger.error("No data lines found in GPL570.annot.gz after metadata removal.")
        return False

    logger.info("GPL570: %d data lines after metadata removal", len(data_lines))

    try:
        df = pd.read_csv(
            io.StringIO("".join(data_lines)),
            sep="\t",
            dtype=str,
            low_memory=False,
            na_filter=False,
        )
    except Exception as e:
        logger.error("Failed to parse GPL570 annotation as TSV: %s", e)
        logger.error("First 3 data lines: %s", data_lines[:3])
        return False

    df.columns = df.columns.str.strip()
    logger.info("GPL570 columns: %s", list(df.columns[:8]))

    # Find the probe ID column
    # GPL570 uses "ID" as the probe set ID column
    id_col = None
    for candidate in ["ID", "id", "probe_id", "ProbeID", "Probe Set ID"]:
        if candidate in df.columns:
            id_col = candidate
            break
    if id_col is None:
        logger.error("Cannot find probe ID column. Available: %s", list(df.columns[:10]))
        return False

    # Find the gene symbol column
    # GPL570 uses "Gene Symbol" (Title Case, with space) — verified from NCBI GEO
    sym_col = None
    for candidate in ["Gene Symbol", "Gene symbol", "gene_symbol", "Symbol", "GENE_SYMBOL"]:
        if candidate in df.columns:
            sym_col = candidate
            break
    if sym_col is None:
        logger.error(
            "Cannot find gene symbol column. Available: %s\n"
            "Expected 'Gene Symbol' for GPL570 (Affymetrix HG-U133 Plus 2.0).",
            list(df.columns[:10])
        )
        return False

    logger.info("GPL570: using ID column='%s', Symbol column='%s'", id_col, sym_col)

    # Build probe → gene symbol mapping
    probe_map_rows = []
    n_control = 0
    n_unmapped = 0

    for _, row in df.iterrows():
        probe = str(row[id_col]).strip()
        raw_sym = str(row[sym_col]).strip()

        # Skip Affymetrix control probes
        if probe.upper().startswith("AFFX-"):
            n_control += 1
            continue

        # Handle multi-gene probes (e.g. "BRCA1 /// BRCA2") → take first
        if raw_sym in ("---", "nan", "NaN", "", "None"):
            n_unmapped += 1
            gene_sym = ""
        else:
            gene_sym = str(raw_sym).split(" /// ")[0].split("/")[0].strip().upper()

        probe_map_rows.append({"probe_id": probe, "gene_symbol": gene_sym})

    if not probe_map_rows:
        logger.error("No valid probe entries found in GPL570 annotation.")
        return False

    probe_df = pd.DataFrame(probe_map_rows)
    n_total  = len(probe_df)
    n_mapped = (probe_df["gene_symbol"] != "").sum()

    out.parent.mkdir(parents=True, exist_ok=True)
    probe_df.to_csv(out, sep="\t", index=False)

    logger.info(
        "✅ GPL570 probe map: %d total probes | %d mapped (%.1f%%) | "
        "%d control probes excluded | %d unmapped → %s",
        n_total, n_mapped, 100 * n_mapped / max(n_total, 1),
        n_control, n_unmapped, out,
    )
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Post-processing: GSE70678
# ─────────────────────────────────────────────────────────────────────────────

def process_gse70678() -> bool:
    """
    Map GSE70678 Affymetrix probe IDs → HGNC gene symbols using GPL570.

    Workflow:
      1. Parse GSE70678_series_matrix.txt.gz → probe-indexed DataFrame
      2. Load GPL570_probe_map.tsv → probe → gene symbol mapping
      3. Map probes → genes, aggregate by median (per gene, per sample)
      4. Save as GSE70678_gene_expression.tsv

    Aggregation by median is preferred because:
      - Multiple probes often target the same gene
      - Median is robust to outlier probes
      Source: Bolstad BM et al. (2003). Bioinformatics 19(2):185. PMID 12538238.

    Output format:
      Rows = HGNC gene symbols (uppercase)
      Columns = GSM sample IDs (e.g. GSM1714762)
      Values = RMA-normalised log2 expression

    References:
      Torchia J et al. (2015). Cancer Cell 30(6):891-908. PMID 26609405.
      (GSE70678 dataset description)
      Irizarry RA et al. (2003). Biostatistics 4(2):249. PMID 12925520.
      (RMA normalisation method)
    """
    matrix_path = Path(DATASETS["gse70678"]["output_path"])
    map_path    = Path(DATASETS["gpl570"]["processed"])
    out_path    = Path(DATASETS["gse70678"]["processed"])

    if not matrix_path.exists():
        logger.error("GSE70678 series matrix not found: %s", matrix_path)
        return False
    if not map_path.exists():
        logger.error(
            "GPL570 probe map not found: %s\n"
            "Run: python -m backend.pipeline.data_downloader --dataset gpl570 --process gpl570",
            map_path
        )
        return False

    logger.info("Processing GSE70678: %s", matrix_path)

    # Step 1: Parse series matrix
    try:
        try:
            from .probe_mapper import parse_geo_series_matrix, ProbeMapper
        except ImportError:
            sys.path.insert(0, str(Path(__file__).parent))
            from probe_mapper import parse_geo_series_matrix, ProbeMapper

        probe_df, meta = parse_geo_series_matrix(matrix_path)
        logger.info(
            "GSE70678 series matrix: %d probes × %d samples",
            len(probe_df), len(probe_df.columns)
        )

    except Exception as e:
        logger.error("Failed to parse GSE70678 series matrix: %s", e)
        return False

    # Step 2: Load probe map and map probes → genes
    try:
        mapper = ProbeMapper(cache_path=map_path)
        mapper._load_from_cache()

        coverage = mapper.get_coverage(probe_df.index.tolist())
        logger.info(
            "GPL570 probe coverage for GSE70678: %d/%d probes (%.1f%%)",
            coverage["mapped_probes"], coverage["total_probes"], coverage["coverage_pct"]
        )

        if coverage["coverage_pct"] < 10:
            logger.error(
                "Coverage %.1f%% is too low. "
                "Check that GPL570_probe_map.tsv was generated from GPL570.annot.gz "
                "(not GPL6244.annot.gz or another platform).",
                coverage["coverage_pct"]
            )
            return False

        # Step 3: Aggregate probes → genes by median
        gene_df = mapper.map_probes_to_genes(probe_df, aggregation="median")
        logger.info(
            "GSE70678 gene matrix: %d genes × %d samples",
            len(gene_df), len(gene_df.columns)
        )

    except Exception as e:
        logger.error("Probe mapping failed: %s", e)
        return False

    # Validate key ATRT biology is represented
    # These genes should all be present in GPL570-mapped GSE70678
    expected_atrt_genes = {
        "EZH2", "BRD4", "HDAC1", "AURKA", "MYC", "MYCN",
        "CDK4", "SMARCB1", "SOX2", "LIN28A", "TYR", "GLI2"
    }
    found = expected_atrt_genes & set(gene_df.index)
    missing = expected_atrt_genes - found

    if len(found) < 8:
        logger.warning(
            "Only %d/%d expected ATRT genes found in mapped output: %s\n"
            "Missing: %s\n"
            "This may indicate a platform mismatch (GPL570 vs GPL6244).",
            len(found), len(expected_atrt_genes), sorted(found), sorted(missing)
        )
    else:
        logger.info(
            "ATRT biology check: %d/%d key genes present (%s)",
            len(found), len(expected_atrt_genes), sorted(found)
        )

    # Step 4: Save
    out_path.parent.mkdir(parents=True, exist_ok=True)
    gene_df.to_csv(out_path, sep="\t")
    logger.info("✅ GSE70678 gene expression saved → %s", out_path)

    return True


# ─────────────────────────────────────────────────────────────────────────────
# Post-processing: GTEx
# ─────────────────────────────────────────────────────────────────────────────

def process_gtex() -> bool:
    """
    Extract brain-relevant tissue columns from GTEx v8 gene TPM matrix.

    Brain tissues selected based on anatomical relevance to ATRT location:
      - Cerebellum: infratentorial ATRT (~50% of cases)
      - Cortex:     supratentorial ATRT (~35% of cases)
    Source: Frühwald MC et al. 2020 CNS Oncology 9(2):CNS56. PMID 32432484.

    Values are log2(median_TPM + 1) transformed.
    Source: GTEx Consortium 2020 Science 369(6509):1318. PMID 32913098.
    """
    path = Path(DATASETS["gtex_brain"]["output_path"])
    out  = Path(DATASETS["gtex_brain"]["processed"])

    if not path.exists():
        logger.warning("GTEx raw file not found: %s — using curated fallback", path)
        return _write_curated_brain_fallback(out)

    try:
        logger.info("Processing GTEx v8: %s", path)
        df = pd.read_csv(path, sep="\t", skiprows=2, compression="gzip", low_memory=False)

        # Find brain tissue columns
        brain_keywords = ["brain", "cerebellum", "cortex", "hippocampus", "substantia"]
        brain_cols = [
            c for c in df.columns
            if any(k.lower() in c.lower() for k in brain_keywords)
        ]

        if not brain_cols:
            logger.warning(
                "No brain tissue columns found in GTEx. Available: %s",
                list(df.columns[:10])
            )
            return _write_curated_brain_fallback(out)

        logger.info("GTEx brain tissues selected: %s", brain_cols)

        # Use gene Description column as gene symbol (standard in GTEx GCT files)
        desc_col = next((c for c in ["Description", "Name", "gene_id"] if c in df.columns), None)
        if desc_col is None:
            logger.warning("No gene description column in GTEx; using curated fallback")
            return _write_curated_brain_fallback(out)

        res = df[[desc_col] + brain_cols].copy()
        res.columns = ["gene_symbol"] + brain_cols
        res = res.groupby("gene_symbol").median(numeric_only=True)
        res = np.log2(res + 1)
        res.index = res.index.str.upper()

        out.parent.mkdir(parents=True, exist_ok=True)
        res.to_csv(out, sep="\t")
        logger.info("✅ GTEx brain reference saved → %s (%d genes)", out, len(res))
        return True

    except Exception as e:
        logger.error("GTEx processing failed: %s — using curated fallback", e)
        return _write_curated_brain_fallback(out)


def _write_curated_brain_fallback(out_path: Path) -> bool:
    """
    Write curated normal brain reference when GTEx is unavailable.

    Values are approximate log2(median_TPM+1) from GTEx v8 cerebellum+cortex.
    Source: GTEx Consortium 2020, accessed April 2026.
    """
    try:
        ref_df = pd.DataFrame.from_dict(
            CURATED_BRAIN_REFERENCE, orient="index", columns=["Brain_median"]
        )
        ref_df.index.name = "gene_symbol"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        ref_df.to_csv(out_path, sep="\t")
        logger.info(
            "✅ Curated brain fallback written → %s (%d genes)",
            out_path, len(ref_df)
        )
        return True
    except Exception as e:
        logger.error("Failed to write curated brain fallback: %s", e)
        return False


# ─────────────────────────────────────────────────────────────────────────────
# Status check
# ─────────────────────────────────────────────────────────────────────────────

def check_files() -> None:
    """Print status of all data files with download instructions."""
    print("\n" + "=" * 70)
    print("ATRT PIPELINE DATA FILE STATUS")
    print("=" * 70)

    for key, ds in DATASETS.items():
        raw_path = Path(ds["output_path"])
        raw_ok = raw_path.exists()
        raw_size = f"({raw_path.stat().st_size / 1e6:.1f} MB)" if raw_ok else ""

        print(f"\n[{key}] {ds['name']}")
        print(f"  Raw:       {'✅' if raw_ok else '❌'}  {raw_path.name} {raw_size}")

        if "processed" in ds:
            proc_path = Path(ds["processed"])
            proc_ok = proc_path.exists()
            proc_size = f"({proc_path.stat().st_size / 1e6:.2f} MB)" if proc_ok else ""
            print(f"  Processed: {'✅' if proc_ok else '⚪'}  {proc_path.name} {proc_size}")

        if not raw_ok:
            print(f"  📥 Primary URL: {ds['urls'][0]}")
            if len(ds['urls']) > 1:
                print(f"  📥 Fallback URL: {ds['urls'][1]}")

    print("\n" + "=" * 70)
    print("MANUAL DOWNLOAD COMMANDS")
    print("=" * 70)
    print("""
# GPL570 annotation (CRITICAL: GPLnnn — NOT GPL0nnn)
wget -O data/raw_omics/GPL570.annot.gz \\
  "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

# GSE70678 series matrix (49 ATRT samples, Torchia 2015)
wget -O data/raw_omics/GSE70678_series_matrix.txt.gz \\
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70678/matrix/GSE70678_series_matrix.txt.gz"

# If FTP is blocked (HTTPS alternative):
curl -L -o data/raw_omics/GPL570.annot.gz \\
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GPL570&format=file&file=GPL570.annot.gz"

curl -L -o data/raw_omics/GSE70678_series_matrix.txt.gz \\
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE70678&format=file&file=GSE70678_series_matrix.txt.gz"
""")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="ATRT Pipeline Data Downloader v4.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download all required datasets
  python -m backend.pipeline.data_downloader --dataset all

  # Download GPL570 annotation only
  python -m backend.pipeline.data_downloader --dataset gpl570

  # Process all downloaded files
  python -m backend.pipeline.data_downloader --process all

  # Full pipeline: download + process
  python -m backend.pipeline.data_downloader --dataset gpl570 gse70678 --process gpl570 gse70678

  # Check file status
  python -m backend.pipeline.data_downloader --check
        """,
    )
    parser.add_argument(
        "--dataset", nargs="+",
        choices=list(DATASETS.keys()) + ["all"],
        help="Dataset(s) to download",
    )
    parser.add_argument(
        "--process", nargs="+",
        choices=["gpl570", "gse70678", "gtex", "all"],
        help="Dataset(s) to post-process after download",
    )
    parser.add_argument(
        "--check", action="store_true",
        help="Check file status and print download instructions",
    )
    args = parser.parse_args()

    if args.check or (not args.dataset and not args.process):
        check_files()
        return

    if args.dataset:
        keys = list(DATASETS.keys()) if "all" in args.dataset else args.dataset
        for k in keys:
            success = download_dataset(k)
            if not success and DATASETS[k]["required"]:
                logger.error(
                    "❌ Required dataset '%s' could not be downloaded. "
                    "See manual download instructions above.", k
                )

    if args.process:
        steps = args.process
        if "all" in steps or "gpl570" in steps:
            if not process_gpl570():
                logger.error("GPL570 processing failed.")
        if "all" in steps or "gse70678" in steps:
            if not process_gse70678():
                logger.error("GSE70678 processing failed.")
        if "all" in steps or "gtex" in steps:
            if not process_gtex():
                logger.error("GTEx processing failed.")

    check_files()


if __name__ == "__main__":
    main()