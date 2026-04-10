"""
data_downloader.py  (v2.0)
==========================
Automated data acquisition for the ATRT Drug Repurposing Pipeline.

KEY CHANGES FROM v1.0
----------------------
1. GPL6244 probe annotation download (required to map probe IDs → gene symbols)
2. GTEx v8 normal brain download (required for differential expression baseline)
3. Cleaner status reporting with per-file size and data-readiness check

RUN ORDER
----------
Step 1 (required)  : --dataset depmap         → manual (portal auth needed)
Step 2 (required)  : --dataset gpl6244         → auto download (~5 MB)
Step 3 (required)  : --dataset gse70678        → auto download (~50 MB)
Step 4 (required)  : --dataset gtex_brain      → auto download (~300 MB)
Step 5 (optional)  : --dataset gse106982       → auto download (~30 MB)
Step 6 (optional)  : cbtn ATRT                 → manual (DCC access needed)

WHY GTEx FOR NORMAL BRAIN?
----------------------------
GSE70678 only contains 49 ATRT tumour samples — no matched normal brain controls.
The original Torchia 2015 paper compared against normal brain from a separate
GEO submission. GTEx v8 provides the highest-quality human brain expression
reference (n=209 cerebellum, n=255 cortex donors) and is the current community
standard for normal tissue expression.

GTEx data is at gene level (TPM) — compatible with our probe-mapped ATRT data
after log2(TPM + 1) transformation.

Reference: GTEx Consortium (2020) Science 369(6509). PMID 32913098.
"""

import argparse
import gzip
import logging
import os
import shutil
import time
import urllib.request
from pathlib import Path
from typing import Dict, Optional, Tuple

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

    # ── GPL6244 probe annotation ──────────────────────────────────────────────
    "gpl6244": {
        "name":        "GPL570 — Affymetrix HG-U133 Plus 2.0 probe annotation",
        "url":         "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL0nnn/GPL570/annot/GPL570.annot.gz",
        "output_path": "data/raw_omics/GPL570.annot.gz",
        "processed":   "data/raw_omics/GPL570_probe_map.tsv",
        "compressed":  False,
        "size_mb":     5,
        "required":    True,
        "pmid":        None,
        "note":        "Corrected for GSE70678 ATRT samples.",
    },

    # ── GSE70678 ATRT bulk RNA-seq ────────────────────────────────────────────
    "gse70678": {
        "name":        "GSE70678 — Torchia 2015 ATRT bulk RNA-seq (49 tumour samples)",
        "url":         (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/"
            "GSE70nnn/GSE70678/matrix/GSE70678_series_matrix.txt.gz"
        ),
        "output_path": "data/raw_omics/GSE70678_series_matrix.txt.gz",
        "processed":   "data/raw_omics/GSE70678_gene_expression.tsv",
        "compressed":  False,   # Keep as .gz; tissue_expression.py handles it
        "size_mb":     50,
        "required":    True,
        "pmid":        "26609405",
        "note": (
            "Contains Affymetrix probe IDs — must be processed by probe_mapper.py "
            "(ProbeMapper class) before use. Run process_gse70678() after download."
        ),
    },

    # ── GTEx v8 brain normal reference ────────────────────────────────────────
    "gtex_brain": {
        "name":        "GTEx v8 — Gene median TPM (all tissues, ~300 MB)",
        "url":         (
            "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
            "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
        ),
        "output_path": "data/raw_omics/GTEx_v8_gene_median_tpm.gct.gz",
        "processed":   "data/raw_omics/GTEx_brain_normal_reference.tsv",
        "compressed":  False,
        "size_mb":     300,
        "required":    True,
        "pmid":        "32913098",
        "note": (
            "GTEx v8 median TPM per gene per tissue. "
            "After download run process_gtex() to extract brain tissues and "
            "save a compact reference file (only Brain - Cerebellum + Brain - Cortex). "
            "Full file is 300 MB but the processed output is <5 MB."
        ),
    },

    # ── GSE106982 methylation subgroups (optional) ────────────────────────────
    "gse106982": {
        "name":        "GSE106982 — Johann 2016 ATRT methylation subgroups (optional)",
        "url":         (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/"
            "GSE106nnn/GSE106982/matrix/GSE106982_series_matrix.txt.gz"
        ),
        "output_path": "data/raw_omics/GSE106982_series_matrix.txt.gz",
        "processed":   None,
        "compressed":  False,
        "size_mb":     30,
        "required":    False,
        "pmid":        "26923874",
        "note":        "Methylation-based subgroup labels for 150 ATRT. Optional.",
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Generic download helpers
# ─────────────────────────────────────────────────────────────────────────────

def _progress_hook(count: int, block_size: int, total_size: int) -> None:
    if total_size > 0:
        pct = min(int(count * block_size * 100 / total_size), 100)
        print(f"\r  Downloading: {pct:3d}%", end="", flush=True)


def download_file(
    url: str,
    output_path: str,
    decompress_to: Optional[str] = None,
    max_retries: int = 3,
    timeout: int = 120,
) -> bool:
    """
    Download a file from URL with retry logic.

    Parameters
    ----------
    url          : Remote URL
    output_path  : Local destination path (may be .gz)
    decompress_to: If set, decompress .gz to this path after downloading.
                   Set to None to keep the compressed file.
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")

    for attempt in range(1, max_retries + 1):
        try:
            logger.info("  Attempt %d/%d — %s", attempt, max_retries, url)
            urllib.request.urlretrieve(url, tmp, reporthook=_progress_hook)
            print()  # newline after progress bar
            tmp.rename(out)

            if decompress_to:
                _decompress(out, Path(decompress_to))

            size_mb = out.stat().st_size / 1_048_576
            logger.info("  ✅ %s (%.1f MB)", out, size_mb)
            return True

        except Exception as exc:
            logger.warning("  Attempt %d failed: %s", attempt, exc)
            if tmp.exists():
                tmp.unlink()
            if attempt < max_retries:
                time.sleep(5 * attempt)

    logger.error("  ❌ All %d attempts failed for %s", max_retries, url)
    return False


def _decompress(gz_path: Path, out_path: Path) -> None:
    """Decompress a .gz file to out_path."""
    logger.info("  Decompressing → %s", out_path)
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


# ─────────────────────────────────────────────────────────────────────────────
# Post-processing helpers
# ─────────────────────────────────────────────────────────────────────────────

def process_gse70678() -> bool:
    """
    Convert GSE70678 from probe IDs to gene symbols and save as TSV.

    Requires:
      data/raw_omics/GSE70678_series_matrix.txt.gz   (downloaded above)
      data/raw_omics/GPL6244.annot.gz                (downloaded above)

    Produces:
      data/raw_omics/GSE70678_gene_expression.tsv
        Rows  = HGNC gene symbols
        Cols  = sample GSM IDs (49 ATRT samples)
        Values = log2-scale Affymetrix expression (probe-median aggregated)
    """
    try:
        import sys
        sys.path.insert(0, str(Path(__file__).parent))
        from probe_mapper import ProbeMapper, parse_geo_series_matrix

        series_path = Path("data/raw_omics/GSE70678_series_matrix.txt.gz")
        out_path    = Path("data/raw_omics/GSE70678_gene_expression.tsv")

        if not series_path.exists():
            logger.error(
                "GSE70678 series matrix not found at %s. "
                "Run: python -m backend.pipeline.data_downloader --dataset gse70678",
                series_path,
            )
            return False

        # 1. Parse the series matrix
        logger.info("Processing GSE70678 series matrix → gene-level TSV...")
        probe_df, meta = parse_geo_series_matrix(series_path)

        # 2. Map probes → gene symbols using GPL6244 annotation
        mapper = ProbeMapper()
        mapper.load()
        coverage = mapper.get_coverage(probe_df.index.tolist())
        logger.info(
            "Probe coverage: %d/%d mapped (%.1f%%)",
            coverage["mapped_probes"], coverage["total_probes"],
            coverage["coverage_pct"],
        )

        gene_df = mapper.map_probes_to_genes(probe_df, aggregation="median")

        # 3. Save gene-level expression
        out_path.parent.mkdir(parents=True, exist_ok=True)
        gene_df.to_csv(out_path, sep="\t")
        logger.info(
            "✅ Saved gene expression matrix → %s (%d genes × %d samples)",
            out_path, len(gene_df), len(gene_df.columns),
        )

        # 4. Save sample metadata summary
        meta_path = Path("data/raw_omics/GSE70678_sample_metadata.tsv")
        if "Sample_title" in meta:
            import pandas as pd
            titles = meta["Sample_title"]
            gsm_ids = meta.get("Sample_geo_accession", [""] * len(titles))
            pd.DataFrame({"gsm_id": gsm_ids, "title": titles}).to_csv(
                meta_path, sep="\t", index=False
            )
            logger.info("Sample metadata saved → %s", meta_path)

        return True

    except Exception as exc:
        logger.error("process_gse70678() failed: %s", exc)
        import traceback
        traceback.print_exc()
        return False


def process_gtex() -> bool:
    """
    Extract brain tissue columns from the GTEx v8 median TPM file and save
    a compact reference TSV containing only the brain tissues needed for
    differential expression vs ATRT.

    Requires:
      data/raw_omics/GTEx_v8_gene_median_tpm.gct.gz   (~300 MB)

    Produces:
      data/raw_omics/GTEx_brain_normal_reference.tsv
        Rows  = HGNC gene symbols
        Cols  = selected brain tissue names
        Values = log2(median_TPM + 1)

    Brain tissues extracted:
      Brain - Cerebellum           (~50% of ATRT tumours are infratentorial)
      Brain - Cerebellar Hemisphere
      Brain - Cortex               (~35% of ATRT tumours are supratentorial)
      Brain - Frontal Cortex (BA9)
    """
    import gzip as _gz
    import pandas as pd

    gz_path  = Path("data/raw_omics/GTEx_v8_gene_median_tpm.gct.gz")
    out_path = Path("data/raw_omics/GTEx_brain_normal_reference.tsv")

    if not gz_path.exists():
        logger.error(
            "GTEx file not found at %s. "
            "Run: python -m backend.pipeline.data_downloader --dataset gtex_brain",
            gz_path,
        )
        return False

    logger.info("Processing GTEx v8 (extracting brain tissues)...")

    try:
        # GCT format:
        #   Line 1: "#1.2"
        #   Line 2: "n_genes\tn_samples"
        #   Line 3: "Name\tDescription\tsample1\tsample2..."
        #   Lines 4+: data
        with _gz.open(gz_path, "rt", encoding="utf-8") as fh:
            lines = fh.readlines()

        # Skip GCT header lines
        header_idx = 0
        for i, ln in enumerate(lines):
            if ln.startswith("Name\t"):
                header_idx = i
                break

        if header_idx == 0:
            raise ValueError("Could not find GCT column header in GTEx file")

        logger.info("GTEx GCT header found at line %d", header_idx)

        # Parse in chunks to manage memory — file is ~300 MB uncompressed
        gct_text = "".join(lines[header_idx:])
        df = pd.read_csv(
            pd.io.common.StringIO(gct_text),
            sep="\t",
            low_memory=False,
        )

        # Rename columns
        df = df.rename(columns={"Name": "ensembl_id", "Description": "gene_symbol"})

        # The GTEx GCT gene column contains Ensembl IDs with version suffix,
        # e.g., "ENSG00000111276.14". We prefer the Description = gene symbol.
        if "gene_symbol" not in df.columns:
            df["gene_symbol"] = df["ensembl_id"].str.split(".").str[0]

        # Extract brain tissues
        brain_tissue_patterns = [
            "Brain - Cerebellum",
            "Brain - Cerebellar Hemisphere",
            "Brain - Cortex",
            "Brain - Frontal Cortex (BA9)",
        ]

        selected_cols = ["gene_symbol"]
        for pattern in brain_tissue_patterns:
            matches = [c for c in df.columns if pattern.lower() in c.lower()]
            selected_cols.extend(matches)
            if not matches:
                logger.warning("GTEx tissue not found: '%s'", pattern)

        brain_df = df[selected_cols].copy()
        brain_df = brain_df.dropna(subset=["gene_symbol"])
        brain_df["gene_symbol"] = (
            brain_df["gene_symbol"].astype(str).str.upper().str.strip()
        )
        brain_df = brain_df[brain_df["gene_symbol"] != ""].copy()

        # Collapse duplicate gene symbols (take median)
        tissue_cols = [c for c in brain_df.columns if c != "gene_symbol"]
        for col in tissue_cols:
            brain_df[col] = pd.to_numeric(brain_df[col], errors="coerce")

        brain_df = brain_df.groupby("gene_symbol")[tissue_cols].median()

        # Log2(TPM + 1) transform for comparability with Affymetrix log2 values
        import numpy as np
        brain_df = brain_df.apply(lambda col: np.log2(col + 1))

        out_path.parent.mkdir(parents=True, exist_ok=True)
        brain_df.to_csv(out_path, sep="\t")
        logger.info(
            "✅ GTEx brain reference saved → %s (%d genes × %d tissues)",
            out_path, len(brain_df), len(tissue_cols),
        )
        logger.info("   Tissues: %s", tissue_cols)
        return True

    except Exception as exc:
        logger.error("process_gtex() failed: %s", exc)
        import traceback
        traceback.print_exc()
        return False


# ─────────────────────────────────────────────────────────────────────────────
# Status check
# ─────────────────────────────────────────────────────────────────────────────

def check_data_status() -> Dict[str, bool]:
    """Print a detailed data readiness report."""
    checks = {
        "GPL6244 annotation (.gz)":            Path("data/raw_omics/GPL6244.annot.gz"),
        "GPL6244 probe map (processed TSV)":   Path("data/raw_omics/GPL6244_probe_map.tsv"),
        "GSE70678 series matrix (.gz)":        Path("data/raw_omics/GSE70678_series_matrix.txt.gz"),
        "GSE70678 gene expression (processed)": Path("data/raw_omics/GSE70678_gene_expression.tsv"),
        "GTEx v8 full file (.gz)":             Path("data/raw_omics/GTEx_v8_gene_median_tpm.gct.gz"),
        "GTEx brain reference (processed)":    Path("data/raw_omics/GTEx_brain_normal_reference.tsv"),
        "DepMap CRISPRGeneEffect.csv":         Path("data/depmap/CRISPRGeneEffect.csv"),
        "DepMap Model.csv":                    Path("data/depmap/Model.csv"),
        "CMap scores JSON":                    Path("data/cmap_query/atrt_cmap_scores.json"),
        "GSE106982 methylation (optional)":    Path("data/raw_omics/GSE106982_series_matrix.txt.gz"),
    }

    print("\n" + "═" * 62)
    print("ATRT Pipeline  —  Data Readiness Report")
    print("═" * 62)

    status: Dict[str, bool] = {}
    for label, path in checks.items():
        present  = path.exists()
        opt_tag  = " [optional]" if "optional" in label else ""
        size_str = ""
        if present and path.is_file():
            mb = path.stat().st_size / 1_048_576
            size_str = f"  ({mb:.1f} MB)"
        icon = "✅" if present else ("⬛" if opt_tag else "❌")
        print(f"  {icon}  {label}{opt_tag}{size_str}")
        status[label] = present

    # Readiness summary
    core_ready = all([
        checks["GPL6244 probe map (processed TSV)"].exists(),
        checks["GSE70678 gene expression (processed)"].exists(),
        checks["GTEx brain reference (processed)"].exists(),
        checks["DepMap CRISPRGeneEffect.csv"].exists(),
        checks["DepMap Model.csv"].exists(),
    ])

    print()
    if core_ready:
        print("  ✅ FULLY READY — all core data present")
    else:
        missing = []
        if not checks["GPL6244 annotation (.gz)"].exists():
            missing.append("GPL6244 (run --dataset gpl6244)")
        if not checks["GSE70678 series matrix (.gz)"].exists():
            missing.append("GSE70678 (run --dataset gse70678)")
        if not checks["GTEx v8 full file (.gz)"].exists():
            missing.append("GTEx brain (run --dataset gtex_brain)")
        if not checks["DepMap CRISPRGeneEffect.csv"].exists():
            missing.append("DepMap CSVs (manual download)")
        print("  ❌ NOT READY — missing:", ", ".join(missing) if missing else "processed files")
        if (
            checks["GPL6244 annotation (.gz)"].exists()
            and not checks["GPL6244 probe map (processed TSV)"].exists()
        ):
            print("     → Run: python -m backend.pipeline.data_downloader --process gpl6244")
        if (
            checks["GSE70678 series matrix (.gz)"].exists()
            and not checks["GSE70678 gene expression (processed)"].exists()
        ):
            print("     → Run: python -m backend.pipeline.data_downloader --process gse70678")
        if (
            checks["GTEx v8 full file (.gz)"].exists()
            and not checks["GTEx brain reference (processed)"].exists()
        ):
            print("     → Run: python -m backend.pipeline.data_downloader --process gtex")

    print("═" * 62 + "\n")
    return status


# ─────────────────────────────────────────────────────────────────────────────
# Manual download instructions
# ─────────────────────────────────────────────────────────────────────────────

DEPMAP_INSTRUCTIONS = """
══════════════════════════════════════════════════════════════
DepMap 24Q4  —  Manual Download Required (portal auth needed)
══════════════════════════════════════════════════════════════
1. Go to: https://depmap.org/portal/download/all/
2. Select release: DepMap Public 24Q4 (Oct 2024)
3. Download: CRISPRGeneEffect.csv  (~500 MB)
             Model.csv             (~5 MB)
4. Place in: data/depmap/

   mkdir -p data/depmap
   mv ~/Downloads/CRISPRGeneEffect.csv data/depmap/
   mv ~/Downloads/Model.csv            data/depmap/

Expected ATRT/rhabdoid lines in Model.csv:
  BT16 (ACH-000725)  BT37 (ACH-000881)
  G401 (ACH-000039)  A204 (ACH-000658)
══════════════════════════════════════════════════════════════
"""


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ATRT Pipeline Data Downloader v2.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--dataset",
        choices=["gpl6244", "gse70678", "gtex_brain", "gse106982", "all"],
        help="Dataset to download",
    )
    parser.add_argument(
        "--process",
        choices=["gpl6244", "gse70678", "gtex", "all"],
        help="Post-process a downloaded dataset (probe mapping / extraction)",
    )
    parser.add_argument(
        "--check",  action="store_true",
        help="Show data readiness report",
    )
    parser.add_argument(
        "--instructions", action="store_true",
        help="Print manual download instructions (DepMap, CMap, CBTN)",
    )
    args = parser.parse_args()

    if args.check or (not args.dataset and not args.process and not args.instructions):
        check_data_status()
        return

    if args.instructions:
        print(DEPMAP_INSTRUCTIONS)
        return

    # ── Downloads ─────────────────────────────────────────────────────────────
    if args.dataset:
        keys = list(DATASETS.keys()) if args.dataset == "all" else [args.dataset]
        for key in keys:
            ds = DATASETS[key]
            logger.info("═" * 55)
            logger.info("Downloading: %s", ds["name"])
            if ds.get("pmid"):
                logger.info("  Reference: PMID %s", ds["pmid"])
            logger.info("  Note: %s", ds.get("note", ""))
            download_file(ds["url"], ds["output_path"])

    # ── Post-processing ───────────────────────────────────────────────────────
    if args.process:
        procs = ["gpl6244", "gse70678", "gtex"] if args.process == "all" else [args.process]

        if "gpl6244" in procs:
            # GPL6244 is parsed by ProbeMapper on first use — just trigger it
            logger.info("Triggering GPL6244 annotation parsing via ProbeMapper...")
            import sys; sys.path.insert(0, str(Path(__file__).parent))
            from probe_mapper import ProbeMapper
            m = ProbeMapper()
            gz = Path("data/raw_omics/GPL6244.annot.gz")
            if gz.exists():
                # Temporarily set cache to the gz location for ProbeMapper
                m.cache_path = Path("data/raw_omics/GPL6244_probe_map.tsv")
                m._download_and_parse()  # Will read from our downloaded .gz
            else:
                m.load()  # Will download itself if needed

        if "gse70678" in procs:
            logger.info("Processing GSE70678 (probe → gene mapping)...")
            ok = process_gse70678()
            if ok:
                logger.info("✅ GSE70678 processing complete")
            else:
                logger.error("❌ GSE70678 processing failed — see above")

        if "gtex" in procs:
            logger.info("Processing GTEx (extracting brain tissues)...")
            ok = process_gtex()
            if ok:
                logger.info("✅ GTEx processing complete")
            else:
                logger.error("❌ GTEx processing failed — see above")

    check_data_status()


if __name__ == "__main__":
    main()