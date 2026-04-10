"""
data_downloader.py
==================
Automated data download helpers for the ATRT Drug Repurposing Pipeline.

Run this script FIRST before executing save_results.py.

Usage:
    python -m backend.pipeline.data_downloader --all
    python -m backend.pipeline.data_downloader --dataset gse70678
    python -m backend.pipeline.data_downloader --dataset depmap
    python -m backend.pipeline.data_downloader --check   # just verify what's present

DATA SOURCES AND PROVENANCE
----------------------------
All datasets are open-access (no login required except CBTN).

1. GSE70678 — Torchia 2015 ATRT RNA-seq cohort
   Reference: Torchia J et al. Cancer Cell 2015; 30(6):891-908. PMID 26609405.
   URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70678
   Content: 49 ATRT tumors + normal brain samples; subgroup annotations (TYR/SHH/MYC)
   Platform: Affymetrix HuGene 1.0 ST (GPL6244)
   Size: ~50 MB compressed
   License: GEO open access

2. DepMap 24Q4 — Broad Institute CRISPR Chronos scores
   Reference: Behan FM et al. Nature 2019; 568:511. PMID 30971826.
   URL: https://depmap.org/portal/download/all/
   Files: CRISPRGeneEffect.csv (~500 MB), Model.csv (~5 MB)
   Release: 24Q4 (October 2024) — latest as of April 2026
   License: CC BY 4.0
   ATRT/rhabdoid lines included: BT16, BT37, G401, A204 (verify in Model.csv)

3. clue.io CMap L1000 — Transcriptomic reversal scores
   Reference: Subramanian A et al. Cell 2017; 171(6):1437. PMID 29195078.
   URL: https://clue.io (free academic account)
   Process: Run prepare_cmap_query.py → submit gene lists → download GCT → integrate
   Steps documented in prepare_cmap_query.py

4. OpenTargets API — Drug-target associations (already integrated, no download)
   Reference: Ochoa D et al. Nucleic Acids Res 2021; 49(D1):D1302. PMID 33290552.
   URL: https://api.platform.opentargets.org/api/v4/graphql (live API)

5. CBTN ATRT cohort (OPTIONAL — controlled access)
   URL: https://portal.kidsfirstdrc.org
   Study: PBTA ATRT samples (~30-40 samples)
   Access: Requires DCC agreement (same as PNOC/PBTA)
   Files: cna.txt, rna_zscores.txt → place in data/validation/cbtn_genomics/atrt/
   If absent: pipeline uses GSE70678 fallback (fully functional)
"""

import argparse
import gzip
import io
import json
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
# Dataset definitions
# ─────────────────────────────────────────────────────────────────────────────

DATASETS = {
    "gse70678": {
        "name":        "GSE70678 — Torchia 2015 ATRT RNA-seq",
        "url":         "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70678/matrix/GSE70678_series_matrix.txt.gz",
        "output_path": "data/raw_omics/GSE70678_ATRT_expression.txt",
        "compressed":  True,
        "size_mb":     50,
        "pmid":        "26609405",
        "license":     "GEO open access",
        "instructions": (
            "After download, the file will be decompressed automatically.\n"
            "The series matrix contains expression data for 49 ATRT + normal samples.\n"
            "Pipeline reads this file directly — no additional processing needed."
        ),
    },
    "gse106982": {
        "name":        "GSE106982 — Johann 2016 ATRT methylation subgroups (optional)",
        "url":         "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106982/matrix/GSE106982_series_matrix.txt.gz",
        "output_path": "data/raw_omics/GSE106982_ATRT_methylation.txt",
        "compressed":  True,
        "size_mb":     30,
        "pmid":        "26923874",
        "license":     "GEO open access",
        "optional":    True,
        "instructions": (
            "Optional: provides methylation-based subgroup labels (TYR/SHH/MYC) for 150 ATRT.\n"
            "Pipeline uses prevalence priors from Johann 2016 if this file is absent."
        ),
    },
}

# DepMap must be downloaded manually due to portal authentication
DEPMAP_MANUAL_INSTRUCTIONS = """
═══════════════════════════════════════════════════════════════
DepMap 24Q4 — Manual Download Required
═══════════════════════════════════════════════════════════════
DepMap requires a free account to download data files.

Steps:
1. Go to: https://depmap.org/portal/download/all/
2. Select release: "24Q4" (October 2024 — latest as of April 2026)
3. Download these two files:
   - CRISPRGeneEffect.csv  (~500 MB)
   - Model.csv             (~5 MB)
4. Place both files in: data/depmap/

mkdir -p data/depmap
# Then move downloaded files:
mv ~/Downloads/CRISPRGeneEffect.csv data/depmap/
mv ~/Downloads/Model.csv data/depmap/

Reference: Behan FM et al. Nature 2019; 568:511. PMID 30971826.
License: CC BY 4.0

ATRT/rhabdoid cell lines to look for in Model.csv:
  BT16  (ACH-000725) — CNS ATRT, ATRT-MYC, SMARCB1-null
  BT37  (ACH-000881) — CNS ATRT, ATRT-TYR, SMARCB1-null
  G401  (ACH-000039) — Renal rhabdoid, SMARCB1-null (Knutson 2013 validation line)
  A204  (ACH-000658) — Rhabdoid, SMARCB1-null
  BT12  (ACH-001082) — CNS ATRT, SMARCB1-null (verify ID in your Model.csv)
═══════════════════════════════════════════════════════════════
"""

CMAP_MANUAL_INSTRUCTIONS = """
═══════════════════════════════════════════════════════════════
clue.io CMap L1000 — Semi-automated (requires free account)
═══════════════════════════════════════════════════════════════
Steps:
1. First prepare the ATRT gene signature:
   python -m backend.pipeline.prepare_cmap_query

2. This creates:
   data/cmap_query/atrt_up_genes.txt   (top upregulated genes)
   data/cmap_query/atrt_down_genes.txt (top downregulated genes)

3. Go to: https://clue.io → create free academic account
4. Click "Query" → "L1000 Query"
5. Paste contents of atrt_up_genes.txt into UP box
6. Paste contents of atrt_down_genes.txt into DOWN box
7. Name the query: "ATRT_SMARCB1_signature"
8. Submit → results ready in ~10 minutes

9. Download results GCT file → save as:
   data/cmap_query/query_result.gct

10. Integrate results:
    python -m backend.pipeline.integrate_cmap_results

Reference: Subramanian A et al. Cell 2017; 171(6):1437. PMID 29195078.
License: clue.io free academic

Note: If GSE70678 is not yet downloaded, prepare_cmap_query.py will
use curated ATRT signature genes from literature as a starting point.
═══════════════════════════════════════════════════════════════
"""

CBTN_MANUAL_INSTRUCTIONS = """
═══════════════════════════════════════════════════════════════
CBTN ATRT Cohort — Controlled Access (Optional)
═══════════════════════════════════════════════════════════════
This is OPTIONAL — pipeline works fully with GSE70678 fallback.

If you want patient-level ATRT genomics:
1. Go to: https://portal.kidsfirstdrc.org
2. Register and apply for data access (DCC agreement)
3. Search for PBTA ATRT samples (~30-40 samples)
4. Download: cna.txt, rna_zscores.txt, mutations.txt
5. Place files in: data/validation/cbtn_genomics/atrt/

mkdir -p data/validation/cbtn_genomics/atrt/

Reference: Rokita JL et al. Nat Cancer 2021. PMID 33842919.
═══════════════════════════════════════════════════════════════
"""


# ─────────────────────────────────────────────────────────────────────────────
# Download functions
# ─────────────────────────────────────────────────────────────────────────────

def _progress_hook(count: int, block_size: int, total_size: int) -> None:
    """Simple download progress indicator."""
    if total_size > 0:
        pct = min(int(count * block_size * 100 / total_size), 100)
        print(f"\r  Progress: {pct}%", end="", flush=True)


def download_file(
    url: str,
    output_path: str,
    compressed: bool = False,
    max_retries: int = 3,
) -> bool:
    """Download a file from URL to output_path."""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, max_retries + 1):
        try:
            logger.info("Downloading %s (attempt %d/%d)...", url, attempt, max_retries)
            tmp_path = out.with_suffix(out.suffix + ".tmp")

            urllib.request.urlretrieve(url, tmp_path, reporthook=_progress_hook)
            print()  # newline after progress

            if compressed and url.endswith(".gz"):
                logger.info("Decompressing...")
                with gzip.open(tmp_path, "rb") as f_in:
                    with open(out, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                tmp_path.unlink()
            else:
                tmp_path.rename(out)

            size_mb = out.stat().st_size / (1024 * 1024)
            logger.info("✅ Downloaded: %s (%.1f MB)", out, size_mb)
            return True

        except Exception as e:
            logger.warning("Attempt %d failed: %s", attempt, e)
            if attempt < max_retries:
                time.sleep(5 * attempt)

    logger.error("❌ All %d download attempts failed for %s", max_retries, url)
    return False


def check_data_status() -> Dict[str, bool]:
    """Check which data files are present."""
    checks = {
        "GSE70678 (tissue expression)": Path("data/raw_omics/GSE70678_ATRT_expression.txt"),
        "GSE106982 (methylation subgroups, optional)": Path("data/raw_omics/GSE106982_ATRT_methylation.txt"),
        "DepMap CRISPRGeneEffect.csv": Path("data/depmap/CRISPRGeneEffect.csv"),
        "DepMap Model.csv": Path("data/depmap/Model.csv"),
        "CMap scores JSON": Path("data/cmap_query/atrt_cmap_scores.json"),
        "CBTN ATRT genomics (optional)": Path("data/validation/cbtn_genomics/atrt/"),
    }

    print("\n" + "═" * 60)
    print("ATRT Pipeline Data Status")
    print("═" * 60)

    status = {}
    for label, path in checks.items():
        present = path.exists()
        size_str = ""
        if present and path.is_file():
            size_mb = path.stat().st_size / (1024 * 1024)
            size_str = f"  ({size_mb:.1f} MB)"
        icon = "✅" if present else "❌"
        opt  = " [optional]" if "optional" in label else ""
        print(f"  {icon}  {label}{opt}{size_str}")
        status[label] = present

    # Minimum viable check
    viable = (
        checks["DepMap CRISPRGeneEffect.csv"].exists()
        and checks["DepMap Model.csv"].exists()
    )
    gse_present = checks["GSE70678 (tissue expression)"].exists()

    print()
    if viable and gse_present:
        print("  ✅ Pipeline is FULLY operational")
    elif viable:
        print("  ⚠️  Pipeline is FUNCTIONAL (GSE70678 absent — will use curated scores)")
        print("     Run: python -m backend.pipeline.data_downloader --dataset gse70678")
    else:
        print("  ❌ Pipeline NOT operational — DepMap files required")
        print("     See DepMap download instructions above")

    print("═" * 60 + "\n")
    return status


def download_gse70678() -> bool:
    """Download GSE70678 from NCBI GEO."""
    ds = DATASETS["gse70678"]
    logger.info("Downloading %s", ds["name"])
    logger.info("Reference: PMID %s", ds["pmid"])
    success = download_file(ds["url"], ds["output_path"], compressed=ds["compressed"])
    if success:
        logger.info(ds["instructions"])
    return success


def download_gse106982() -> bool:
    """Download GSE106982 from NCBI GEO (optional)."""
    ds = DATASETS["gse106982"]
    logger.info("Downloading %s", ds["name"])
    return download_file(ds["url"], ds["output_path"], compressed=ds["compressed"])


def print_all_instructions() -> None:
    """Print all manual download instructions."""
    print(DEPMAP_MANUAL_INSTRUCTIONS)
    print(CMAP_MANUAL_INSTRUCTIONS)
    print(CBTN_MANUAL_INSTRUCTIONS)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ATRT Pipeline Data Downloader",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--dataset",
        choices=["gse70678", "gse106982", "all"],
        help="Which dataset to download",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check which data files are present",
    )
    parser.add_argument(
        "--instructions",
        action="store_true",
        help="Print manual download instructions for DepMap/CMap/CBTN",
    )
    args = parser.parse_args()

    if args.check or (not args.dataset and not args.instructions):
        check_data_status()

    if args.instructions:
        print_all_instructions()
        return

    if args.dataset == "gse70678" or args.dataset == "all":
        download_gse70678()

    if args.dataset == "gse106982" or args.dataset == "all":
        download_gse106982()

    if args.dataset == "all":
        print(DEPMAP_MANUAL_INSTRUCTIONS)
        print(CMAP_MANUAL_INSTRUCTIONS)

    check_data_status()


if __name__ == "__main__":
    main()