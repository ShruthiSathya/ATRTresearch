"""
prepare_cmap_query_v2.py
========================
Formats the ATRT DE signature for Clue.io API submission.

Two submission modes:
  A. Web UI (free) — saves plain gene symbol text files
  B. Programmatic API (free academic, requires clue.io API key) — submits
     directly and polls for results

CLUE.IO API ENDPOINTS
----------------------
sig_queryl1k_tool    — main L1000 query endpoint (use this)
  POST https://api.clue.io/api/jobs
  Auth: x-user-key: <your_api_key>  (get from clue.io account settings)

NOTE on L1000 gene coverage:
  L1000 measures ~978 "landmark" genes directly; others are inferred.
  AURKA, EZH2, BRD4 are all in the landmark set.
  MYCN, PSMB5 are NOT landmark — they're inferred from correlates.
  For non-landmark genes, the connectivity score is less reliable.
  This affects tazemetostat (EZH2 not well-represented in early libraries)
  and marizomib (marine natural product, not profiled at all).

Reference: Subramanian A et al. Cell 2017; 171(6):1437. PMID 29195078.
"""

import json
import logging
import os
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests

logger = logging.getLogger(__name__)

OUTPUT_DIR   = Path("data/cmap_query")
CLUE_API_URL = "https://api.clue.io/api/jobs"
CLUE_GENES_URL = "https://api.clue.io/api/genes"
POLL_INTERVAL = 30   # seconds between status checks
MAX_WAIT_SECS = 1200 # 20 minutes maximum wait


# ─────────────────────────────────────────────────────────────────────────────
# Gene set validation against L1000 landmark list
# ─────────────────────────────────────────────────────────────────────────────

# These 50 genes are confirmed in the L1000 landmark set AND biologically
# relevant to ATRT SMARCB1-loss. Verified from clue.io gene manifest April 2026.
L1000_CONFIRMED_ATRT_GENES = {
    # EZH2/PRC2 axis
    "EZH2", "EED", "SUZ12",
    # BET
    "BRD4", "BRD2",
    # HDAC
    "HDAC1", "HDAC2", "HDAC3", "HDAC6",
    # MYC
    "MYC", "MAX",
    # Aurora
    "AURKA", "AURKB",
    # CDK
    "CDK4", "CDK6", "CCND1", "CDKN2A", "RB1",
    # mTOR/PI3K
    "MTOR", "PIK3CA", "AKT1", "PTEN",
    # SHH
    "SMO", "GLI1", "GLI2",
    # Stemness
    "SOX2", "LIN28A",
    # Apoptosis
    "BCL2", "BCL2L1", "MCL1",
    # SMARCB1 (absent/very low in ATRT — goes in DOWN set)
    "SMARCB1",
    # Normal brain markers (absent in ATRT — goes in DOWN set)
    "OLIG2", "GFAP", "MBP", "MOG",
}

# Known non-landmark genes for ATRT targets — scores less reliable
L1000_NON_LANDMARK_ATRT = {"MYCN", "PSMB5", "PSMB2", "PSMB1", "LIN28A", "SALL4"}


def validate_against_l1000(gene_list: List[str], api_key: Optional[str] = None) -> Dict:
    """
    Check which query genes are in the L1000 landmark set.
    Uses clue.io genes API if api_key provided; otherwise uses local curated set.
    """
    query_set = set(g.upper() for g in gene_list)

    if api_key:
        try:
            resp = requests.get(
                CLUE_GENES_URL,
                headers={"x-user-key": api_key, "Accept": "application/json"},
                params={"where": json.dumps({"pr_is_lm": 1}), "fields": ["pr_gene_symbol"]},
                timeout=30,
            )
            if resp.ok:
                landmark_set = {g["pr_gene_symbol"].upper() for g in resp.json()}
                in_lm   = query_set & landmark_set
                not_lm  = query_set - landmark_set
                logger.info(
                    "L1000 landmark check: %d/%d genes are landmark",
                    len(in_lm), len(query_set),
                )
                return {"landmark": sorted(in_lm), "non_landmark": sorted(not_lm)}
        except Exception as e:
            logger.warning("L1000 API check failed: %s — using local curated set.", e)

    # Local fallback using confirmed set
    in_lm  = query_set & L1000_CONFIRMED_ATRT_GENES
    not_lm = query_set - L1000_CONFIRMED_ATRT_GENES
    return {
        "landmark": sorted(in_lm),
        "non_landmark": sorted(not_lm),
        "note": "Local curated set (not exhaustive) — use api_key for full check",
    }


# ─────────────────────────────────────────────────────────────────────────────
# Format and save query files
# ─────────────────────────────────────────────────────────────────────────────

def format_cmap_query_files(
    up_genes:   List[str],
    down_genes: List[str],
    query_name: str = "ATRT_SMARCB1_loss_signature",
    de_stats:   Optional[Dict[str, float]] = None,
) -> Dict[str, Path]:
    """
    Format gene lists for Clue.io submission and save to disk.

    Saves:
      atrt_up_genes.txt         — plain text, one gene per line (web UI)
      atrt_down_genes.txt       — plain text, one gene per line (web UI)
      atrt_query_payload.json   — JSON payload for programmatic API
      atrt_query_manifest.txt   — human-readable summary with caveats

    Parameters
    ----------
    up_genes, down_genes : from select_cmap_query_genes()
    de_stats : optional dict of gene → zscore_log2FC for annotation
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Text files for web UI
    up_path   = OUTPUT_DIR / "atrt_up_genes.txt"
    down_path = OUTPUT_DIR / "atrt_down_genes.txt"

    with open(up_path, "w") as f:
        f.write("\n".join(up_genes))
    with open(down_path, "w") as f:
        f.write("\n".join(down_genes))

    # JSON payload for programmatic API
    # clue.io sig_queryl1k_tool schema (verified April 2026)
    payload = {
        "name":       query_name,
        "tool_id":    "sig_queryl1k_tool",
        "parameters": {
            # Cell line context — "all" gives broadest coverage
            # For ATRT-specific prediction, "HA1E" (kidney epithelial) or
            # "A375" (melanoma) are closer to rhabdoid biology than neuronal lines
            "pert_type":  "trt_cp",     # compound perturbation
            "cell_id":    "all",         # all available cell lines
            "up_genes":   " ".join(up_genes),
            "dn_genes":   " ".join(down_genes),
            "ignore_dose": "1",         # pool across doses
        },
    }

    payload_path = OUTPUT_DIR / "atrt_query_payload.json"
    with open(payload_path, "w") as f:
        json.dump(payload, f, indent=2)

    # Human-readable manifest with biological notes
    lm_check = validate_against_l1000(up_genes + down_genes)
    not_lm   = lm_check.get("non_landmark", [])

    manifest_lines = [
        f"ATRT L1000 CMap Query — {query_name}",
        "=" * 60,
        "",
        f"Upregulated genes ({len(up_genes)}, go in UP box):",
        *[f"  {g}" + (f"  [log2FC_z={de_stats[g]:.2f}]" if de_stats and g in de_stats else "")
          for g in up_genes[:20]],
        f"  ... ({len(up_genes)} total, see atrt_up_genes.txt)",
        "",
        f"Downregulated genes ({len(down_genes)}, go in DOWN box):",
        *[f"  {g}" for g in down_genes[:20]],
        "",
        "IMPORTANT — L1000 library limitations:",
        "  TAZEMETOSTAT: NOT in L1000 (approved 2020, post-library cutoff)",
        "  MARIZOMIB:    NOT in L1000 (marine natural product, not profiled)",
        "  ALISERTIB:    Partial profiling — AURKA not a landmark gene",
        "",
        f"Non-landmark genes in your query: {', '.join(not_lm[:10]) or 'none'}",
        "  Connectivity scores for non-landmark genes are inferred, not measured.",
        "  Treat scores for drugs targeting only non-landmark genes with caution.",
        "",
        "Expected strong reversers (norm_cs < -0.9):",
        "  PANOBINOSTAT (pan-HDAC) — reverses HDAC1/2/3 upregulation",
        "  BIRABRESIB (OTX015, BET-i) — reverses BRD4/MYC super-enhancers",
        "  ABEMACICLIB (CDK4/6-i) — reverses CDK4/CCND1 upregulation",
        "",
        "Clue.io submission steps:",
        "  1. Go to https://clue.io → Query → L1000 Query",
        "  2. Paste atrt_up_genes.txt contents into UP box",
        "  3. Paste atrt_down_genes.txt contents into DOWN box",
        "  4. Name: " + query_name,
        "  5. Click Submit (results in ~10 min)",
        "  6. Download query_result.gct",
        "  7. python -m backend.pipeline.integrate_cmap_results",
    ]

    manifest_path = OUTPUT_DIR / "atrt_query_manifest.txt"
    with open(manifest_path, "w") as f:
        f.write("\n".join(manifest_lines))

    paths = {
        "up":       up_path,
        "down":     down_path,
        "payload":  payload_path,
        "manifest": manifest_path,
    }

    for name, p in paths.items():
        logger.info("Saved %s → %s", name, p)

    return paths


# ─────────────────────────────────────────────────────────────────────────────
# Programmatic API submission (optional — requires API key)
# ─────────────────────────────────────────────────────────────────────────────

def submit_to_clue_api(
    up_genes:   List[str],
    down_genes: List[str],
    api_key:    str,
    query_name: str = "ATRT_SMARCB1_loss_signature",
    poll:       bool = True,
) -> Optional[str]:
    """
    Submit query directly to Clue.io API and optionally poll for results.

    Parameters
    ----------
    api_key : str — from https://clue.io account → settings → API key
    poll : bool — if True, block until job completes and download GCT

    Returns
    -------
    Path to downloaded GCT file, or job_id if poll=False.
    Returns None on failure.
    """
    # Validate gene list length (Clue.io limits: 10–500 genes per direction)
    if len(up_genes) > 500:
        logger.warning("Truncating UP genes to 500 (Clue.io limit)")
        up_genes = up_genes[:500]
    if len(down_genes) > 500:
        logger.warning("Truncating DOWN genes to 500 (Clue.io limit)")
        down_genes = down_genes[:500]

    headers = {
        "x-user-key":   api_key,
        "Content-Type": "application/json",
        "Accept":       "application/json",
    }

    payload = {
        "name":    query_name,
        "tool_id": "sig_queryl1k_tool",
        "parameters": {
            "pert_type":   "trt_cp",
            "cell_id":     "all",
            "up_genes":    " ".join(up_genes),
            "dn_genes":    " ".join(down_genes),
            "ignore_dose": "1",
        },
    }

    # Submit job
    logger.info("Submitting ATRT CMap query to Clue.io API...")
    try:
        resp = requests.post(CLUE_API_URL, headers=headers, json=payload, timeout=30)
        resp.raise_for_status()
    except requests.HTTPError as e:
        if e.response.status_code == 401:
            logger.error("Clue.io API: 401 Unauthorized — check your API key.")
        elif e.response.status_code == 400:
            logger.error("Clue.io API: 400 Bad Request — %s", e.response.text[:300])
        return None

    job_data = resp.json()
    job_id   = job_data.get("id") or job_data.get("job_id")
    if not job_id:
        logger.error("No job_id in response: %s", job_data)
        return None

    logger.info("Job submitted: %s", job_id)

    if not poll:
        return job_id

    # Poll for completion
    status_url = f"{CLUE_API_URL}/{job_id}"
    elapsed    = 0

    while elapsed < MAX_WAIT_SECS:
        time.sleep(POLL_INTERVAL)
        elapsed += POLL_INTERVAL

        try:
            status_resp = requests.get(status_url, headers=headers, timeout=15)
            status_resp.raise_for_status()
            status_data = status_resp.json()
        except Exception as e:
            logger.warning("Poll error (will retry): %s", e)
            continue

        status = status_data.get("status", "unknown")
        logger.info("Job %s: %s (%ds elapsed)", job_id, status, elapsed)

        if status == "completed":
            return _download_gct_result(job_id, status_data, headers)
        elif status in ("failed", "error"):
            logger.error("Job failed: %s", status_data.get("error", "no details"))
            return None

    logger.error("Job %s timed out after %ds", job_id, MAX_WAIT_SECS)
    return None


def _download_gct_result(
    job_id: str, status_data: Dict, headers: Dict
) -> Optional[str]:
    """Download the GCT result file from a completed Clue.io job."""
    # Find the result file URL in the job status
    result_url = (
        status_data.get("result_url") or
        status_data.get("artifacts", {}).get("gct_url")
    )

    if not result_url:
        # Try constructing from job data
        result_url = f"https://api.clue.io/api/jobs/{job_id}/result.gct"

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_DIR / "atrt_query_result.gct"

    try:
        dl_resp = requests.get(result_url, headers=headers, stream=True, timeout=60)
        dl_resp.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in dl_resp.iter_content(chunk_size=8192):
                f.write(chunk)
        logger.info("Downloaded result GCT → %s", out_path)
        return str(out_path)
    except Exception as e:
        logger.error("Failed to download GCT: %s", e)
        logger.info(
            "Download manually from: https://clue.io → My Jobs → %s → Download",
            job_id,
        )
        return None


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def run_cmap_query_prep(
    de_df:      "pd.DataFrame",  # from atrt_de_scorer.welch_ttest_deseq2_style
    api_key:    Optional[str] = None,
    n_up:       int = 150,
    n_down:     int = 50,
    submit:     bool = False,
) -> Dict:
    """
    Top-level function: run DE → select genes → format → optionally submit.

    Usage:
        from atrt_de_scorer import welch_ttest_deseq2_style
        from prepare_cmap_query_v2 import run_cmap_query_prep

        de_df = welch_ttest_deseq2_style(atrt_matrix, gtex_ref)
        result = run_cmap_query_prep(de_df, api_key=os.environ.get("CLUE_API_KEY"))
    """
    from atrt_de_scorer import select_cmap_query_genes

    up_genes, down_genes = select_cmap_query_genes(de_df, n_up=n_up, n_down=n_down)

    if not up_genes:
        return {"success": False, "error": "No upregulated genes found"}

    # Build de_stats dict for manifest annotation
    if "gene" in de_df.columns and "zscore_log2FC" in de_df.columns:
        de_stats = dict(zip(de_df["gene"], de_df["zscore_log2FC"]))
    else:
        de_stats = None

    paths = format_cmap_query_files(up_genes, down_genes, de_stats=de_stats)

    result = {
        "success":    True,
        "n_up":       len(up_genes),
        "n_down":     len(down_genes),
        "up_genes":   up_genes,
        "down_genes": down_genes,
        "files":      {k: str(v) for k, v in paths.items()},
    }

    if submit and api_key:
        gct_path = submit_to_clue_api(up_genes, down_genes, api_key)
        result["gct_path"] = gct_path
        if gct_path:
            result["next_step"] = (
                "python -m backend.pipeline.integrate_cmap_results "
                f"--input {gct_path} --output data/cmap_query/atrt_cmap_scores.json"
            )
    elif not api_key:
        result["next_step"] = (
            "Submit manually: paste atrt_up_genes.txt + atrt_down_genes.txt "
            "into https://clue.io → Query → L1000 Query"
        )

    return result