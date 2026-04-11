"""
prepare_cmap_query.py  (v2.1 — import fix)
==========================================
Fixed: run_cmap_query_prep() now uses try/except relative import
instead of bare `from atrt_de_scorer import ...` which fails when
running as part of the backend.pipeline package.
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
POLL_INTERVAL = 30
MAX_WAIT_SECS = 1200

L1000_CONFIRMED_ATRT_GENES = {
    "EZH2", "EED", "SUZ12",
    "BRD4", "BRD2",
    "HDAC1", "HDAC2", "HDAC3", "HDAC6",
    "MYC", "MAX",
    "AURKA", "AURKB",
    "CDK4", "CDK6", "CCND1", "CDKN2A", "RB1",
    "MTOR", "PIK3CA", "AKT1", "PTEN",
    "SMO", "GLI1", "GLI2",
    "SOX2", "LIN28A",
    "BCL2", "BCL2L1", "MCL1",
    "SMARCB1",
    "OLIG2", "GFAP", "MBP", "MOG",
}

L1000_NON_LANDMARK_ATRT = {"MYCN", "PSMB5", "PSMB2", "PSMB1", "LIN28A", "SALL4"}


def validate_against_l1000(gene_list: List[str], api_key: Optional[str] = None) -> Dict:
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
                logger.info("L1000 landmark: %d/%d genes are landmark", len(in_lm), len(query_set))
                return {"landmark": sorted(in_lm), "non_landmark": sorted(not_lm)}
        except Exception as e:
            logger.warning("L1000 API check failed: %s — using local curated set.", e)

    in_lm  = query_set & L1000_CONFIRMED_ATRT_GENES
    not_lm = query_set - L1000_CONFIRMED_ATRT_GENES
    return {
        "landmark": sorted(in_lm),
        "non_landmark": sorted(not_lm),
        "note": "Local curated set — use api_key for full check",
    }


def format_cmap_query_files(
    up_genes:   List[str],
    down_genes: List[str],
    query_name: str = "ATRT_SMARCB1_loss_signature",
    de_stats:   Optional[Dict[str, float]] = None,
) -> Dict[str, Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    up_path   = OUTPUT_DIR / "atrt_up_genes.txt"
    down_path = OUTPUT_DIR / "atrt_down_genes.txt"
    with open(up_path, "w") as f:
        f.write("\n".join(up_genes))
    with open(down_path, "w") as f:
        f.write("\n".join(down_genes))

    payload = {
        "name":       query_name,
        "tool_id":    "sig_queryl1k_tool",
        "parameters": {
            "pert_type":   "trt_cp",
            "cell_id":     "all",
            "up_genes":    " ".join(up_genes),
            "dn_genes":    " ".join(down_genes),
            "ignore_dose": "1",
        },
    }

    payload_path = OUTPUT_DIR / "atrt_query_payload.json"
    with open(payload_path, "w") as f:
        json.dump(payload, f, indent=2)

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


def submit_to_clue_api(
    up_genes:   List[str],
    down_genes: List[str],
    api_key:    str,
    query_name: str = "ATRT_SMARCB1_loss_signature",
    poll:       bool = True,
) -> Optional[str]:
    if len(up_genes) > 500:
        up_genes = up_genes[:500]
    if len(down_genes) > 500:
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

    logger.info("Submitting ATRT CMap query to Clue.io API...")
    try:
        resp = requests.post(CLUE_API_URL, headers=headers, json=payload, timeout=30)
        resp.raise_for_status()
    except requests.HTTPError as e:
        if e.response.status_code == 401:
            logger.error("Clue.io API: 401 Unauthorized — check API key.")
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
        logger.info("Job %s: %s (%ds)", job_id, status, elapsed)
        if status == "completed":
            return _download_gct_result(job_id, status_data, headers)
        elif status in ("failed", "error"):
            logger.error("Job failed: %s", status_data.get("error"))
            return None

    logger.error("Job %s timed out after %ds", job_id, MAX_WAIT_SECS)
    return None


def _download_gct_result(job_id, status_data, headers):
    result_url = (
        status_data.get("result_url") or
        status_data.get("artifacts", {}).get("gct_url") or
        f"https://api.clue.io/api/jobs/{job_id}/result.gct"
    )
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
        return None


def run_cmap_query_prep(
    de_df,
    api_key:    Optional[str] = None,
    n_up:       int = 150,
    n_down:     int = 50,
    submit:     bool = False,
) -> Dict:
    """
    Top-level function: run DE → select genes → format → optionally submit.

    FIX v2.1: use try/except import to handle both package and direct module execution.
    """
    # ── FIX: safe relative/absolute import ──────────────────────────────────
    try:
        from .atrt_de_scorer import select_cmap_query_genes
    except ImportError:
        try:
            from atrt_de_scorer import select_cmap_query_genes
        except ImportError:
            logger.error(
                "Cannot import select_cmap_query_genes from atrt_de_scorer. "
                "Run this from the repo root: python scripts/01_prepare_cmap.py"
            )
            return {"success": False, "error": "atrt_de_scorer not importable"}

    up_genes, down_genes = select_cmap_query_genes(de_df, n_up=n_up, n_down=n_down)
    if not up_genes:
        return {"success": False, "error": "No upregulated genes found"}

    de_stats = None
    if "gene" in de_df.columns and "zscore_log2FC" in de_df.columns:
        de_stats = dict(zip(de_df["gene"], de_df["zscore_log2FC"]))

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