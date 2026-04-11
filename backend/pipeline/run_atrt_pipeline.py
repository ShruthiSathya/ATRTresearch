#!/usr/bin/env python3
"""
scripts/run_pipeline.py
=======================
Complete ATRT Drug Repurposing Pipeline runner.

USAGE
-----
# Full pipeline (all drugs)
python scripts/run_pipeline.py

# Generic drugs only (FDA-approved generics — cheaper, more accessible)
python scripts/run_pipeline.py --generic-only

# Specific subgroup
python scripts/run_pipeline.py --subgroup MYC
python scripts/run_pipeline.py --subgroup TYR
python scripts/run_pipeline.py --subgroup SHH

# With location (affects BBB penalty)
python scripts/run_pipeline.py --location infratentorial

# Combined
python scripts/run_pipeline.py --subgroup MYC --location supratentorial --generic-only

WHAT THIS RUNS
--------------
1. Data check  — verifies required files exist
2. OpenTargets — fetches CNS/oncology drug candidates (~500+ drugs)
3. DepMap      — CRISPR essentiality scores (ATRT/rhabdoid lines)
4. GSE70678    — tissue expression scores vs GTEx normal brain
5. Escape bypass — SMARCB1-loss resistance map
6. PPI network — STRING-DB proximity
7. Composite score — weighted combination
8. ATRT boosts — EZH2 ×1.40, AURKA ×1.15–1.30, SMO ×1.25 (SHH)
9. BBB filter  — location-aware penalty
10. Generic filter — (if --generic-only)
11. IC50 annotation — published ATRT cell-line data
12. Hypothesis — top triple combination
13. Figures    — 4 publication figures
14. Save results → results/atrt_pipeline_results.json

PREREQUISITES
-------------
pip install pandas numpy scipy matplotlib requests aiohttp fastapi uvicorn h5py

Data files (must exist before running):
  data/depmap/CRISPRGeneEffect.csv   ← download from depmap.org/portal
  data/depmap/Model.csv              ← same download

Optional (improves scoring — auto-downloaded if missing):
  data/raw_omics/GSE70678_series_matrix.txt.gz  ← or use pre-processed TSV
  data/raw_omics/GSE70678_gene_expression.tsv   ← (produced by --process gse70678)
  data/raw_omics/GTEx_brain_normal_reference.tsv ← (produced by --process gtex)
  data/cmap_query/atrt_cmap_scores.json         ← from clue.io query

HOW TO GET CMAP SCORES (optional but improves ranking)
------------------------------------------------------
1. Run: python scripts/01_prepare_cmap.py
        → creates data/cmap_query/atrt_up_genes.txt
                    data/cmap_query/atrt_down_genes.txt
2. Go to: https://clue.io → Tools → L1000 Query
3. Paste gene lists → Submit → wait ~10 min → Download query_result.gct
4. Run: python -m backend.pipeline.integrate_cmap_results
        → creates data/cmap_query/cmap_scores_pipeline.json
        (rename to atrt_cmap_scores.json for ATRT-specific results)
"""

import argparse
import asyncio
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "backend"))
sys.path.insert(0, str(REPO_ROOT / "backend" / "pipeline"))


# ─────────────────────────────────────────────────────────────────────────────
# Step 0 — Data check
# ─────────────────────────────────────────────────────────────────────────────

def check_data_files() -> bool:
    """Verify required data files exist. Print instructions if missing."""
    required = {
        "data/depmap/CRISPRGeneEffect.csv":
            "Download from https://depmap.org/portal/download/all/ (DepMap Public 24Q4)",
        "data/depmap/Model.csv":
            "Download from https://depmap.org/portal/download/all/ (same package)",
    }
    optional = {
        "data/raw_omics/GSE70678_gene_expression.tsv":
            "Run: python -m backend.pipeline.data_downloader --dataset gse70678 gpl6244 --process all",
        "data/raw_omics/GTEx_brain_normal_reference.tsv":
            "Run: python -m backend.pipeline.data_downloader --dataset gtex_brain --process gtex",
        "data/cmap_query/atrt_cmap_scores.json":
            "Run: python scripts/01_prepare_cmap.py → submit to clue.io → integrate_cmap_results.py",
    }

    all_ok = True
    print("\n" + "=" * 60)
    print("DATA FILE CHECK")
    print("=" * 60)

    for path_str, instructions in required.items():
        path = REPO_ROOT / path_str
        if path.exists():
            size_mb = path.stat().st_size / 1e6
            print(f"  ✅ {path_str} ({size_mb:.1f} MB)")
        else:
            print(f"  ❌ MISSING: {path_str}")
            print(f"     → {instructions}")
            all_ok = False

    print()
    for path_str, instructions in optional.items():
        path = REPO_ROOT / path_str
        if path.exists():
            size_mb = path.stat().st_size / 1e6
            print(f"  ✅ {path_str} ({size_mb:.1f} MB)")
        else:
            print(f"  ⚠️  OPTIONAL (missing): {path_str}")
            print(f"     → {instructions}")

    print("=" * 60)

    if not all_ok:
        print("\n❌ Required files missing. Download them before running the pipeline.")
        print("   Pipeline will still run using fallback values, but DepMap scores")
        print("   will be from published literature rather than your local CSV.")
    else:
        print("\n✅ All required files present.")

    return all_ok


# ─────────────────────────────────────────────────────────────────────────────
# Optional data download
# ─────────────────────────────────────────────────────────────────────────────

def download_and_process_optional_data():
    """Download and process GSE70678 + GTEx if not already present."""
    gene_expr = REPO_ROOT / "data/raw_omics/GSE70678_gene_expression.tsv"
    gtex_ref  = REPO_ROOT / "data/raw_omics/GTEx_brain_normal_reference.tsv"

    if gene_expr.exists() and gtex_ref.exists():
        logger.info("Optional RNA data already present — skipping download.")
        return

    logger.info("Downloading optional RNA data (GSE70678 + GTEx)...")
    import subprocess
    cmds = []

    if not (REPO_ROOT / "data/raw_omics/GSE70678_series_matrix.txt.gz").exists():
        cmds.append([
            sys.executable, "-m", "backend.pipeline.data_downloader",
            "--dataset", "gpl6244", "gse70678", "gtex_brain",
        ])

    if not gene_expr.exists() or not gtex_ref.exists():
        cmds.append([
            sys.executable, "-m", "backend.pipeline.data_downloader",
            "--process", "all",
        ])

    for cmd in cmds:
        result = subprocess.run(cmd, cwd=str(REPO_ROOT), capture_output=False)
        if result.returncode != 0:
            logger.warning("Data download/process step returned non-zero exit code.")


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline run
# ─────────────────────────────────────────────────────────────────────────────

async def run_pipeline(
    subgroup:     str  = None,
    location:     str  = "unknown_location",
    top_n:        int  = 20,
    generic_only: bool = False,
    output_path:  Path = None,
):
    from backend.pipeline.discovery_pipeline import ProductionPipeline

    output_path = output_path or REPO_ROOT / "results/atrt_pipeline_results.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("")
    logger.info("=" * 65)
    logger.info("ATRT Drug Repurposing Pipeline v1.0")
    logger.info("  Subgroup    : %s", subgroup or "pan-ATRT")
    logger.info("  Location    : %s", location)
    logger.info("  Generic only: %s", generic_only)
    logger.info("  Top N       : %d", top_n)
    logger.info("  Output      : %s", output_path)
    logger.info("=" * 65)

    pipeline = ProductionPipeline(generic_only=generic_only)
    await pipeline.initialize()

    result = await pipeline.run(
        disease_name  = "atrt",
        top_k         = top_n,
        subgroup      = subgroup,
        location      = location,
        generic_only  = generic_only,
    )

    # Print results table
    candidates = result.get("top_candidates", [])
    stats      = result.get("stats", {})

    print("\n" + "=" * 80)
    print(f"TOP {min(len(candidates), 12)} ATRT DRUG CANDIDATES")
    if generic_only:
        print("(GENERIC DRUGS ONLY — FDA-approved generics)")
    print("=" * 80)
    print(f"  {'Rank':<4}  {'Drug':<26}  {'Score':>6}  {'BBB':<10}  "
          f"{'EZH2↑':>5}  {'AURKA↑':>6}  {'Generic':>8}  IC50 (µM)")
    print("-" * 80)

    for i, c in enumerate(candidates[:12], 1):
        ic50_str  = f"{c.get('ic50_um')} ({c.get('ic50_cell_line', '')})" \
                    if c.get("ic50_validated") else "—"
        ezh2_str  = "YES" if c.get("ezh2_boosted")  else "—"
        aurka_str = "YES" if c.get("aurka_boosted") else "—"
        gen_str   = "✓" if c.get("has_generic") else "—"
        print(
            f"  {i:<4}  {c.get('name','?'):<26}  {c.get('score',0):>6.3f}  "
            f"{c.get('bbb_penetrance','?'):<10}  "
            f"{ezh2_str:>5}  {aurka_str:>6}  {gen_str:>8}  {ic50_str}"
        )

    # Hypothesis
    hyps = result.get("hypotheses", [])
    if hyps:
        h  = hyps[0]
        bd = h.get("confidence_breakdown", {})
        print(f"\n{'=' * 65}")
        print("TOP COMBINATION HYPOTHESIS")
        print(f"{'=' * 65}")
        print(f"  Combo : {h.get('drug_or_combo', '?')}")
        print(f"  Conf  : {bd.get('confidence_range', h.get('confidence', '?'))}")
        print(f"  EZH2  : {bd.get('ezh2_combo_note', 'N/A')[:60]}")

    # Stats
    print(f"\n{'=' * 65}")
    print("PIPELINE STATISTICS")
    print(f"{'=' * 65}")
    print(f"  Drugs screened   : {stats.get('n_screened', '?')}")
    print(f"  SMARCB1 loss     : {stats.get('smarcb1_loss_count', 0)}/{stats.get('total_samples', 0)}")
    print(f"  RNA upregulated  : {stats.get('rna_upregulated_genes', 0)} genes")
    print(f"  EZH2 boosted     : {stats.get('n_ezh2_boosted', 0)}")
    print(f"  Escape bypass    : {stats.get('escape_bypass_mode', 'N/A')}")
    print(f"  Generic only     : {stats.get('generic_only', False)}")

    # Save JSON
    def _safe(obj):
        import math
        if isinstance(obj, set):  return sorted(obj)
        if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)): return None
        raise TypeError(type(obj))

    with open(output_path, "w") as f:
        json.dump({
            "run_timestamp":   datetime.now().isoformat(),
            "disease":         "atrt",
            "subgroup":        subgroup or "pan-ATRT",
            "location":        location,
            "generic_only":    generic_only,
            "stats":           stats,
            "top_candidates":  candidates,
            "hypotheses":      hyps,
        }, f, indent=2, default=_safe)

    logger.info("\n✅ Results saved → %s", output_path)

    # Optional: generate figures
    try:
        from backend.pipeline.generate_figures import (
            fig1_smarcb1_biology, fig2_drug_rankings,
            fig3_score_scatter, fig4_confidence,
        )
        data_for_figs = {
            "stats": {
                "smarcb1_loss_count":   stats.get("smarcb1_loss_count", 0),
                "smarca4_loss_count":   stats.get("smarca4_loss_count", 0),
                "total_atrt_samples":   stats.get("total_samples", 0),
                "n_drugs_screened":     stats.get("n_screened", 0),
                "n_ezh2_boosted":       stats.get("n_ezh2_boosted", 0),
                "n_aurka_boosted":      stats.get("n_aurka_boosted", 0),
            },
            "top_candidates": candidates,
            "hypotheses":     hyps,
            "confidence_breakdown": hyps[0].get("confidence_breakdown") if hyps else {},
        }
        fig1_smarcb1_biology(data_for_figs)
        fig2_drug_rankings(data_for_figs, top_n=min(len(candidates), 12))
        fig3_score_scatter(data_for_figs)
        fig4_confidence(data_for_figs)
        logger.info("✅ Figures saved → figures/atrt/")
    except Exception as e:
        logger.warning("Figures skipped: %s", e)

    return result


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="ATRT Drug Repurposing Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full pan-ATRT run
  python scripts/run_pipeline.py

  # Generic drugs only (valproic acid, sirolimus, etc.)
  python scripts/run_pipeline.py --generic-only

  # Subgroup-specific
  python scripts/run_pipeline.py --subgroup MYC --location supratentorial

  # Check data files only
  python scripts/run_pipeline.py --check-data
        """,
    )
    parser.add_argument("--subgroup",     default=None, choices=["TYR", "SHH", "MYC"],
                        help="ATRT molecular subgroup (default: pan-ATRT)")
    parser.add_argument("--location",     default="unknown_location",
                        choices=["infratentorial", "supratentorial", "unknown_location"],
                        help="Tumour location — affects BBB penalty")
    parser.add_argument("--top-n",        default=20, type=int,
                        help="Number of top candidates to report (default: 20)")
    parser.add_argument("--generic-only", action="store_true",
                        help="Keep only drugs with FDA-approved generic formulations")
    parser.add_argument("--output",       default=None,
                        help="Output JSON path (default: results/atrt_pipeline_results.json)")
    parser.add_argument("--check-data",   action="store_true",
                        help="Check data files and exit")
    parser.add_argument("--download",     action="store_true",
                        help="Download optional RNA data (GSE70678 + GTEx) before running")
    args = parser.parse_args()

    # Data check
    check_data_files()

    if args.check_data:
        sys.exit(0)

    if args.download:
        download_and_process_optional_data()

    output = Path(args.output) if args.output else None
    asyncio.run(run_pipeline(
        subgroup     = args.subgroup,
        location     = args.location,
        top_n        = args.top_n,
        generic_only = args.generic_only,
        output_path  = output,
    ))


if __name__ == "__main__":
    main()