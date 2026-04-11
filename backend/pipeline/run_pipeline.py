#!/usr/bin/env python3
"""
run_pipeline.py
===============
ATRT Drug Repurposing Pipeline — Standalone Runner v2.0

USAGE
-----
# Full pan-ATRT run (all drugs)
python run_pipeline.py

# Generic drugs only (FDA-approved generics — more accessible)
python run_pipeline.py --generic-only

# Subgroup-specific
python run_pipeline.py --subgroup MYC
python run_pipeline.py --subgroup TYR
python run_pipeline.py --subgroup SHH

# With location (affects BBB penalty)
python run_pipeline.py --location infratentorial

# Combined
python run_pipeline.py --subgroup MYC --location supratentorial --generic-only

# Check data files only
python run_pipeline.py --check-data

WHAT THIS NEEDS
---------------
Required (for live DepMap scoring):
  data/depmap/CRISPRGeneEffect.csv
  data/depmap/Model.csv

Optional (auto-detected, significantly improves tissue scoring):
  data/raw_omics/GSE70678_gene_expression.tsv  ← YOUR FILE
  data/raw_omics/GTEx_brain_normal_reference.tsv ← YOUR FILE
  data/raw_omics/GPL570_probe_map.tsv           ← YOUR FILE
  data/raw_omics/GPL570.annot.gz
  data/raw_omics/GSE70678_series_matrix.txt.gz ← YOUR FILE

If any optional file is missing, the pipeline uses curated literature-based
scores and clearly logs which data source is active.

OUTPUT
------
  results/atrt_pipeline_results.json
  results/atrt_pipeline_output.txt (log)
  Printed ranked table to stdout
"""

import argparse
import asyncio
import json
import logging
import math
import sys
from datetime import datetime
from pathlib import Path

# ─── Path setup: add backend/ to sys.path so imports work ────────────────────
_THIS_DIR = Path(__file__).resolve().parent
_BACKEND  = _THIS_DIR / "backend"
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

# ─── Logging ─────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s — %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("atrt_runner")


def _safe(obj):
    if isinstance(obj, set):
        return sorted(obj)
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    raise TypeError(f"Not serialisable: {type(obj)}")


# ─────────────────────────────────────────────────────────────────────────────
# Data file check
# ─────────────────────────────────────────────────────────────────────────────

def check_data_files(repo_root: Path) -> bool:
    required = {
        "data/depmap/CRISPRGeneEffect.csv": (
            "https://depmap.org/portal/download/all/ → DepMap Public 24Q4"
        ),
        "data/depmap/Model.csv": "Same download as CRISPRGeneEffect.csv",
    }
    optional = {
        "data/raw_omics/GSE70678_gene_expression.tsv":  "GEO GSE70678 processed matrix",
        "data/raw_omics/GTEx_brain_normal_reference.tsv":"GTEx v8 brain tissues",
        "data/raw_omics/GPL570_probe_map.tsv":           "GPL570 Affymetrix probe map",
        "data/raw_omics/GPL570.annot.gz":                "GPL570 annotation (raw)",
        "data/raw_omics/GSE70678_series_matrix.txt.gz":  "GEO GSE70678 series matrix",
        "data/cmap_query/atrt_cmap_scores.json":          "clue.io CMap scores (optional)",
    }

    print("\n" + "=" * 65)
    print("DATA FILE STATUS")
    print("=" * 65)

    all_required_ok = True
    for rel, instructions in required.items():
        p = repo_root / rel
        ok = p.exists()
        size_str = f"({p.stat().st_size / 1e6:.1f} MB)" if ok else ""
        print(f"  {'✅' if ok else '❌'} {'REQUIRED':8} {rel} {size_str}")
        if not ok:
            print(f"              → {instructions}")
            all_required_ok = False

    print()
    for rel, desc in optional.items():
        p = repo_root / rel
        ok = p.exists()
        size_str = f"({p.stat().st_size / 1e6:.2f} MB)" if ok else ""
        print(f"  {'✅' if ok else '⚪'} {'OPTIONAL':8} {rel} {size_str}")
        if not ok:
            print(f"              → {desc}")

    print("=" * 65)
    if not all_required_ok:
        print(
            "\n⚠️  Required DepMap files missing.\n"
            "   Pipeline will run using VERIFIED FALLBACK Chronos values\n"
            "   (sourced from Knutson 2013, Geoerger 2017, Sredni 2017, Lin 2019).\n"
            "   Download DepMap CSV files to get live cell-line-specific scores.\n"
        )
    else:
        print("\n✅ Required files present — live DepMap scoring enabled.\n")
    return all_required_ok


# ─────────────────────────────────────────────────────────────────────────────
# Main run
# ─────────────────────────────────────────────────────────────────────────────

async def run_pipeline(
    subgroup:     Optional[str] = None,
    location:     str  = "unknown_location",
    top_n:        int  = 20,
    generic_only: bool = False,
    output_path:  Path = None,
    repo_root:    Path = None,
) -> Dict:
    from pipeline.discovery_pipeline import ProductionPipeline

    output_path = output_path or (repo_root / "results/atrt_pipeline_results.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 65)
    print("ATRT Drug Repurposing Pipeline v2.0")
    print(f"  Subgroup     : {subgroup or 'pan-ATRT'}")
    print(f"  Location     : {location}")
    print(f"  Generic only : {generic_only}")
    print(f"  Top N        : {top_n}")
    print(f"  Timestamp    : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 65)

    pipeline = ProductionPipeline(generic_only=generic_only)
    await pipeline.initialize()

    result = await pipeline.run(
        disease_name  = "atrt",
        top_k         = top_n,
        subgroup      = subgroup,
        location      = location,
        generic_only  = generic_only,
    )

    candidates = result.get("top_candidates", [])
    stats      = result.get("stats", {})
    hyps       = result.get("hypotheses", [])

    # ── Print ranked table ───────────────────────────────────────────────────
    print(f"\n{'='*80}")
    if generic_only:
        print("TOP GENERIC DRUG CANDIDATES FOR ATRT")
        print("(Restricted to FDA-approved generic formulations)")
    else:
        print(f"TOP {min(len(candidates), top_n)} DRUG CANDIDATES FOR ATRT")
    print(f"{'='*80}")
    print(
        f"  {'Rk':<3}  {'Drug':<24}  {'Score':>6}  {'BBB':<10}  "
        f"{'Generic':>8}  {'EZH2↑':>5}  {'AURKA↑':>6}  IC50 (µM)"
    )
    print("-" * 80)

    for i, c in enumerate(candidates[:top_n], 1):
        ic50_str = (
            f"{c.get('ic50_um')} ({c.get('ic50_cell_line', '')})"
            if c.get("ic50_validated") else "—"
        )
        gen_str  = "✓ GENERIC" if c.get("has_generic") else "—"
        ezh2_str = "✓" if c.get("ezh2_boosted") else "—"
        aurk_str = "✓" if c.get("aurka_boosted") else "—"
        print(
            f"  {i:<3}  {c.get('name','?'):<24}  {c.get('score',0):>6.3f}  "
            f"{c.get('bbb_penetrance','?'):<10}  "
            f"{gen_str:>8}  {ezh2_str:>5}  {aurk_str:>6}  {ic50_str}"
        )

    # ── Score breakdown ──────────────────────────────────────────────────────
    print(f"\n{'='*65}")
    print("SCORE COMPONENT BREAKDOWN (top 8)")
    print(f"{'='*65}")
    print(f"  {'Drug':<24}  {'Total':>6}  {'Tissue':>6}  {'DepMap':>6}  {'Escape':>6}  {'PPI':>4}")
    print("-" * 65)
    for c in candidates[:8]:
        print(
            f"  {c.get('name','?'):<24}  {c.get('score',0):>6.3f}  "
            f"{c.get('tissue_expression_score',0):>6.3f}  "
            f"{c.get('depmap_score',0):>6.3f}  "
            f"{c.get('escape_bypass_score',0):>6.3f}  "
            f"{c.get('ppi_score',0):>4.2f}"
        )

    # ── Top hypothesis ───────────────────────────────────────────────────────
    if hyps:
        h  = hyps[0]
        bd = h.get("confidence_breakdown", {})
        print(f"\n{'='*65}")
        print("TOP COMBINATION HYPOTHESIS")
        print(f"{'='*65}")
        print(f"  Combo          : {h.get('drug_or_combo', '?')}")
        print(f"  Confidence     : {bd.get('confidence_range', h.get('confidence', '?'))}")
        print(f"  Priority       : {h.get('priority', '?')}")
        print(f"  SMARCB1 note   : {h.get('statistical_note', 'N/A')[:60]}")
        if bd.get("ezh2_inhibitors_in_combo"):
            print(f"  EZH2 in combo  : {bd['ezh2_inhibitors_in_combo']} ← POSITIVE in ATRT")
        if bd.get("toxicity_flag"):
            print(f"  Toxicity       : {bd.get('toxicity_flag')} — {bd.get('toxicity_note','')[:50]}")

    # ── Pipeline statistics ───────────────────────────────────────────────────
    print(f"\n{'='*65}")
    print("PIPELINE STATISTICS")
    print(f"{'='*65}")
    print(f"  Drugs screened       : {stats.get('n_screened', '?')}")
    print(f"  SMARCB1 loss         : {stats.get('smarcb1_loss_count',0)}/{stats.get('total_samples',0)}")
    print(f"  RNA upregulated genes: {stats.get('rna_upregulated_genes',0)}")
    print(f"  EZH2 inhibitors boost: {stats.get('n_ezh2_boosted',0)}")
    print(f"  Generic in top-{top_n}    : {stats.get('n_generic_in_top_k',0)}")
    print(f"  DepMap source        : {stats.get('depmap_source','?')}")
    print(f"  Tissue source        : {stats.get('tissue_source','?')}")
    print(f"  Escape bypass mode   : {stats.get('escape_bypass_mode','?')}")
    print(f"  Genomics method      : {stats.get('calling_method','?')}")

    # ── Generic drug summary ─────────────────────────────────────────────────
    generic_top = [c for c in candidates[:top_n] if c.get("has_generic")]
    if generic_top:
        print(f"\n{'='*65}")
        print(f"GENERIC DRUGS IN TOP {top_n} (n={len(generic_top)})")
        print(f"{'='*65}")
        for c in generic_top:
            print(
                f"  {c.get('name','?'):<24}  score={c.get('score',0):.3f}  "
                f"BBB={c.get('bbb_penetrance','?')}  "
                f"mechanism: {(c.get('mechanism','')[:40])}"
            )

    # ── Save JSON ────────────────────────────────────────────────────────────
    output = {
        "run_timestamp":    datetime.now().isoformat(),
        "disease":          "atrt",
        "subgroup":         subgroup or "pan-ATRT",
        "location":         location,
        "generic_only":     generic_only,
        "pipeline_version": "v2.0",
        "stats":            stats,
        "top_candidates":   candidates,
        "hypotheses":       hyps,
        "generic_candidates": [
            {
                "rank":       i + 1,
                "name":       c.get("name"),
                "score":      c.get("score"),
                "bbb":        c.get("bbb_penetrance"),
                "mechanism":  c.get("mechanism"),
                "ic50_um":    c.get("ic50_um"),
                "ic50_cell":  c.get("ic50_cell_line"),
                "depmap":     c.get("depmap_score"),
                "tissue":     c.get("tissue_expression_score"),
                "escape":     c.get("escape_bypass_score"),
            }
            for i, c in enumerate(generic_top)
        ],
    }

    with open(str(output_path), "w") as f:
        json.dump(output, f, indent=2, default=_safe)

    logger.info("✅ Results saved → %s", output_path)
    print(f"\n✅ Results saved to: {output_path}")

    return result


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

from typing import Optional, Dict  # noqa: E402


def main():
    parser = argparse.ArgumentParser(
        description="ATRT Drug Repurposing Pipeline v2.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_pipeline.py                          # Full pan-ATRT run
  python run_pipeline.py --generic-only           # Generic drugs only
  python run_pipeline.py --subgroup MYC           # MYC subgroup
  python run_pipeline.py --check-data             # Verify data files
  python run_pipeline.py --generic-only --top-n 15 --subgroup SHH
        """,
    )
    parser.add_argument("--subgroup",     default=None, choices=["TYR", "SHH", "MYC"])
    parser.add_argument("--location",     default="unknown_location",
                        choices=["infratentorial", "supratentorial", "unknown_location"])
    parser.add_argument("--top-n",        default=20, type=int)
    parser.add_argument("--generic-only", action="store_true",
                        help="Restrict to FDA-approved generic formulations only")
    parser.add_argument("--output",       default=None,
                        help="Output JSON path (default: results/atrt_pipeline_results.json)")
    parser.add_argument("--check-data",   action="store_true",
                        help="Check data files and exit")
    args = parser.parse_args()

    repo_root = _THIS_DIR
    check_data_files(repo_root)
    if args.check_data:
        sys.exit(0)

    output = Path(args.output) if args.output else None
    asyncio.run(run_pipeline(
        subgroup     = args.subgroup,
        location     = args.location,
        top_n        = args.top_n,
        generic_only = args.generic_only,
        output_path  = output,
        repo_root    = repo_root,
    ))


if __name__ == "__main__":
    main()