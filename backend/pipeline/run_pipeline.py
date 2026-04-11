#!/usr/bin/env python3
"""
run_pipeline.py
===============
ATRT Drug Repurposing Pipeline — Standalone Runner v2.1

USAGE
-----
python run_pipeline.py                                  # Full pan-ATRT
python run_pipeline.py --generic-only                   # Generic drugs only
python run_pipeline.py --subgroup MYC                   # MYC subgroup
python run_pipeline.py --subgroup SHH --location infratentorial
python run_pipeline.py --check-data                     # Verify data files

WHAT THIS NEEDS
---------------
Required (for live DepMap scoring):
  data/depmap/CRISPRGeneEffect.csv
  data/depmap/Model.csv

Optional (auto-detected):
  data/raw_omics/GSE70678_gene_expression.tsv
  data/raw_omics/GTEx_brain_normal_reference.tsv
  data/raw_omics/GPL570_probe_map.tsv
  data/cmap_query/atrt_cmap_scores.json  ← run scripts/01_prepare_cmap.py first
"""

import argparse
import asyncio
import json
import logging
import math
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional  # ← must be at top

# ── Path setup ────────────────────────────────────────────────────────────────
_THIS_DIR = Path(__file__).resolve().parent
_BACKEND  = _THIS_DIR / "backend"
for _p in [str(_BACKEND), str(_THIS_DIR)]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

# ── Logging ───────────────────────────────────────────────────────────────────
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


# ── Data file check ───────────────────────────────────────────────────────────

def check_data_files(repo_root: Path) -> bool:
    required = {
        "data/depmap/CRISPRGeneEffect.csv": (
            "https://depmap.org/portal/download/all/ → DepMap Public 24Q4"
        ),
        "data/depmap/Model.csv": "Same download as CRISPRGeneEffect.csv",
    }
    optional = {
        "data/raw_omics/GSE70678_gene_expression.tsv":   "GEO GSE70678 processed (already present ✅)",
        "data/raw_omics/GTEx_brain_normal_reference.tsv": "GTEx v8 brain tissues (already present ✅)",
        "data/raw_omics/GPL570_probe_map.tsv":            "GPL570 probe map (already present ✅)",
        "data/raw_omics/GPL570.annot.gz":                 "GPL570 annotation (already present ✅)",
        "data/raw_omics/GSE70678_series_matrix.txt.gz":   "GEO GSE70678 series matrix (already present ✅)",
        "data/cmap_query/atrt_cmap_scores.json": (
            "Run: python scripts/01_prepare_cmap.py → submit to clue.io → "
            "python -m backend.pipeline.integrate_cmap_results"
        ),
    }

    print("\n" + "=" * 70)
    print("DATA FILE STATUS")
    print("=" * 70)

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

    print("=" * 70)
    if not all_required_ok:
        print(
            "\n⚠️  Required DepMap files missing.\n"
            "   Pipeline will run using VERIFIED FALLBACK Chronos values.\n"
            "   Download DepMap CSV files for live cell-line-specific scores.\n"
        )
    else:
        print("\n✅ Required files present — live DepMap scoring enabled.\n")

    # CMap status
    cmap_path = repo_root / "data/cmap_query/atrt_cmap_scores.json"
    if cmap_path.exists():
        print("✅ CMap scores present — transcriptomic reversal scoring enabled.\n")
    else:
        print(
            "⚪ CMap scores not yet generated.\n"
            "   All drugs will receive neutral CMap prior (0.50).\n"
            "   To generate: python scripts/01_prepare_cmap.py\n"
            "   Then submit gene lists to clue.io and run integrate_cmap_results.py\n"
        )
    return all_required_ok


# ── Main run ──────────────────────────────────────────────────────────────────

async def run_pipeline(
    subgroup:     Optional[str] = None,
    location:     str  = "unknown_location",
    top_n:        int  = 20,
    generic_only: bool = False,
    output_path:  Optional[Path] = None,
    repo_root:    Optional[Path] = None,
) -> Dict:
    from pipeline.discovery_pipeline import ProductionPipeline

    repo_root   = repo_root or _THIS_DIR
    output_path = output_path or (repo_root / "results/atrt_pipeline_results.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 70)
    print("ATRT Drug Repurposing Pipeline v2.1")
    print(f"  Subgroup     : {subgroup or 'pan-ATRT'}")
    print(f"  Location     : {location}")
    print(f"  Generic only : {generic_only}")
    print(f"  Top N        : {top_n}")
    print(f"  Timestamp    : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

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
    print(f"\n{'='*85}")
    header = "TOP GENERIC DRUG CANDIDATES FOR ATRT" if generic_only else f"TOP {min(len(candidates), top_n)} DRUG CANDIDATES FOR ATRT"
    print(header)
    print(f"{'='*85}")
    print(
        f"  {'Rk':<3}  {'Drug':<24}  {'Score':>6}  {'BBB':<10}  "
        f"{'Generic':>8}  {'EZH2↑':>5}  {'AURKA↑':>6}  IC50 (µM)"
    )
    print("-" * 85)

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
    print(f"\n{'='*70}")
    print("SCORE COMPONENT BREAKDOWN (top 8)")
    print(f"{'='*70}")
    print(f"  {'Drug':<24}  {'Total':>6}  {'Tissue':>6}  {'DepMap':>6}  {'Escape':>6}  {'PPI':>4}  {'CMap':>5}")
    print("-" * 70)
    for c in candidates[:8]:
        cmap_str = f"{c.get('cmap_score', 0):.2f}" if c.get("cmap_score") is not None else "N/A"
        print(
            f"  {c.get('name','?'):<24}  {c.get('score',0):>6.3f}  "
            f"{c.get('tissue_expression_score',0):>6.3f}  "
            f"{c.get('depmap_score',0):>6.3f}  "
            f"{c.get('escape_bypass_score',0):>6.3f}  "
            f"{c.get('ppi_score',0):>4.2f}  "
            f"{cmap_str:>5}"
        )

    # ── Top hypothesis ───────────────────────────────────────────────────────
    if hyps:
        h  = hyps[0]
        bd = h.get("confidence_breakdown", {})
        print(f"\n{'='*70}")
        print("TOP COMBINATION HYPOTHESIS")
        print(f"{'='*70}")
        print(f"  Combo          : {h.get('drug_or_combo', '?')}")
        print(f"  Confidence     : {bd.get('confidence_range', h.get('confidence', '?'))}")
        print(f"  Priority       : {h.get('priority', '?')}")
        smarcb1_note = h.get('statistical_note', 'N/A')
        print(f"  SMARCB1 note   : {smarcb1_note[:70]}")
        if bd.get("ezh2_inhibitors_in_combo"):
            print(f"  EZH2 in combo  : {bd['ezh2_inhibitors_in_combo']} ← POSITIVE in ATRT (synthetic lethality)")
        if bd.get("toxicity_flag"):
            tox_note = bd.get('toxicity_note', '')[:55]
            print(f"  Toxicity       : {bd.get('toxicity_flag')} — {tox_note}")

    # ── Pipeline statistics ───────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("PIPELINE STATISTICS")
    print(f"{'='*70}")
    print(f"  Drugs screened         : {stats.get('n_screened', '?')}")
    print(f"  SMARCB1 loss           : {stats.get('smarcb1_loss_count',0)}/{stats.get('total_samples',0)}")
    print(f"  RNA upregulated genes  : {stats.get('rna_upregulated_genes',0)}")
    print(f"  EZH2 inhibitors boosted: {stats.get('n_ezh2_boosted',0)}")
    print(f"  Generic in top-{top_n:<2}       : {stats.get('n_generic_in_top_k',0)}")
    print(f"  DepMap source          : {stats.get('depmap_source','?')}")
    print(f"  Tissue source          : {stats.get('tissue_source','?')}")
    print(f"  Escape bypass mode     : {stats.get('escape_bypass_mode','?')}")
    print(f"  Genomics calling method: {stats.get('calling_method','?')}")
    print(f"  CMap data              : {'live scores' if stats.get('cmap_loaded') else 'neutral prior (0.50) — run 01_prepare_cmap.py'}")

    # ── Generic drug summary ─────────────────────────────────────────────────
    generic_top = [c for c in candidates[:top_n] if c.get("has_generic")]
    if generic_top:
        print(f"\n{'='*70}")
        print(f"GENERIC DRUGS IN TOP {top_n} (n={len(generic_top)}) — More accessible for pediatric patients")
        print(f"{'='*70}")
        for c in generic_top:
            mech = (c.get('mechanism', '')[:45] or '—')
            print(
                f"  {c.get('name','?'):<26}  score={c.get('score',0):.3f}  "
                f"BBB={c.get('bbb_penetrance','?')}  {mech}"
            )

    # ── Save JSON ────────────────────────────────────────────────────────────
    output = {
        "run_timestamp":    datetime.now().isoformat(),
        "disease":          "atrt",
        "subgroup":         subgroup or "pan-ATRT",
        "location":         location,
        "generic_only":     generic_only,
        "pipeline_version": "v2.1",
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
    print("   Next step: python -m backend.pipeline.generate_figures\n")

    return result


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="ATRT Drug Repurposing Pipeline v2.1",
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