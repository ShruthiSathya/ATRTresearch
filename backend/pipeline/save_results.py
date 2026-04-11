"""
save_results.py  (v3.0)
========================
ATRT Drug Repurposing Pipeline — Save Results to JSON

FIXES FROM v2.1
---------------
1. Stats key consistency: accepts both "n_screened" and "n_drugs_screened"
   from discovery_pipeline.py (uses whichever is present).
2. "total_samples" and "total_atrt_samples" both accepted.
3. Improved error handling around pipeline initialisation.
4. Output txt log now uses a separate FileHandler so it doesn't conflict
   with the main log stream.

Run:
  python -m backend.pipeline.save_results --disease atrt --top_n 20
  python -m backend.pipeline.save_results --disease atrt --generic-only

Outputs:
  results/atrt_pipeline_results.json
  results/atrt_pipeline_output.txt
"""

import asyncio
import json
import logging
import argparse
import math
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

try:
    from .pipeline_config import BBB as BBB_CONFIG, COMPOSITE_WEIGHTS
except ImportError:
    from pipeline_config import BBB as BBB_CONFIG, COMPOSITE_WEIGHTS


def setup_logging(log_path: Path) -> None:
    """Configure logging to both stdout and a file."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # Clear existing handlers to avoid duplicates on re-run
    root_logger.handlers.clear()

    fmt = logging.Formatter("%(message)s")

    # Console handler
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(fmt)
    root_logger.addHandler(console)

    # File handler
    try:
        fh = logging.FileHandler(str(log_path), mode="w")
        fh.setFormatter(fmt)
        root_logger.addHandler(fh)
    except Exception as e:
        logging.warning("Could not create log file %s: %s", log_path, e)


def _safe(obj):
    """JSON serialisation helper — handles sets, NaN, Inf."""
    if isinstance(obj, set):
        return sorted(obj)
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    raise TypeError(f"Not JSON serialisable: {type(obj)}")


def _get_stat(stats: Dict, *keys, default=0):
    """Get first matching key from stats dict."""
    for k in keys:
        if k in stats:
            return stats[k]
    return default


def _normalise_candidate(c: dict) -> dict:
    """Normalise candidate dict to a consistent output schema."""
    return {
        "name":                    c.get("name") or c.get("drug_name") or "?",
        "targets":                 list(c.get("targets") or []),
        "mechanism":               c.get("mechanism") or c.get("drug_class", ""),
        "score":                   round(float(c.get("score", 0)), 4),
        "tissue_expression_score": round(float(c.get("tissue_expression_score", 0)), 4),
        "depmap_score":            round(float(c.get("depmap_score", 0)), 4),
        "ppi_score":               round(float(c.get("ppi_score", 0)), 4),
        "escape_bypass_score":     round(float(c.get("escape_bypass_score", 0)), 4),
        "bbb_penetrance":          c.get("bbb_penetrance", "UNKNOWN"),
        "bbb_score":               round(float(c.get("bbb_score", 0)), 4),
        "clinical_failure":        bool(c.get("clinical_failure", False)),
        # Generic drug fields
        "has_generic":             bool(c.get("has_generic", False)),
        "generic_source":          c.get("generic_source", ""),
        "cost_category":           c.get("cost_category", ""),
        # ATRT-specific fields
        "ezh2_boosted":            bool(c.get("ezh2_boosted", False)),
        "aurka_boosted":           bool(c.get("aurka_boosted", False)),
        "atrt_boosts_applied":     c.get("atrt_boosts_applied", []),
        # IC50 validation
        "ic50_validated":          bool(c.get("ic50_validated", False)),
        "ic50_um":                 c.get("ic50_um"),
        "ic50_cell_line":          c.get("ic50_cell_line"),
        "ic50_source":             c.get("ic50_source"),
        "ic50_has_smarcb1_data":   bool(c.get("ic50_has_smarcb1_data", False)),
        # CMap (neutral prior if not loaded)
        "cmap_score":              c.get("cmap_score", 0.50),
        "is_reverser":             bool(c.get("is_reverser", False)),
        # Notes
        "depmap_note":             c.get("depmap_note", ""),
        "escape_note":             c.get("escape_note", ""),
    }


def _extract_confidence_breakdown(hypotheses: list) -> Optional[Dict]:
    if not hypotheses:
        return None
    h  = hypotheses[0]
    cb = h.get("confidence_breakdown") or {}
    return {
        "drug_combo":               h.get("drug_or_combo", ""),
        "confidence":               round(float(h.get("confidence", 0) or 0), 4),
        "confidence_optimistic":    round(float(h.get("confidence_optimistic", 0) or 0), 4),
        "confidence_range":         h.get("confidence_range", ""),
        "priority":                 h.get("priority", ""),
        "statistical_significance": h.get("statistical_significance", ""),
        "depmap_essentiality":      round(float(cb.get("depmap_essentiality", 0) or 0), 4),
        "bbb_penetrance":           round(float(cb.get("bbb_penetrance", 0) or 0), 4),
        "mechanistic_diversity":    round(float(cb.get("mechanistic_diversity", 0) or 0), 4),
        "toxicity_flag":            cb.get("toxicity_flag", ""),
        "toxicity_note":            cb.get("toxicity_note", ""),
        "ezh2_boosted_in_combo":    cb.get("ezh2_inhibitors_in_combo", []),
    }


def _build_summary_table(candidates: list) -> list:
    rows = []
    for i, c in enumerate(candidates, 1):
        rows.append({
            "rank":           i,
            "drug":           c.get("name", "?"),
            "score":          round(float(c.get("score", 0)), 3),
            "bbb":            c.get("bbb_penetrance", "?"),
            "depmap":         round(float(c.get("depmap_score", 0)), 3),
            "tissue":         round(float(c.get("tissue_expression_score", 0)), 3),
            "cmap_score":     c.get("cmap_score", 0.50),
            "ic50_um":        c.get("ic50_um"),
            "ic50_cell_line": c.get("ic50_cell_line"),
            "ezh2_boosted":   bool(c.get("ezh2_boosted", False)),
            "aurka_boosted":  bool(c.get("aurka_boosted", False)),
            "has_generic":    bool(c.get("has_generic", False)),
        })
    return rows


async def main(
    disease:      str,
    top_n:        int,
    output_path:  Path,
    subgroup:     Optional[str] = None,
    location:     str  = "unknown_location",
    generic_only: bool = False,
) -> None:
    setup_logging(output_path.parent / "atrt_pipeline_output.txt")
    logger = logging.getLogger(__name__)

    logger.info("=" * 65)
    logger.info("ATRT Drug Repurposing Pipeline — save_results.py v3.0")
    logger.info("  disease      : %s", disease)
    logger.info("  top_n        : %d", top_n)
    logger.info("  subgroup     : %s", subgroup or "pan-ATRT")
    logger.info("  location     : %s", location)
    logger.info("  generic_only : %s", generic_only)
    logger.info("  timestamp    : %s", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("=" * 65)

    # Dynamic import to handle both package and direct run
    try:
        from backend.pipeline.discovery_pipeline import ProductionPipeline
        from backend.pipeline.cmap_query import CMAPQuery
    except ImportError:
        try:
            from pipeline.discovery_pipeline import ProductionPipeline
            from pipeline.cmap_query import CMAPQuery
        except ImportError:
            from discovery_pipeline import ProductionPipeline
            from cmap_query import CMAPQuery

    pipeline = ProductionPipeline(generic_only=generic_only)
    await pipeline.initialize(disease=disease)
    raw = await pipeline.run(
        disease_name=disease,
        top_k=top_n,
        subgroup=subgroup,
        location=location,
        generic_only=generic_only,
    )

    hypotheses = raw.get("hypotheses", [])
    stats      = raw.get("stats", {})
    candidates = raw.get("top_candidates", [])

    # Resolve stats keys — accept both naming conventions
    n_screened = _get_stat(
        stats, "n_screened", "n_drugs_screened", default=len(candidates)
    )
    total_samples = _get_stat(
        stats, "total_samples", "total_atrt_samples", default=0
    )

    # CMap status
    cmap = CMAPQuery()
    cmap_loaded = cmap.has_precomputed_scores()

    results = {
        "run_timestamp":    datetime.now().isoformat(),
        "disease":          disease,
        "subgroup":         subgroup or "pan-ATRT",
        "location":         location,
        "generic_only":     generic_only,
        "pipeline_version": "ATRT v3.0",

        "stats": {
            "smarcb1_loss_count":    _get_stat(stats, "smarcb1_loss_count"),
            "smarca4_loss_count":    _get_stat(stats, "smarca4_loss_count"),
            "total_atrt_samples":    total_samples,
            "total_samples":         total_samples,
            "rna_upregulated_genes": _get_stat(stats, "rna_upregulated_genes"),
            "p_value_label":         stats.get("p_value_label", "N/A"),
            "statistical_note":      stats.get("statistical_note", ""),
            "n_drugs_screened":      n_screened,
            "n_screened":            n_screened,
            "n_ezh2_boosted":        _get_stat(stats, "n_ezh2_boosted"),
            "n_aurka_boosted":       _get_stat(stats, "n_aurka_boosted"),
            "n_generic_in_top_k":    _get_stat(stats, "n_generic_in_top_k"),
            "escape_bypass_mode":    stats.get("escape_bypass_mode", "curated fallback"),
            "composite_weights":     COMPOSITE_WEIGHTS,
            "calling_method":        stats.get("calling_method", "none"),
            "cmap_loaded":           cmap_loaded,
            "depmap_source":         stats.get("depmap_source", "unknown"),
            "tissue_source":         stats.get("tissue_source", "unknown"),
            "cmap_note": (
                "Live CMap scores loaded"
                if cmap_loaded
                else (
                    "Neutral prior (0.50) applied — run scripts/01_prepare_cmap.py "
                    "then submit to https://clue.io to generate real scores"
                )
            ),
            "data_streams_active": [
                "OpenTargets API (EFO_0002915/0000543 — live; falls back to curated list if blocked)",
                "DepMap CRISPR (Broad Institute — ATRT/rhabdoid lines: BT16, BT37, G401, A204)",
                "GSE70678 bulk RNA-seq (Torchia 2015 Cancer Cell PMID 26609405 — 49 ATRT tumours)",
                "GTEx v8 normal brain reference (cerebellum + cortex)",
                "FDA/RxNorm generic drug annotation (live; seed fallback if blocked)",
                "STRING-DB PPI network (live + curated)",
                "ATRT published IC50 validation (BT16/BT37/G401/A204)",
                f"CMap L1000 transcriptomic reversal: {'live scores' if cmap_loaded else 'neutral prior (0.50)'}",
            ],
        },

        "top_candidates":      [_normalise_candidate(c) for c in candidates[:top_n]],
        "confidence_breakdown": _extract_confidence_breakdown(hypotheses),
        "hypotheses":           hypotheses,
        "atrt_summary":         _build_summary_table(candidates[:8]),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(str(output_path), "w") as f:
        json.dump(results, f, indent=2, default=_safe)

    logger.info("\n✅ Results saved to : %s", output_path)
    logger.info("✅ n_drugs_screened : %d", n_screened)
    logger.info("✅ CMap status      : %s", "live scores" if cmap_loaded else "neutral prior (0.50)")

    # ── Print candidate table ─────────────────────────────────────────────────
    logger.info("\n── SMARCB1 Statistics ───────────────────────────────────")
    logger.info(
        "  SMARCB1 loss     : %d/%d samples",
        _get_stat(stats, "smarcb1_loss_count"),
        total_samples,
    )
    logger.info("  Statistical note : %s", stats.get("p_value_label", "N/A"))

    logger.info("\n── Top 8 Candidates ─────────────────────────────────────────")
    logger.info(
        "  %-26s %6s  %-10s  %8s  %6s  %7s  IC50 (µM)",
        "Drug", "Score", "BBB", "Generic", "EZH2↑", "AURKA↑"
    )
    for c in candidates[:8]:
        ic50_str = (
            f"{c.get('ic50_um', 'N/A')} ({c.get('ic50_cell_line', '')})"
            if c.get("ic50_validated") else "—"
        )
        gen_str  = "✓ GENERIC" if c.get("has_generic") else "—"
        logger.info(
            "  %-26s %6.3f  %-10s  %8s  %6s  %7s  %s",
            c.get("name", "")[:26],
            c.get("score", 0),
            c.get("bbb_penetrance", "?"),
            gen_str,
            "YES" if c.get("ezh2_boosted") else "—",
            "YES" if c.get("aurka_boosted") else "—",
            ic50_str,
        )

    logger.info("\n  Total drugs screened: %d", n_screened)
    logger.info("  Depmap source: %s", stats.get("depmap_source", "unknown"))
    logger.info("  Tissue source: %s", stats.get("tissue_source", "unknown"))
    logger.info("\nNext step: python -m backend.pipeline.generate_figures")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="ATRT Pipeline — Save Results v3.0"
    )
    parser.add_argument("--disease",      default="atrt")
    parser.add_argument("--top_n",        default=20, type=int)
    parser.add_argument("--output",       default="results/atrt_pipeline_results.json")
    parser.add_argument("--subgroup",     default=None, choices=["TYR", "SHH", "MYC"])
    parser.add_argument("--location",     default="unknown_location",
                        choices=["infratentorial", "supratentorial", "unknown_location"])
    parser.add_argument("--generic-only", action="store_true",
                        help="Only include FDA-approved generic formulations")
    args = parser.parse_args()
    asyncio.run(main(
        disease      = args.disease,
        top_n        = args.top_n,
        output_path  = Path(args.output),
        subgroup     = args.subgroup,
        location     = args.location,
        generic_only = getattr(args, "generic_only", False),
    ))