"""
save_results.py
===============
ATRT Drug Repurposing Pipeline — Save Results to JSON

Run:
  python -m backend.pipeline.save_results --disease atrt --top_n 20

Outputs:
  results/atrt_pipeline_results.json   ← main results
  results/atrt_pipeline_output.txt     ← run log
"""

import asyncio
import json
import logging
import argparse
import math
import sys
from datetime import datetime
from pathlib import Path

try:
    from .pipeline_config import BBB as BBB_CONFIG, COMPOSITE_WEIGHTS
except ImportError:
    from pipeline_config import BBB as BBB_CONFIG, COMPOSITE_WEIGHTS


def setup_logging(log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handlers = [
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_path, mode="w"),
    ]
    logging.basicConfig(
        level=logging.INFO, format="%(message)s", handlers=handlers
    )


def _safe(obj):
    if isinstance(obj, set):
        return sorted(obj)
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    raise TypeError(f"Not serialisable: {type(obj)}")


def _normalise_candidate(c: dict) -> dict:
    """Normalise candidate dict to consistent output schema."""
    return {
        "name":                    c.get("name") or c.get("drug_name") or "?",
        "targets":                 list(c.get("targets") or []),
        "score":                   round(float(c.get("score", 0)), 4),
        "tissue_expression_score": round(float(c.get("tissue_expression_score", 0)), 4),
        "depmap_score":            round(float(c.get("depmap_score", 0)), 4),
        "ppi_score":               round(float(c.get("ppi_score", 0)), 4),
        "escape_bypass_score":     round(float(c.get("escape_bypass_score", 0)), 4),
        "bbb_penetrance":          c.get("bbb_penetrance", "UNKNOWN"),
        "bbb_score":               round(float(c.get("bbb_score", 0)), 4),
        "clinical_failure":        bool(c.get("clinical_failure", False)),
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
        # CMap
        "cmap_score":              c.get("cmap_score"),
        "is_reverser":             c.get("is_reverser"),
        # Scoring notes
        "depmap_note":             c.get("depmap_note", ""),
        "sc_context":              c.get("sc_context", ""),
        "escape_note":             c.get("escape_note", ""),
        "mechanism":               c.get("mechanism") or c.get("drug_class", ""),
    }


def _extract_confidence_breakdown(hypotheses: list) -> dict | None:
    if not hypotheses:
        return None
    h  = hypotheses[0]
    cb = h.get("confidence_breakdown") or {}
    return {
        "drug_combo":               h.get("drug_or_combo", ""),
        "confidence":               round(float(h.get("confidence", 0)), 4),
        "confidence_optimistic":    round(float(h.get("confidence_optimistic", 0)), 4),
        "confidence_range":         h.get("confidence_range", ""),
        "priority":                 h.get("priority", ""),
        "statistical_significance": h.get("statistical_significance", ""),
        "depmap_essentiality":      round(float(cb.get("depmap_essentiality", 0)), 4),
        "bbb_penetrance":           round(float(cb.get("bbb_penetrance", 0)), 4),
        "mechanistic_diversity":    round(float(cb.get("mechanistic_diversity", 0)), 4),
        "toxicity_flag":            cb.get("toxicity_flag", ""),
        "toxicity_note":            cb.get("toxicity_note", ""),
        "ezh2_boosted_in_combo":    cb.get("ezh2_inhibitors_in_combo", []),
    }


async def main(
    disease: str,
    top_n: int,
    output_path: Path,
    subgroup: str | None = None,
    location: str = "unknown_location",
) -> None:
    setup_logging(output_path.parent / "atrt_pipeline_output.txt")
    logger = logging.getLogger(__name__)

    logger.info("=" * 65)
    logger.info("ATRT Drug Repurposing Pipeline — save_results.py")
    logger.info(f"  disease   : {disease}")
    logger.info(f"  top_n     : {top_n}")
    logger.info(f"  subgroup  : {subgroup or 'pan-ATRT'}")
    logger.info(f"  location  : {location}")
    logger.info(f"  timestamp : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 65)

    from backend.pipeline.discovery_pipeline import ProductionPipeline

    pipeline = ProductionPipeline()
    await pipeline.initialize(disease=disease)
    raw = await pipeline.run(
        disease_name=disease,
        top_k=top_n,
        subgroup=subgroup,
        location=location,
    )

    hypotheses = raw.get("hypotheses", [])
    stats      = raw.get("stats", {})
    candidates = raw.get("top_candidates", [])

    n_screened = (
        stats.get("n_screened")
        or stats.get("n_drugs_screened")
        or len(candidates)
    )

    results = {
        "run_timestamp":  datetime.now().isoformat(),
        "disease":        disease,
        "subgroup":       subgroup or "pan-ATRT",
        "location":       location,
        "pipeline_version": "ATRT v1.0",

        "stats": {
            "smarcb1_loss_count":       stats.get("smarcb1_loss_count", 0),
            "smarca4_loss_count":       stats.get("smarca4_loss_count", 0),
            "total_atrt_samples":       stats.get("total_samples", 0),
            "rna_upregulated_genes":    stats.get("rna_upregulated_genes", 0),
            "p_value_label":            stats.get("p_value_label", "N/A"),
            "statistical_note":         stats.get("statistical_note", ""),
            "n_drugs_screened":         n_screened,
            "n_ezh2_boosted":           stats.get("n_ezh2_boosted", 0),
            "n_aurka_boosted":          stats.get("n_aurka_boosted", 0),
            "escape_bypass_mode":       stats.get("escape_bypass_mode", "curated fallback"),
            "composite_weights":        COMPOSITE_WEIGHTS,
            "calling_method":           stats.get("calling_method", "none"),
            "data_streams_active": [
                "DepMap CRISPR (Broad Institute — ATRT/rhabdoid lines: BT16, BT37, G401, A204)",
                "GSE70678 bulk RNA-seq (Torchia 2015 — 49 ATRT tumours)",
                "OpenTargets API (CNS/Oncology drugs)",
                "STRING-DB PPI network",
                "ATRT published IC50 validation (BT16/BT37/G401/A204)",
            ],
        },

        "top_candidates": [_normalise_candidate(c) for c in candidates[:top_n]],

        "top_combinations": [],   # populated if synergy predictor runs

        "confidence_breakdown": _extract_confidence_breakdown(hypotheses),
        "hypotheses":           hypotheses,

        "reports": {
            "novelty":          raw.get("novelty_report", ""),
            "polypharmacology": raw.get("poly_report", ""),
        },

        # ATRT-specific summary table (mirrors README format)
        "atrt_summary": _build_summary_table(candidates[:8]),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, default=_safe)

    logger.info(f"\n✅ Results saved to : {output_path}")
    logger.info(f"✅ n_drugs_screened : {n_screened}")
    logger.info(f"\nNext step: python -m backend.pipeline.generate_figures")

    # Summary table
    logger.info(f"\n── SMARCB1 Statistics ──────────────────────────────")
    logger.info(f"  SMARCB1 loss     : {stats.get('smarcb1_loss_count', 0)}/{stats.get('total_samples', 0)} samples")
    logger.info(f"  SMARCA4 loss     : {stats.get('smarca4_loss_count', 0)} samples")
    logger.info(f"  Statistical note : {stats.get('p_value_label', 'N/A')}")

    logger.info(f"\n── Top 8 Candidates ──────────────────────────────────")
    logger.info(f"  {'Drug':<26} {'Score':>6}  {'BBB':<10}  {'EZH2↑':>6}  {'AURKA↑':>7}  IC50 (µM)")
    for c in candidates[:8]:
        ic50_str = f"{c.get('ic50_um', 'N/A')} ({c.get('ic50_cell_line', '')})" if c.get("ic50_validated") else "—"
        logger.info(
            f"  {c.get('name',''):<26} {c.get('score',0):>6.3f}  "
            f"{c.get('bbb_penetrance','?'):<10}  "
            f"{'YES' if c.get('ezh2_boosted') else '—':>6}  "
            f"{'YES' if c.get('aurka_boosted') else '—':>7}  "
            f"{ic50_str}"
        )
    logger.info(f"\n  Total drugs screened (OpenTargets): {n_screened}")


def _build_summary_table(candidates: list) -> list:
    """Build summary table for README/report output."""
    rows = []
    for i, c in enumerate(candidates, 1):
        rows.append({
            "rank":       i,
            "drug":       c.get("name", "?"),
            "score":      round(float(c.get("score", 0)), 3),
            "bbb":        c.get("bbb_penetrance", "?"),
            "depmap":     round(float(c.get("depmap_score", 0)), 3),
            "tissue":     round(float(c.get("tissue_expression_score", 0)), 3),
            "cmap_norm_cs": c.get("cmap_score"),
            "ic50_um":    c.get("ic50_um"),
            "ic50_cell_line": c.get("ic50_cell_line"),
            "ezh2_boosted":  c.get("ezh2_boosted", False),
            "aurka_boosted": c.get("aurka_boosted", False),
        })
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ATRT Pipeline — Save Results")
    parser.add_argument("--disease",  default="atrt")
    parser.add_argument("--top_n",    default=20, type=int)
    parser.add_argument("--output",   default="results/atrt_pipeline_results.json")
    parser.add_argument("--subgroup", default=None, choices=["TYR", "SHH", "MYC"])
    parser.add_argument("--location", default="unknown_location",
                        choices=["infratentorial", "supratentorial", "unknown_location"])
    args = parser.parse_args()
    asyncio.run(main(
        args.disease, args.top_n, Path(args.output),
        subgroup=args.subgroup, location=args.location,
    ))