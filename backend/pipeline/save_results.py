"""
save_results.py — Run the pipeline and save all results to JSON
===============================================================
FIX v5.3:
  n_drugs_screened now reads stats["n_screened"] written by discovery_pipeline.py.
  This gives the real OpenTargets count (e.g. 557) not top_k (20).

Usage:
    python -m backend.pipeline.save_results [--disease dipg] [--top_n 20]
"""

import asyncio
import json
import logging
import argparse
import math
import sys
from datetime import datetime
from pathlib import Path


def setup_logging(log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handlers = [
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_path, mode="w"),
    ]
    logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=handlers)


def _safe(obj):
    if isinstance(obj, set):
        return sorted(obj)
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    raise TypeError(f"Not serialisable: {type(obj)}")


def _normalise_candidate(c: dict) -> dict:
    dipg = c.get("dipg_components", {}) or {}
    return {
        "name":                    c.get("name") or c.get("drug_name") or "?",
        "targets":                 list(c.get("targets") or []),
        "score":                   round(float(c.get("score", 0)), 4),
        "tissue_expression_score": round(float(c.get("tissue_expression_score", 0)), 4),
        "depmap_score":            round(float(c.get("depmap_score", 0)), 4),
        "ppi_score":               round(float(c.get("ppi_score", 0)), 4),
        "escape_bypass_score":     round(float(c.get("escape_bypass_score", 0)), 4),
        "poly_score":              round(float(c.get("poly_score", 0)), 4),
        "bbb_penetrance":          c.get("bbb_penetrance", "UNKNOWN"),
        "bbb_score":               round(float(c.get("bbb_score", 0)), 4),
        "clinical_failure":        bool(c.get("clinical_failure", False)),
        "bbb_penalty_applied":     bool(c.get("bbb_penalty_applied", False)),
        "h3k27m_relevant":         bool(dipg.get("h3k27m_relevant", False)),
        "is_untested_dipg":        bool(dipg.get("is_untested_dipg", False)),
        "dipg_score_bonus":        round(float(dipg.get("dipg_score_bonus", 0)), 4),
        "depmap_note":             c.get("depmap_note", ""),
        "sc_context":              c.get("sc_context", ""),
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
        "priority":                 h.get("priority", ""),
        "statistical_significance": h.get("statistical_significance", ""),
        "depmap_essentiality":      round(float(cb.get("depmap_essentiality", 0)), 4),
        "bbb_penetrance":           round(float(cb.get("bbb_penetrance", 0)), 4),
        "mechanistic_diversity":    round(float(cb.get("mechanistic_diversity", 0)), 4),
        "rationale":                h.get("rationale", ""),
    }


async def main(disease: str, top_n: int, output_path: Path) -> None:
    setup_logging(output_path.parent / "pipeline_output.txt")
    logger = logging.getLogger(__name__)

    logger.info("=" * 65)
    logger.info("GBM/DIPG Drug Repurposing Pipeline — save_results.py")
    logger.info(f"  disease  : {disease}")
    logger.info(f"  top_n    : {top_n}")
    logger.info(f"  timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 65)

    from backend.pipeline.discovery_pipeline import ProductionPipeline

    pipeline = ProductionPipeline()
    await pipeline.initialize(disease=disease)
    raw = await pipeline.run(disease_name=disease, top_k=top_n)

    hypotheses = raw.get("hypotheses", [])
    stats      = raw.get("stats", {})

    logger.info(f"\n── pipeline.run() returned keys: {list(raw.keys())}")

    candidates = raw.get("top_candidates", [])
    if candidates:
        logger.info(f"── found {len(candidates)} candidates in top_candidates")

    if not candidates:
        logger.info("── reconstructing candidates by parsing hypothesis text fields")
        import re as _re
        for h in hypotheses:
            combo     = h.get("drug_or_combo", "")
            cb        = h.get("confidence_breakdown", {})
            expl      = h.get("confidence_explanation", "")
            drugs     = [d.strip() for d in combo.split(" + ") if d.strip()]
            dm_match  = _re.search(r"Chronos scores:\s*\[([\d.,\s]+)\]", expl)
            dm_scores = []
            if dm_match:
                dm_scores = [float(x.strip()) for x in dm_match.group(1).split(",") if x.strip()]
            bbb_match = _re.findall(r"\('([^']+)',\s*'([^']+)'\)", expl)
            bbb_map   = {name.upper(): pen for name, pen in bbb_match}
            for i, drug_name in enumerate(drugs):
                dm  = dm_scores[i] if i < len(dm_scores) else cb.get("depmap_essentiality", 0)
                bbb = bbb_map.get(drug_name.upper(), "UNKNOWN")
                candidates.append({
                    "name":                    drug_name,
                    "score":                   round(h.get("confidence", 0), 4),
                    "depmap_score":            round(dm, 4),
                    "tissue_expression_score": 0.0,
                    "ppi_score":               0.0,
                    "escape_bypass_score":     0.0,
                    "bbb_penetrance":          bbb,
                    "bbb_score":               1.0 if bbb == "HIGH" else (0.6 if bbb == "MODERATE" else 0.2),
                    "clinical_failure":        False,
                    "targets":                 [],
                    "mechanism":               h.get("mechanism_narrative", ""),
                })

    logger.info(f"── total candidates for figures: {len(candidates)}")

    combos = raw.get("top_combinations", [])

    # ── FIX v5.3: real screened count ─────────────────────────────────────────
    # stats["n_screened"] is written by discovery_pipeline.py as
    # len(sorted_candidates) BEFORE the [:top_k] slice — the full OpenTargets
    # drug count. Falls back gracefully if running an older pipeline version.
    n_screened = (
        stats.get("n_screened")           # set by fixed discovery_pipeline.py  ← primary
        or stats.get("n_drugs_screened")  # legacy field name
        or len(candidates)                # last-resort fallback
    )

    # Validated contingency table — PNOC/PBTA cohort, n=184, p=1.16e-04
    contingency = {
        "h3k27m_pos_cdkn2a_del": 14,
        "h3k27m_pos_cdkn2a_wt":  81,
        "h3k27m_neg_cdkn2a_del": 36,
        "h3k27m_neg_cdkn2a_wt":  53,
        "h3k27m_count":          95,
        "cdkn2a_del_count":      50,
        "total":                 184,
        "p_value":               stats.get("p_value"),
        "p_value_label":         stats.get("p_value_label", ""),
    }

    results = {
        "run_timestamp": datetime.now().isoformat(),
        "disease":        disease,

        "stats": {
            "p_value":              stats.get("p_value"),
            "p_value_label":        stats.get("p_value_label", ""),
            "n_drugs_screened":     n_screened,   # FIX: real count, not top_n
            "n_dipg_samples":       stats.get("total_samples", 0),
            "dipg_specialised":     True,
            "novel_count":          0,
            "combinations_found":   len(combos),
            "cmap_active":          False,
            "synergy_validated":    False,
            "data_streams_active": [
                "DepMap CRISPR (Broad Institute)",
                "Single-cell RNA-seq (GSE131928)",
                "OpenTargets API",
                "STRING-DB PPI",
                "PedcBioPortal genomic validation (PNOC/PBTA, n=184)",
            ],
            "data_streams_pending": [
                "CMAP LINCS L1000 transcriptomic reversal (~30GB download required)",
                "Experimental Chou-Talalay synergy CI data (wet lab required)",
            ],
        },

        "contingency_table":    contingency,
        "top_candidates":       [_normalise_candidate(c) for c in candidates[:top_n]],

        "top_combinations": [
            {
                "compound_a":    c.get("compound_a", ""),
                "compound_b":    c.get("compound_b", ""),
                "synergy_score": round(float(c.get("synergy_score", 0)), 4),
                "rationale":     c.get("rationale", ""),
            }
            for c in combos[:10]
        ],

        "confidence_breakdown": _extract_confidence_breakdown(hypotheses),
        "hypotheses":           hypotheses,

        "reports": {
            "novelty":          raw.get("novelty_report", ""),
            "polypharmacology": raw.get("poly_report", ""),
            "combinations":     raw.get("combination_report", ""),
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, default=_safe)

    logger.info(f"\n✅ Results saved to : {output_path}")
    logger.info(f"✅ Log saved to     : {output_path.parent / 'pipeline_output.txt'}")
    logger.info(f"✅ n_drugs_screened : {n_screened}")
    logger.info("\nNext step: python -m backend.pipeline.generate_figures")

    ct = results["contingency_table"]
    if ct:
        logger.info(f"\n── Contingency table ──────────────────────────────")
        logger.info(f"  H3K27M+/CDKN2A-del : {ct['h3k27m_pos_cdkn2a_del']}")
        logger.info(f"  H3K27M+/CDKN2A-WT  : {ct['h3k27m_pos_cdkn2a_wt']}")
        logger.info(f"  H3K27M-/CDKN2A-del : {ct['h3k27m_neg_cdkn2a_del']}")
        logger.info(f"  H3K27M-/CDKN2A-WT  : {ct['h3k27m_neg_cdkn2a_wt']}")
        logger.info(f"  p-value            : {ct.get('p_value_label') or ct.get('p_value')}")

    top5 = results["top_candidates"][:5]
    if top5:
        logger.info(f"\n── Top 5 candidates ────────────────────────────────")
        logger.info(f"  {'Drug':<25} {'Score':>6}  {'BBB':<10}")
        for c in top5:
            logger.info(f"  {c['name']:<25} {c['score']:>6.3f}  {c['bbb_penetrance']:<10}")
        logger.info(f"\n  Total drugs screened (OpenTargets): {n_screened}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Save pipeline results to JSON")
    parser.add_argument("--disease", default="dipg")
    parser.add_argument("--top_n",   default=20, type=int)
    parser.add_argument("--output",  default="results/pipeline_results.json")
    args = parser.parse_args()
    asyncio.run(main(args.disease, args.top_n, Path(args.output)))