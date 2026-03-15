import asyncio
import json
import logging
import argparse
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def _is_dipg_or_gbm(disease_name: str) -> bool:
    keywords = (
        "dipg", "diffuse intrinsic pontine glioma",
        "glioblastoma", "gbm", "h3k27m", "high-grade glioma",
        "diffuse midline glioma",
    )
    return any(k in disease_name.lower() for k in keywords)


async def run_dipg_pipeline(
    disease_name:         str   = "dipg",
    exclude_low_bbb:      bool  = False,
    apply_bbb_penalty:    bool  = True,
    top_n:                int   = 20,
    predict_combinations: bool  = True,
    run_quantile_sensitivity: bool = True,
) -> Dict:
    """
    Run the full DIPG/GBM scoring stack.

    v5.5 additions:
    - GSC quantile sensitivity analysis (p75/p80/p85/p90)
    - ACVR1 subgroup stratification report
    - EZH2 inhibitor penalty reporting
    - Toxicity confidence range [conservative, optimistic]
    - RNA-confirmed vs curated escape bypass mode surfaced in stats
    """
    from backend.pipeline.discovery_pipeline import ProductionPipeline
    from backend.pipeline.dipg_specialization import (
        DIPGSpecializedScorer,
        DIPG_CORE_GENES,
        get_dipg_disease_data_supplement,
    )
    from backend.pipeline.polypharmacology import PolypharmacologyScorer
    from backend.pipeline.synergy_predictor import SynergyPredictor

    is_dipg = _is_dipg_or_gbm(disease_name)

    logger.info("=" * 70)
    logger.info("DIPG/GBM Pipeline v5.5 — %s", disease_name)
    logger.info("=" * 70)

    # ── Step 1: Run base pipeline ─────────────────────────────────────────────
    logger.info("[1/6] Running base ProductionPipeline...")
    pipeline = ProductionPipeline()
    await pipeline.initialize(disease=disease_name)
    results = await pipeline.run(disease_name=disease_name, top_k=top_n)

    hypotheses = results.get("hypotheses", [])
    stats      = results.get("stats", {})

    logger.info(
        "      Base pipeline complete — %d hypotheses, p=%s, escape mode=%s",
        len(hypotheses),
        stats.get("p_value_label", "N/A"),
        stats.get("escape_bypass_mode", "unknown"),
    )

    # ── Step 2: Fetch candidates ───────────────────────────────────────────────
    logger.info("[2/6] Fetching candidates for DIPG scoring...")
    candidates = await pipeline._data_fetcher.fetch_approved_drugs()

    if not candidates:
        logger.warning("No candidates returned — using fallback library")
        candidates = [
            {"name": "ONC201",       "targets": ["DRD2", "CLPB"]},
            {"name": "Panobinostat", "targets": ["HDAC1", "HDAC2"]},
            {"name": "Abemaciclib",  "targets": ["CDK4", "CDK6"]},
            {"name": "Marizomib",    "targets": ["PSMB5", "PSMB2"]},
            {"name": "Tazemetostat", "targets": ["EZH2"]},
        ]

    # ── Step 3: DIPG specialization ───────────────────────────────────────────
    if is_dipg:
        logger.info("[3/6] Applying DIPG H3K27M/ACVR1 specialization v5.5...")
        dipg_scorer = DIPGSpecializedScorer(
            apply_bbb_penalty=apply_bbb_penalty,
            novelty_bonus=0.08,
            h3k27m_bonus=0.12,
            acvr1_bonus=0.10,
        )
        candidates = dipg_scorer.score_batch(candidates)

        novel_candidates = [
            c for c in candidates
            if c.get("dipg_components", {}).get("is_untested_dipg")
            and c.get("score", 0) > 0.40
        ]
        ezh2_penalised = [
            c for c in candidates
            if c.get("dipg_components", {}).get("is_ezh2_inhibitor")
        ]
        novelty_report    = dipg_scorer.generate_novelty_report(candidates, top_n=10)
        acvr1_report      = dipg_scorer.generate_acvr1_subgroup_report(candidates)

        n_h3k27m  = sum(1 for c in candidates if c.get("dipg_components", {}).get("h3k27m_relevant"))
        n_novel   = len(novel_candidates)
        n_ezh2pen = len(ezh2_penalised)

        logger.info(
            "      H3K27M-relevant: %d | Novel: %d | EZH2 inhibitors penalised: %d",
            n_h3k27m, n_novel, n_ezh2pen,
        )
        if ezh2_penalised:
            logger.info(
                "      EZH2 penalty applied to: %s",
                [c.get("name", "?") for c in ezh2_penalised],
            )
    else:
        novel_candidates  = []
        novelty_report    = "Non-DIPG disease — novelty report not generated."
        acvr1_report      = "Non-DIPG disease — ACVR1 subgroup report not generated."
        ezh2_penalised    = []
        logger.info("[3/6] Non-DIPG disease — skipping H3K27M specialization")

    # ── Step 4: GSC quantile sensitivity analysis ──────────────────────────────
    quantile_sensitivity = {}
    if run_quantile_sensitivity and is_dipg:
        logger.info("[4/6] Running GSC quantile sensitivity analysis (p75/p80/p85/p90)...")
        quantile_sensitivity = await pipeline._tissue.run_quantile_sensitivity(candidates)
        if quantile_sensitivity.get("stable") is not None:
            stability = "STABLE ✅" if quantile_sensitivity["stable"] else "UNSTABLE ⚠️"
            logger.info("      GSC quantile sensitivity: %s", stability)
            logger.info("      %s", quantile_sensitivity.get("note", ""))
    else:
        logger.info("[4/6] Quantile sensitivity skipped (non-DIPG or disabled)")

    # ── Step 5: Polypharmacology scoring ──────────────────────────────────────
    logger.info("[5/6] Polypharmacology scoring...")
    disease_genes = DIPG_CORE_GENES if is_dipg else ["EGFR", "PTEN", "TP53", "CDK4"]

    poly_scorer = PolypharmacologyScorer(disease=disease_name, dipg_mode=is_dipg)
    top_50      = sorted(candidates, key=lambda c: c.get("score", 0), reverse=True)[:50]
    top_50      = poly_scorer.score_batch(top_50, disease_targets=disease_genes)
    poly_report = poly_scorer.generate_poly_report(top_50)

    # ── Step 6: Combination synergy ───────────────────────────────────────────
    top_combinations  = []
    combination_report = (
        "⚠️  Synergy module: biological prior estimates (Grasso 2015). "
        "Experimental CI validation required.\n"
    )

    if predict_combinations:
        logger.info("[6/6] Drug combination synergy prediction...")
        syn_predictor    = SynergyPredictor()
        top_combinations = syn_predictor.predict_top_combinations(top_50[:30])
        if top_combinations:
            combination_report = (
                "⚠️  NOTE: Synergy scores are biological priors (Grasso 2015), "
                "not experimentally validated combination indices.\n\n"
                + "\n".join(
                    f"  {c.get('compound_a','?')} + {c.get('compound_b','?')}: "
                    f"score={c.get('synergy_score',0):.2f} | {c.get('rationale','')}"
                    for c in top_combinations[:5]
                )
            )
    else:
        logger.info("[6/6] Combination prediction skipped")

    # ── Final summary ─────────────────────────────────────────────────────────
    top_candidates = sorted(candidates, key=lambda c: c.get("score", 0), reverse=True)[:top_n]

    # IC50 validation report
    try:
        from backend.pipeline.published_ic50_validation import (
            annotate_candidates_with_ic50,
            generate_ic50_validation_report,
        )
        top_candidates = annotate_candidates_with_ic50(top_candidates)
        ic50_report    = generate_ic50_validation_report(top_candidates)
    except Exception as e:
        ic50_report = f"IC50 validation report unavailable: {e}"

    logger.info("\n" + "=" * 70)
    logger.info("DIPG v5.5 RESULTS — Top 10")
    logger.info("%-30s %6s %8s %8s %8s", "Drug", "Score", "H3K27M", "Novel", "EZH2-pen")
    logger.info("-" * 65)
    for c in top_candidates[:10]:
        name    = (c.get("name") or "?")[:29]
        score   = c.get("score", 0)
        h3k27m  = "Y" if c.get("dipg_components", {}).get("h3k27m_relevant") else "-"
        novel   = "Y" if c.get("dipg_components", {}).get("is_untested_dipg") else "-"
        ezh2pen = "⚠" if c.get("dipg_components", {}).get("is_ezh2_inhibitor") else "-"
        logger.info("%-30s %6.3f %8s %8s %8s", name, score, h3k27m, novel, ezh2pen)

    return {
        "hypotheses":          hypotheses,
        "stats":               {
            **stats,
            "quantile_sensitivity": quantile_sensitivity,
            "ezh2_penalised_count": len(ezh2_penalised),
        },
        "top_candidates":      top_candidates,
        "novel_candidates":    novel_candidates,
        "top_combinations":    top_combinations,
        "novelty_report":      novelty_report,
        "poly_report":         poly_report,
        "combination_report":  combination_report,
        "acvr1_subgroup_report": acvr1_report,
        "ic50_validation_report": ic50_report,
        "pipeline_stats": {
            "total_candidates":    len(candidates),
            "dipg_specialised":    is_dipg,
            "novel_count":         len(novel_candidates),
            "ezh2_penalised":      [c.get("name", "?") for c in ezh2_penalised],
            "combinations_found":  len(top_combinations),
            "escape_bypass_mode":  stats.get("escape_bypass_mode", "unknown"),
            "rna_upregulated_genes": stats.get("rna_upregulated_genes", 0),
            "quantile_stable":     quantile_sensitivity.get("stable"),
            "data_streams_active": [
                "DepMap CRISPR (Broad Institute) [w=0.35, up from 0.30]",
                "Single-cell RNA-seq (GSE102130, Filbin 2018, H3K27M DIPG) [w=0.40]",
                "OpenTargets API",
                "STRING-DB PPI [w=0.05, down from 0.10 — sparse coverage]",
                "PedcBioPortal genomic validation (PNOC/PBTA, n=184)",
                "RNA-confirmed escape bypass (GSE115397 n=5 — note small N)",
            ],
            "v5_5_fixes_applied": [
                "EZH2 inhibitor penalty (0.50×) — consistent with tissue scorer",
                "Marizomib PPI score fixed (PSMB5 neighbors added)",
                "UNKNOWN BBB score: 0.5 → 0.4 (unknown ≠ moderate)",
                "AZD-8055, crizotinib, pazopanib etc. BBB filled from PK literature",
                "PPI weight: 0.10 → 0.05; DepMap weight: 0.30 → 0.35",
                "Escape bypass: RNA-confirmed when n≥10 upregulated genes",
                "Toxicity confidence: range [conservative, optimistic] not point estimate",
                "GSC quantile sensitivity: p75/p80/p85/p90 stability test",
                "ACVR1 subgroup stratification report",
            ],
        },
    }


async def _cli_main(disease: str, output: Optional[str], top_n: int, combinations: bool):
    result = await run_dipg_pipeline(
        disease_name=disease,
        top_n=top_n,
        predict_combinations=combinations,
    )

    print("\n" + "=" * 70)
    print("HYPOTHESES")
    print("=" * 70)
    from backend.pipeline.discovery_pipeline import ProductionPipeline
    tmp = ProductionPipeline()
    print(tmp._hyp_gen.generate_report(result["hypotheses"]))

    print("\n" + "=" * 70)
    print("IC50 VALIDATION (Published DIPG Cell Line Data)")
    print("=" * 70)
    print(result.get("ic50_validation_report", "Not available"))

    print("\n" + "=" * 70)
    print("ACVR1 SUBGROUP REPORT")
    print("=" * 70)
    print(result["acvr1_subgroup_report"])

    print("\n" + "=" * 70)
    print("NOVELTY REPORT")
    print("=" * 70)
    print(result["novelty_report"])

    print("\n" + "=" * 70)
    print("GSC QUANTILE SENSITIVITY")
    print("=" * 70)
    qs = result["stats"].get("quantile_sensitivity", {})
    print(qs.get("note", "Not run"))
    if qs.get("rankings_by_quantile"):
        for q, ranking in qs["rankings_by_quantile"].items():
            print(f"  p{float(q)*100:.0f}: {ranking[:5]}")

    print("\n" + "=" * 70)
    print("v5.5 FIXES APPLIED")
    print("=" * 70)
    for fix in result["pipeline_stats"].get("v5_5_fixes_applied", []):
        print(f"  ✅ {fix}")

    if output:
        with open(output, "w") as f:
            json.dump(result, f, indent=2, default=str)
        print(f"\nResults saved to: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GBM/DIPG Drug Repurposing Pipeline v5.5")
    parser.add_argument("--disease",      default="dipg")
    parser.add_argument("--output",       default=None)
    parser.add_argument("--top_n",        default=20, type=int)
    parser.add_argument("--combinations", action="store_true")
    parser.add_argument("--no-sensitivity", action="store_true",
                        help="Skip GSC quantile sensitivity analysis")
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )
    asyncio.run(_cli_main(args.disease, args.output, args.top_n, args.combinations))