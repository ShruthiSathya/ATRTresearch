import logging
import math
from typing import Dict, List, Optional

from .toxicity_constraint import combination_toxicity_penalty
from .pipeline_config import (
    CONFIDENCE_WEIGHTS,
    DEPMAP_MISSING_PRIOR,
    BBB as BBB_CONFIG,
    HYPOTHESIS as HYP_CONFIG,
)

logger = logging.getLogger(__name__)


def _compute_externally_grounded_confidence(top_3: List[Dict]) -> Dict:
    """
    Compute confidence from three externally-grounded signals,
    then apply toxicity penalty multiplier.

    v5.5: Returns both conservative and optimistic confidence bounds,
    making explicit that the adjusted confidence is a RANGE not a point.
    """
    if not top_3:
        return {"confidence": 0.0, "confidence_raw": 0.0, "explanation": "No candidates"}

    cw = CONFIDENCE_WEIGHTS

    # ── Component A: DepMap CRISPR essentiality ───────────────────────────────
    depmap_scores    = [c.get("depmap_score", HYP_CONFIG["missing_depmap_score"]) for c in top_3]
    depmap_component = sum(depmap_scores) / len(depmap_scores)

    default_score = HYP_CONFIG["depmap_default_score"]
    tolerance     = HYP_CONFIG["depmap_default_tolerance"]
    all_default   = all(abs(s - default_score) < tolerance for s in depmap_scores)
    if all_default:
        depmap_component = DEPMAP_MISSING_PRIOR
        depmap_note = "DepMap data not loaded — using prior 0.30"
    else:
        depmap_note = f"Broad CRISPR Chronos scores: {[round(s, 2) for s in depmap_scores]}"

    # ── Component B: BBB penetrance ───────────────────────────────────────────
    penetrance_scores = BBB_CONFIG["penetrance_scores"]
    bbb_cats      = [c.get("bbb_penetrance", "UNKNOWN") for c in top_3]
    bbb_scores    = [penetrance_scores.get(cat, penetrance_scores["UNKNOWN"]) for cat in bbb_cats]
    bbb_component = sum(bbb_scores) / len(bbb_scores)
    bbb_note = (
        f"BBB penetrance: "
        f"{ list(zip([c.get('name', '?') for c in top_3], bbb_cats)) }"
    )

    # ── Component C: Mechanistic diversity ────────────────────────────────────
    all_target_sets   = [set(c.get("targets", [])) for c in top_3]
    pairwise_overlaps = [
        len(a & b) / max(len(a | b), 1)
        for i, a in enumerate(all_target_sets)
        for j, b in enumerate(all_target_sets)
        if i < j and (a or b)
    ]
    avg_overlap         = sum(pairwise_overlaps) / max(len(pairwise_overlaps), 1)
    diversity_component = 1.0 - avg_overlap
    diversity_note      = f"Target Jaccard overlap: {round(avg_overlap, 3)} (lower = more diverse)"

    # ── Raw confidence (pre-toxicity) ─────────────────────────────────────────
    confidence_raw = round(min(1.0, max(0.0,
        depmap_component    * cw["depmap"]
        + bbb_component     * cw["bbb"]
        + diversity_component * cw["diversity"]
    )), 4)

    # ── Toxicity penalty (single call — not duplicated in discovery_pipeline) ─
    drug_names = [c.get("drug_name", c.get("name", "UNKNOWN")).upper() for c in top_3]
    tox_result = combination_toxicity_penalty(drug_names)

    multiplier     = tox_result["multiplier"]
    multiplier_opt = tox_result["multiplier_optimistic"]

    confidence_final     = round(confidence_raw * multiplier, 4)
    confidence_optimistic = round(confidence_raw * multiplier_opt, 4)

    return {
        "confidence":               confidence_final,       # conservative (lower bound)
        "confidence_optimistic":    confidence_optimistic,  # dose-optimised (upper bound)
        "confidence_raw":           confidence_raw,
        "depmap_component":         round(depmap_component, 4),
        "bbb_component":            round(bbb_component, 4),
        "diversity_component":      round(diversity_component, 4),
        "toxicity_multiplier":      multiplier,
        "toxicity_multiplier_optimistic": multiplier_opt,
        "toxicity_flag":            tox_result["flag"],
        "toxicity_note":            tox_result["note"],
        "toxicity_confidence_note": tox_result["confidence_note"],
        "depmap_note":              depmap_note,
        "bbb_note":                 bbb_note,
        "diversity_note":           diversity_note,
        "explanation": (
            f"Raw = {cw['depmap']}×DepMap({depmap_component:.2f}) "
            f"+ {cw['bbb']}×BBB({bbb_component:.2f}) "
            f"+ {cw['diversity']}×Diversity({diversity_component:.2f}) = {confidence_raw:.2f}. "
            f"Toxicity multiplier: {multiplier:.3f} [conservative] / {multiplier_opt:.3f} [optimistic] "
            f"({tox_result['flag']}). "
            f"Adjusted confidence RANGE: [{confidence_final:.2f}, {confidence_optimistic:.2f}]. "
            f"{tox_result['note']}"
        ),
    }


class HypothesisGenerator:
    """
    v5.5: Hypothesis assembly with externally-grounded, toxicity-adjusted
    confidence scoring. Reports confidence as a RANGE not a point estimate.
    """

    def __init__(self):
        logger.info("✅ Hypothesis Generator v5.5 (Confidence range reporting)")

    def generate(
        self,
        candidates:        List[Dict],
        cmap_results:      List[Dict],
        synergy_combos:    List[Dict],
        differential_cmap: List[Dict],
        genomic_stats:     Optional[Dict] = None,
        p_value:           Optional[float] = None,
    ) -> List[Dict]:

        sorted_candidates = sorted(candidates, key=lambda x: x.get("score", 0), reverse=True)
        if len(sorted_candidates) < 3:
            return []

        # Select top 3 with non-overlapping targets
        top_3:           List[Dict] = []
        covered_targets: set        = set()

        for c in sorted_candidates:
            c_targets = set(c.get("targets", []))
            if not c_targets.intersection(covered_targets):
                top_3.append(c)
                covered_targets.update(c_targets)
            if len(top_3) == HYP_CONFIG["combo_size"]:
                break

        if len(top_3) < 3:
            top_3 = sorted_candidates[:3]

        confidence_data = _compute_externally_grounded_confidence(top_3)

        # p-value handling
        p_value_is_valid = (
            p_value is not None
            and not math.isnan(p_value)
            and math.isfinite(p_value)
        )

        if p_value is None:
            p_str    = "N/A — genomic validation data not loaded"
            priority = "COMPUTATIONAL"
        elif math.isnan(p_value):
            p_str    = "N/A — sample counts insufficient for Fisher's exact test"
            priority = "COMPUTATIONAL"
        elif p_value < HYP_CONFIG["p_value_significance"]:
            p_str    = f"{p_value:.2e} ✅"
            priority = "HIGH"
        else:
            p_str    = f"{p_value:.4f} (not significant)"
            priority = "MODERATE"

        combo_name = " + ".join([
            c.get("drug_name", c.get("name", "Unknown")).upper() for c in top_3
        ])
        combo_targets = []
        for c in top_3:
            combo_targets.extend(c.get("targets", []))
        target_str = " / ".join(list(dict.fromkeys(combo_targets))[:5])

        # FIX v5.5: Report EZH2 inhibitor penalties in hypothesis if any top-3 were penalised
        ezh2_warnings = [
            c.get("name", "?") for c in top_3
            if c.get("dipg_components", {}).get("is_ezh2_inhibitor")
        ]

        # FIX v5.5: Note escape bypass mode (RNA-confirmed vs curated fallback)
        escape_mode = "RNA-confirmed" if (
            genomic_stats and genomic_stats.get("has_rna_data")
        ) else "curated fallback"

        triple_hit = {
            "drug_or_combo":   combo_name,
            "priority":        priority,
            "confidence":      confidence_data["confidence"],             # conservative lower bound
            "confidence_optimistic": confidence_data["confidence_optimistic"],  # upper bound
            "confidence_range": (
                f"[{confidence_data['confidence']:.2f}, "
                f"{confidence_data['confidence_optimistic']:.2f}]"
            ),
            "confidence_breakdown": {
                "confidence_raw":             confidence_data["confidence_raw"],
                "confidence_adjusted":        confidence_data["confidence"],
                "confidence_adjusted_optimistic": confidence_data["confidence_optimistic"],
                "confidence_range":           (
                    f"[{confidence_data['confidence']:.2f}, "
                    f"{confidence_data['confidence_optimistic']:.2f}]"
                ),
                "depmap_essentiality":        confidence_data["depmap_component"],
                "bbb_penetrance":             confidence_data["bbb_component"],
                "mechanistic_diversity":      confidence_data["diversity_component"],
                "toxicity_multiplier":        confidence_data["toxicity_multiplier"],
                "toxicity_multiplier_optimistic": confidence_data["toxicity_multiplier_optimistic"],
                "toxicity_flag":              confidence_data["toxicity_flag"],
                "toxicity_note":              confidence_data["toxicity_note"],
                "toxicity_confidence_note":   confidence_data["toxicity_confidence_note"],
                "weights_used": {
                    "depmap":    CONFIDENCE_WEIGHTS["depmap"],
                    "bbb":       CONFIDENCE_WEIGHTS["bbb"],
                    "diversity": CONFIDENCE_WEIGHTS["diversity"],
                },
                "method": (
                    "Weighted combination: "
                    "(1) DepMap CRISPR Chronos essentiality (Broad Institute), "
                    "(2) BBB penetrance (curated PK literature), "
                    "(3) Target Jaccard diversity. "
                    "v5.5: confidence reported as range [conservative, optimistic] "
                    "accounting for dose-optimisation uncertainty."
                ),
                # FIX v5.5: surface pipeline accuracy notes
                "escape_bypass_mode": escape_mode,
                "ezh2_inhibitors_in_combo": ezh2_warnings,
                "ppi_weight_note": (
                    "PPI weight reduced to 0.05 (from 0.10) in v5.5: "
                    "490/557 drugs scored at floor (0.20) due to sparse curated neighbors. "
                    "DepMap weight increased to 0.35 to compensate."
                ),
            },
            "confidence_explanation": confidence_data["explanation"],
            "supporting_streams":     ["Multi-Omic Integration"],
            "target_context":         f"Multi-node blockade targeting {target_str}",
            "mechanism_narrative": (
                "Computationally derived synergistic combination. "
                "Selected for mechanism diversity, network proximity, and stem-cell eradication."
            ),
            "statistical_significance": p_str,
            "statistical_note": (
                "Fisher's exact test for H3K27M/CDKN2A-del co-occurrence in CBTN cohort."
                if p_value_is_valid
                else "Requires data/validation/cbtn_genomics/ files to be populated."
            ),
            "bypass_status": (
                "HIGH"
                if all(c.get("escape_bypass_score", 0) > HYP_CONFIG["bypass_high_threshold"]
                       for c in top_3)
                else "MODERATE"
            ),
        }

        if genomic_stats and genomic_stats.get("overlap_count"):
            triple_hit["patient_population"] = (
                f"{genomic_stats['overlap_count']} specific samples "
                f"({genomic_stats['prevalence']:.1%} prevalence)"
            )

        if genomic_stats and genomic_stats.get("acvr1_estimated_n"):
            triple_hit["acvr1_subgroup_note"] = (
                f"~{genomic_stats['acvr1_estimated_n']} samples estimated ACVR1-mutant "
                f"(~25% of H3K27M+). ACVR1-relevant candidates in separate subgroup report."
            )

        return [triple_hit]

    def generate_report(self, hypotheses: List[Dict]) -> str:
        lines = ["# GBM/DIPG Unbiased Discovery Report v5.5\n"]
        for h in hypotheses:
            bd       = h.get("confidence_breakdown", {})
            raw      = bd.get("confidence_raw", 0)
            adj      = bd.get("confidence_adjusted", 0)
            adj_opt  = bd.get("confidence_adjusted_optimistic", adj)
            mult     = bd.get("toxicity_multiplier", 1)
            mult_opt = bd.get("toxicity_multiplier_optimistic", 1)
            flag     = bd.get("toxicity_flag", "?")
            wts      = bd.get("weights_used", CONFIDENCE_WEIGHTS)
            conf_range = bd.get("confidence_range", f"[{adj:.2f}, {adj_opt:.2f}]")

            lines += [
                f"## {h['drug_or_combo']}",
                f"- **Priority:** {h['priority']}",
                f"- **Confidence range (adjusted):** {conf_range}",
                f"  - Conservative (additive tox model): {adj:.2f}  "
                f"*(raw: {raw:.2f} × tox: {mult:.3f} — {flag})*",
                f"  - Optimistic (dose-optimised):       {adj_opt:.2f}  "
                f"*(raw: {raw:.2f} × tox: {mult_opt:.3f})*",
                f"  - DepMap essentiality (w={wts.get('depmap', '?')}): {bd.get('depmap_essentiality', 0):.2f}",
                f"  - BBB penetrance     (w={wts.get('bbb', '?')}): {bd.get('bbb_penetrance', 0):.2f}",
                f"  - Target diversity   (w={wts.get('diversity', '?')}): {bd.get('mechanistic_diversity', 0):.2f}",
                f"  - Toxicity note: {bd.get('toxicity_note', '')}",
                f"  - {bd.get('toxicity_confidence_note', '')}",
                f"- **Escape bypass mode:** {bd.get('escape_bypass_mode', 'unknown')}",
                f"- **Targets:** {h['target_context']}",
                f"- **Statistical Significance:** {h.get('statistical_significance', 'N/A')}",
                f"- **Mechanism:** {h['mechanism_narrative']}\n",
            ]

            if h.get("acvr1_subgroup_note"):
                lines.append(f"- **ACVR1 subgroup:** {h['acvr1_subgroup_note']}\n")

        return "\n".join(lines)