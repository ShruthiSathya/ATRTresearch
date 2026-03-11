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

    Weights are read from pipeline_config.CONFIDENCE_WEIGHTS.
    No magic numbers in this function.
    """
    if not top_3:
        return {"confidence": 0.0, "confidence_raw": 0.0, "explanation": "No candidates"}

    cw = CONFIDENCE_WEIGHTS   # shorthand

    # ── Component A: DepMap CRISPR essentiality (Broad Institute) ─────────────
    depmap_scores    = [c.get("depmap_score", HYP_CONFIG["missing_depmap_score"]) for c in top_3]
    depmap_component = sum(depmap_scores) / len(depmap_scores)

    default_score = HYP_CONFIG["depmap_default_score"]
    tolerance     = HYP_CONFIG["depmap_default_tolerance"]
    all_default   = all(abs(s - default_score) < tolerance for s in depmap_scores)
    if all_default:
        depmap_component = DEPMAP_MISSING_PRIOR
        depmap_note = "DepMap data not loaded (CRISPRGeneEffect.csv missing) — using prior"
    else:
        depmap_note = f"Broad CRISPR Chronos scores: {[round(s, 2) for s in depmap_scores]}"

    # ── Component B: BBB penetrance (curated PK literature) ───────────────────
    penetrance_scores = BBB_CONFIG["penetrance_scores"]
    bbb_cats      = [c.get("bbb_penetrance", "UNKNOWN") for c in top_3]
    bbb_scores    = [penetrance_scores.get(cat, penetrance_scores["UNKNOWN"]) for cat in bbb_cats]
    bbb_component = sum(bbb_scores) / len(bbb_scores)
    bbb_note = (
        f"BBB penetrance: "
        f"{ list(zip([c.get('name', '?') for c in top_3], bbb_cats)) }"
    )

    # ── Component C: Mechanistic diversity (target Jaccard) ───────────────────
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
    multiplier = tox_result["multiplier"]

    confidence_final = round(confidence_raw * multiplier, 4)

    return {
        "confidence":          confidence_final,
        "confidence_raw":      confidence_raw,
        "depmap_component":    round(depmap_component, 4),
        "bbb_component":       round(bbb_component, 4),
        "diversity_component": round(diversity_component, 4),
        "toxicity_multiplier": multiplier,
        "toxicity_flag":       tox_result["flag"],
        "toxicity_note":       tox_result["note"],
        "depmap_note":         depmap_note,
        "bbb_note":            bbb_note,
        "diversity_note":      diversity_note,
        "explanation": (
            f"Raw = {cw['depmap']}×DepMap({depmap_component:.2f}) "
            f"+ {cw['bbb']}×BBB({bbb_component:.2f}) "
            f"+ {cw['diversity']}×Diversity({diversity_component:.2f}) = {confidence_raw:.2f}. "
            f"Toxicity multiplier = {multiplier:.3f} ({tox_result['flag']}). "
            f"Adjusted = {confidence_raw:.2f} × {multiplier:.3f} = {confidence_final:.2f}. "
            f"{tox_result['note']}"
        ),
    }


class HypothesisGenerator:
    """
    v5.4: Hypothesis assembly with externally-grounded, toxicity-adjusted
    confidence scoring. All weights read from pipeline_config.
    """

    def __init__(self):
        logger.info("✅ Hypothesis Generator Initialized (Diversity Mode)")

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

        triple_hit = {
            "drug_or_combo": combo_name,
            "priority":      priority,
            "confidence":    confidence_data["confidence"],
            "confidence_breakdown": {
                "confidence_raw":        confidence_data["confidence_raw"],
                "confidence_adjusted":   confidence_data["confidence"],
                "depmap_essentiality":   confidence_data["depmap_component"],
                "bbb_penetrance":        confidence_data["bbb_component"],
                "mechanistic_diversity": confidence_data["diversity_component"],
                "toxicity_multiplier":   confidence_data["toxicity_multiplier"],
                "toxicity_flag":         confidence_data["toxicity_flag"],
                "toxicity_note":         confidence_data["toxicity_note"],
                "weights_used": {
                    "depmap":    CONFIDENCE_WEIGHTS["depmap"],
                    "bbb":       CONFIDENCE_WEIGHTS["bbb"],
                    "diversity": CONFIDENCE_WEIGHTS["diversity"],
                },
                "method": (
                    "Weighted combination of: "
                    "(1) DepMap CRISPR Chronos essentiality (Broad Institute), "
                    "(2) BBB penetrance (curated PK literature), "
                    "(3) Target Jaccard diversity. "
                    "v5.4: multiplied by hematologic toxicity penalty from Phase I AE data."
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

        return [triple_hit]

    def generate_report(self, hypotheses: List[Dict]) -> str:
        lines = ["# GBM/DIPG Unbiased Discovery Report v5.4\n"]
        for h in hypotheses:
            bd   = h.get("confidence_breakdown", {})
            raw  = bd.get("confidence_raw", 0)
            adj  = bd.get("confidence_adjusted", 0)
            mult = bd.get("toxicity_multiplier", 1)
            flag = bd.get("toxicity_flag", "?")
            wts  = bd.get("weights_used", CONFIDENCE_WEIGHTS)

            lines += [
                f"## {h['drug_or_combo']}",
                f"- **Priority:** {h['priority']}",
                f"- **Confidence (adjusted):** {adj:.2f}  *(raw: {raw:.2f} × tox: {mult:.3f} — {flag})*",
                f"  - DepMap essentiality (w={wts.get('depmap', '?')}): {bd.get('depmap_essentiality', 0):.2f}",
                f"  - BBB penetrance     (w={wts.get('bbb', '?')}): {bd.get('bbb_penetrance', 0):.2f}",
                f"  - Target diversity   (w={wts.get('diversity', '?')}): {bd.get('mechanistic_diversity', 0):.2f}",
                f"  - Toxicity note: {bd.get('toxicity_note', '')}",
                f"- **Targets:** {h['target_context']}",
                f"- **Statistical Significance:** {h.get('statistical_significance', 'N/A')}",
                f"- **Mechanism:** {h['mechanism_narrative']}\n",
            ]
        return "\n".join(lines)