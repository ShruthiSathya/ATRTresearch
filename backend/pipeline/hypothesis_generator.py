"""
hypothesis_generator.py
========================
ATRT Drug Repurposing Pipeline — Hypothesis Generator v1.0

Key change from DIPG pipeline:
  - EZH2 inhibitors are flagged as BOOSTED (not penalised) in hypothesis notes
  - Confidence breakdown uses ATRT-aware escape bypass mode labels
  - No H3K27M or CDKN2A references
  - SMARCB1 statistics replace H3K27M/CDKN2A co-occurrence statistics
  - Confidence reported as range [conservative, optimistic]

Sources supporting confidence formula:
  Knutson 2013 PNAS (EZH2 synthetic lethality)
  Monje 2023 Nature Medicine (PBTC-047 toxicity precedent)
  Behan 2019 Nature (DepMap weight justification)
  Fischer 1998 J Med Chem (BBB penetrance)
"""

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

    Returns both conservative and optimistic confidence bounds.
    Conservative = additive toxicity model.
    Optimistic = dose-optimised estimate (60% of additive rate).

    The PBTC-047 precedent (panobinostat at 27.6% DLT with dose modification
    in pediatric DIPG) establishes that dose optimisation is clinically
    feasible — Monje et al. 2023 PMID 37526549.
    """
    if not top_3:
        return {"confidence": 0.0, "confidence_raw": 0.0, "explanation": "No candidates"}

    cw = CONFIDENCE_WEIGHTS

    # Component A: DepMap CRISPR essentiality
    # Source: Behan et al. 2019 Nature — justifies DepMap as primary signal
    depmap_scores    = [c.get("depmap_score", HYP_CONFIG["missing_depmap_score"]) for c in top_3]
    depmap_component = sum(depmap_scores) / len(depmap_scores)

    default_score = HYP_CONFIG["depmap_default_score"]
    tolerance     = HYP_CONFIG["depmap_default_tolerance"]
    all_default   = all(abs(s - default_score) < tolerance for s in depmap_scores)
    if all_default:
        depmap_component = DEPMAP_MISSING_PRIOR
        depmap_note = "DepMap data not loaded — using prior 0.30"
    else:
        depmap_note = (
            f"Broad CRISPR Chronos scores (ATRT/rhabdoid lines): "
            f"{[round(s, 2) for s in depmap_scores]}"
        )

    # Component B: BBB penetrance
    # Source: Fischer 1998 J Med Chem; Pardridge 2003 Mol Interv
    penetrance_scores = BBB_CONFIG["penetrance_scores"]
    bbb_cats      = [c.get("bbb_penetrance", "UNKNOWN") for c in top_3]
    bbb_scores    = [penetrance_scores.get(cat, penetrance_scores["UNKNOWN"]) for cat in bbb_cats]
    bbb_component = sum(bbb_scores) / len(bbb_scores)
    bbb_note = (
        f"BBB penetrance: "
        f"{ list(zip([c.get('name', '?') for c in top_3], bbb_cats)) }"
    )

    # Component C: Mechanistic diversity (target Jaccard)
    # Diverse target combinations reduce likelihood of cross-resistance
    all_target_sets   = [set(c.get("targets", [])) for c in top_3]
    pairwise_overlaps = [
        len(a & b) / max(len(a | b), 1)
        for i, a in enumerate(all_target_sets)
        for j, b in enumerate(all_target_sets)
        if i < j and (a or b)
    ]
    avg_overlap         = sum(pairwise_overlaps) / max(len(pairwise_overlaps), 1)
    diversity_component = 1.0 - avg_overlap
    diversity_note      = (
        f"Target Jaccard overlap: {round(avg_overlap, 3)} (lower = more diverse)"
    )

    # Raw confidence (pre-toxicity)
    confidence_raw = round(min(1.0, max(0.0,
        depmap_component    * cw["depmap"]
        + bbb_component     * cw["bbb"]
        + diversity_component * cw["diversity"]
    )), 4)

    # Toxicity penalty
    drug_names = [c.get("drug_name", c.get("name", "UNKNOWN")).upper() for c in top_3]
    tox_result = combination_toxicity_penalty(drug_names)

    multiplier     = tox_result["multiplier"]
    multiplier_opt = tox_result["multiplier_optimistic"]

    confidence_final      = round(confidence_raw * multiplier, 4)
    confidence_optimistic = round(confidence_raw * multiplier_opt, 4)

    return {
        "confidence":               confidence_final,
        "confidence_optimistic":    confidence_optimistic,
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
    ATRT hypothesis assembly with externally-grounded, toxicity-adjusted
    confidence scoring. Reports confidence as a RANGE not a point estimate.

    EZH2 logic is INVERTED vs DIPG:
      DIPG: EZH2 inhibitors in combo are a WARNING (H3K27M suppresses EZH2 — non-rational)
      ATRT: EZH2 inhibitors in combo are EXPECTED and POSITIVE (SMARCB1 loss synthetic lethality)
    """

    def __init__(self):
        logger.info("HypothesisGenerator v1.0 (ATRT — SMARCB1 synthetic lethality logic)")

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

        # Select top 3 with non-overlapping targets for maximum combination diversity
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

        # p-value handling — ATRT uses SMARCB1 loss statistics, not co-occurrence
        p_value_is_valid = (
            p_value is not None
            and not math.isnan(p_value)
            and math.isfinite(p_value)
        )

        if p_value is None:
            p_str    = "N/A — SMARCB1 loss is defining event in ATRT (no co-occurrence test needed)"
            priority = "COMPUTATIONAL"
        elif math.isnan(p_value):
            p_str    = "N/A — sample counts insufficient"
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

        # ATRT-specific: flag EZH2 inhibitors in combo as a POSITIVE feature
        # (opposite of DIPG where they were a warning)
        ezh2_in_combo = [
            c.get("name", "?") for c in top_3
            if c.get("atrt_components", {}).get("ezh2_boosted")
            or c.get("ezh2_boosted", False)
        ]

        aurka_in_combo = [
            c.get("name", "?") for c in top_3
            if c.get("atrt_components", {}).get("is_aurka_inhibitor")
            or c.get("aurka_boosted", False)
        ]

        # Escape bypass mode label
        escape_mode = "RNA-confirmed" if (
            genomic_stats and genomic_stats.get("has_rna_data")
        ) else "curated fallback"

        # SMARCB1 statistics note
        smarcb1_count = (genomic_stats or {}).get("smarcb1_loss_count", 0)
        total_samples = (genomic_stats or {}).get("total_samples", 0)
        smarcb1_note  = (
            f"SMARCB1 biallelic loss: {smarcb1_count}/{total_samples} samples "
            f"({smarcb1_count / max(total_samples, 1):.0%}). "
            "This is the defining molecular event in ATRT — no co-occurrence test required."
            if total_samples > 0
            else "SMARCB1 loss: estimated ~95% prevalence (Hasselblatt 2011, PMID 20625942)"
        )

        triple_hit = {
            "drug_or_combo":   combo_name,
            "priority":        priority,
            "confidence":      confidence_data["confidence"],
            "confidence_optimistic": confidence_data["confidence_optimistic"],
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
                    "(1) DepMap CRISPR Chronos essentiality in ATRT/rhabdoid lines (Broad Institute), "
                    "(2) BBB penetrance (curated PK literature), "
                    "(3) Target Jaccard diversity. "
                    "Confidence reported as range [conservative, optimistic] "
                    "accounting for dose-optimisation uncertainty (PBTC-047 precedent)."
                ),
                # ATRT-specific fields
                "escape_bypass_mode":       escape_mode,
                "ezh2_inhibitors_in_combo": ezh2_in_combo,
                "ezh2_combo_note": (
                    f"EZH2 inhibitor(s) in combo: {ezh2_in_combo} — "
                    "POSITIVE in ATRT (SMARCB1 loss synthetic lethality, Knutson 2013 PNAS). "
                    "This is the OPPOSITE of DIPG where EZH2 inhibitors are penalised."
                ) if ezh2_in_combo else "",
                "aurka_inhibitors_in_combo": aurka_in_combo,
                "ppi_weight_note": (
                    "PPI weight = 0.05: sparse curated neighbor coverage means low "
                    "discriminative signal for most candidates. DepMap weight = 0.35 compensates."
                ),
                "smarcb1_note": smarcb1_note,
            },
            "confidence_explanation": confidence_data["explanation"],
            "supporting_streams":     ["Multi-Omic Integration (ATRT)"],
            "target_context":         f"Multi-node blockade targeting {target_str}",
            "mechanism_narrative": (
                "Computationally derived synergistic combination based on SMARCB1-loss "
                "dependencies. Selected for mechanism diversity, CRISPR essentiality in "
                "ATRT/rhabdoid cell lines, and CNS penetrance."
            ),
            "statistical_significance": p_str,
            "statistical_note":        smarcb1_note,
            "bypass_status": (
                "HIGH"
                if all(c.get("escape_bypass_score", 0) > HYP_CONFIG["bypass_high_threshold"]
                       for c in top_3)
                else "MODERATE"
            ),
        }

        if genomic_stats and genomic_stats.get("smarcb1_loss_count"):
            triple_hit["smarcb1_cohort_note"] = (
                f"{genomic_stats['smarcb1_loss_count']} SMARCB1-null samples "
                f"({genomic_stats['smarcb1_loss_count'] / max(genomic_stats.get('total_samples', 1), 1):.0%} "
                "of cohort). All are candidates for EZH2-inhibitor synthetic lethality."
            )

        # Subgroup note if available
        if genomic_stats and genomic_stats.get("subgroup_counts"):
            sc = genomic_stats["subgroup_counts"]
            triple_hit["subgroup_note"] = (
                f"Subgroup distribution: TYR={sc.get('TYR', 'N/A')}, "
                f"SHH={sc.get('SHH', 'N/A')}, MYC={sc.get('MYC', 'N/A')} "
                "(Johann 2016, Cancer Cell)"
            )

        return [triple_hit]

    def generate_report(self, hypotheses: List[Dict]) -> str:
        lines = ["# ATRT Drug Repurposing Report v1.0\n"]
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
            ezh2_note  = bd.get("ezh2_combo_note", "")

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
            ]
            if ezh2_note:
                lines.append(f"  - EZH2 note: {ezh2_note}")
            lines += [
                f"- **Escape bypass mode:** {bd.get('escape_bypass_mode', 'unknown')}",
                f"- **Targets:** {h['target_context']}",
                f"- **SMARCB1 statistics:** {h.get('statistical_note', 'N/A')}",
                f"- **Mechanism:** {h['mechanism_narrative']}\n",
            ]
            if h.get("subgroup_note"):
                lines.append(f"- **Subgroup distribution:** {h['subgroup_note']}\n")

        return "\n".join(lines)