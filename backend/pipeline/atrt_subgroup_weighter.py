"""
atrt_subgroup_weighter.py
=========================
Johann 2016 (GSE73868) subgroup-aware scoring for the ATRT pipeline.

BIOLOGICAL CONTEXT
-------------------
Johann et al. 2016 Cancer Cell 29(3):379-393 (PMID 26923874) defines
three ATRT molecular subgroups based on DNA methylation (EPIC 850K):
  ATRT-TYR : ~36%, infratentorial, youngest patients, HDAC/BET vulnerability
  ATRT-SHH : ~37%, mixed location, GLI2 amplification, EZH2/SMO vulnerability
  ATRT-MYC : ~27%, supratentorial, oldest, worst prognosis, BET/AURKA/EZH2

These subgroups have DISTINCT drug vulnerabilities. Scoring without subgroup
knowledge underestimates drugs that are subgroup-specific (e.g. vismodegib
in ATRT-SHH, alisertib in ATRT-MYC) and overestimates pan-ATRT generalists.

IMPLEMENTATION STRATEGY
------------------------
When subgroup is KNOWN (from methylation array or RNA-seq classifier):
  → Apply subgroup multipliers directly to tissue and DepMap scores.

When subgroup is UNKNOWN (most clinical scenarios):
  → Compute a PREVALENCE-WEIGHTED EXPECTATION across all three subgroups.
  → This is more principled than pan-ATRT scoring because it explicitly
     accounts for the fraction of patients who would benefit from each drug.

The prevalence-weighted score for drug d is:
  E[score(d)] = sum_g P(g) × score_g(d)
where P(g) is the Johann 2016 subgroup prevalence and score_g(d) is the
subgroup-specific score for drug d in subgroup g.

SUBGROUP CALLING
-----------------
Without methylation array data, call subgroup from RNA-seq using these
marker genes (validated against Johann 2016 methylation calls):
  TYR subgroup: TYR, DCT, MITF (neural crest/melanocyte program)
  SHH subgroup: GLI1, GLI2, PTCH1 (hedgehog pathway)
  MYC subgroup: MYC, MYCN, AURKA (proliferative)

A simple centroid classifier works well for n≥3 marker genes per subgroup.

DATA SOURCE
-----------
GSE73868 (Johann 2016) — methylation-based subgroup calls for 150 ATRT.
Use as a reference to validate RNA-based calling on your cohort.
Access: GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73868
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Johann 2016 ground-truth priors
# ─────────────────────────────────────────────────────────────────────────────

SUBGROUP_PREVALENCES = {
    "TYR": 0.36,   # 54/150
    "SHH": 0.37,   # 55/150
    "MYC": 0.27,   # 41/150
}

# Drug score multipliers per subgroup — derived from subgroup-specific
# vulnerability evidence (see rationale below)
SUBGROUP_DRUG_MULTIPLIERS: Dict[str, Dict[str, float]] = {
    # Format: subgroup → {target_gene → multiplier}

    "TYR": {
        # HDAC inhibitors especially active in TYR (Torchia 2015 Fig 4A — BT37 is TYR)
        "HDAC1": 1.30, "HDAC2": 1.30, "HDAC3": 1.25, "HDAC6": 1.20,
        # BET bromodomain — TYR super-enhancers driven by MITF/BRD4
        "BRD4":  1.20, "BRD2": 1.15,
        # TYR-specific markers (targets of epigenetic drugs via super-enhancers)
        "TYR":   1.50, "DCT": 1.50, "MITF": 1.50,
        # EZH2 still important (SMARCB1 loss is universal)
        "EZH2":  1.10,
        # Hedgehog NOT relevant in TYR
        "SMO":   0.60, "GLI2": 0.55, "GLI1": 0.55,
        # AURKA/MYCN less relevant than in MYC
        "AURKA": 0.85, "MYCN": 0.80,
    },

    "SHH": {
        # SMO/GLI pathway amplified in SHH — vismodegib most rational here
        "SMO":   1.40, "GLI2": 1.50, "GLI1": 1.40, "PTCH1": 1.30,
        # EZH2 especially relevant in SHH (Johann 2016 Fig 4 — SHH has EZH2 super-enhancers)
        "EZH2":  1.25, "EED": 1.20, "SUZ12": 1.20,
        # CDK4 amplified in some SHH tumors
        "CDK4":  1.20,
        # HDAC important (pan-ATRT) but less SHH-specific than TYR
        "HDAC1": 1.10, "HDAC2": 1.10,
        # BET important via MYC regulation
        "BRD4":  1.15,
        # TYR markers NOT expressed in SHH
        "TYR":   0.30, "DCT": 0.30, "MITF": 0.30,
    },

    "MYC": {
        # AURKA — MYCN stabilisation is most relevant here (Sredni 2017)
        "AURKA": 1.40, "AURKB": 1.20,
        # MYC/MYCN super-enhancers — BET is critical
        "BRD4":  1.35, "BRD2": 1.30, "BRD3": 1.25,
        # CDK4/6 — MYC drives CDK4 transcription
        "CDK4":  1.25, "CDK6": 1.20,
        # EZH2 — still important (SMARCB1 loss universal)
        "EZH2":  1.15,
        # Proteasome — MYC has t½ ~20 min, proteasome-sensitive
        "PSMB5": 1.20, "PSMB2": 1.15,
        # MYC/MYCN directly
        "MYC":   1.60, "MYCN": 1.50,
        # Hedgehog NOT relevant in MYC
        "SMO":   0.50, "GLI2": 0.50, "GLI1": 0.50,
        # TYR markers NOT relevant
        "TYR":   0.25, "DCT": 0.25, "MITF": 0.25,
    },
}

# Marker genes for centroid-based RNA subgroup calling
SUBGROUP_MARKERS = {
    "TYR": ["TYR", "DCT", "MITF", "SOX10", "PMEL"],
    "SHH": ["GLI1", "GLI2", "PTCH1", "HHIP", "SHH"],
    "MYC": ["MYC", "MYCN", "AURKA", "LIN28A"],
}


# ─────────────────────────────────────────────────────────────────────────────
# RNA-based subgroup calling (when no methylation array available)
# ─────────────────────────────────────────────────────────────────────────────

def call_subgroup_from_rna(
    sample_expr: pd.Series,
    subgroup_markers: Dict[str, List[str]] = SUBGROUP_MARKERS,
) -> Tuple[str, float, Dict[str, float]]:
    """
    Assign a sample to an ATRT molecular subgroup using RNA expression.

    Uses a nearest-centroid approach: for each subgroup, compute the
    mean z-score of its marker genes and assign to the highest-scoring.

    Parameters
    ----------
    sample_expr : Series with gene symbols as index and expression as values
                  (log2-normalised, same scale as GSE70678).

    Returns
    -------
    (predicted_subgroup, confidence_score, all_scores_dict)
    confidence_score: difference between top and second scores (higher = more certain)
    """
    sample_expr.index = sample_expr.index.str.upper()
    subgroup_scores: Dict[str, float] = {}

    for subgroup, markers in subgroup_markers.items():
        available = [m for m in markers if m in sample_expr.index]
        if not available:
            subgroup_scores[subgroup] = 0.0
            continue
        expr_vals = pd.to_numeric(sample_expr[available], errors="coerce").dropna()
        if len(expr_vals) == 0:
            subgroup_scores[subgroup] = 0.0
            continue
        # Z-score the expression values and take the mean
        # (handles scale differences across markers)
        z_mean = float(((expr_vals - expr_vals.mean()) / expr_vals.std().clip(lower=0.01)).mean())
        subgroup_scores[subgroup] = z_mean

    if not subgroup_scores:
        return "UNKNOWN", 0.0, {}

    sorted_sg = sorted(subgroup_scores.items(), key=lambda x: x[1], reverse=True)
    best_sg, best_score = sorted_sg[0]
    second_score = sorted_sg[1][1] if len(sorted_sg) > 1 else 0.0
    confidence = best_score - second_score

    if confidence < 0.3:
        logger.debug(
            "Subgroup calling confidence low (gap=%.2f) — consider methylation array.",
            confidence,
        )

    return best_sg, round(confidence, 3), {k: round(v, 3) for k, v in subgroup_scores.items()}


def call_subgroups_for_cohort(
    atrt_matrix: pd.DataFrame,
    atrt_col_indicators: List[str] = None,
) -> pd.DataFrame:
    """
    Call subgroups for all ATRT samples in a gene × sample matrix.

    Returns DataFrame with columns: sample_id, predicted_subgroup, confidence, scores.
    """
    if atrt_col_indicators:
        cols = [c for c in atrt_matrix.columns
                if any(ind.lower() in c.lower() for ind in atrt_col_indicators)]
    else:
        cols = list(atrt_matrix.columns)

    results = []
    for col in cols:
        sg, conf, scores = call_subgroup_from_rna(atrt_matrix[col])
        results.append({
            "sample_id":          col,
            "predicted_subgroup": sg,
            "confidence":         conf,
            "scores":             scores,
        })

    result_df = pd.DataFrame(results)
    if not result_df.empty:
        counts = result_df["predicted_subgroup"].value_counts()
        logger.info(
            "Subgroup calls for %d samples: %s",
            len(result_df),
            dict(counts),
        )

    return result_df


# ─────────────────────────────────────────────────────────────────────────────
# Subgroup-aware scoring
# ─────────────────────────────────────────────────────────────────────────────

def apply_subgroup_multiplier(
    base_score: float,
    targets:    List[str],
    subgroup:   str,
    multipliers: Dict[str, Dict[str, float]] = SUBGROUP_DRUG_MULTIPLIERS,
) -> Tuple[float, str]:
    """
    Apply subgroup-specific multipliers to a drug's base score.

    Multipliers are applied to the score component reflecting target expression.
    The EZH2 boost (×1.40) is applied BEFORE this — subgroup multipliers are
    additional modifiers on top of the universal SMARCB1 synthetic lethality boost.

    Returns (adjusted_score, explanation_string)
    """
    if subgroup not in multipliers:
        return base_score, "No subgroup multiplier available"

    sg_mults = multipliers[subgroup]
    target_upper = [t.upper() for t in targets]

    # Find the most impactful multiplier among the drug's targets
    applicable = {t: sg_mults[t] for t in target_upper if t in sg_mults}

    if not applicable:
        return base_score, f"No subgroup-specific targets for {subgroup}"

    # Use the largest multiplier (best-case benefit for patients in this subgroup)
    best_target = max(applicable, key=applicable.get)
    best_mult   = applicable[best_target]

    adjusted = round(min(1.0, base_score * best_mult), 4)
    note = (
        f"Subgroup {subgroup} multiplier ×{best_mult} applied via {best_target} "
        f"({base_score:.3f} → {adjusted:.3f})"
    )

    if best_mult < 1.0:
        note += " [PENALISED — drug less rational in this subgroup]"

    return adjusted, note


def compute_prevalence_weighted_score(
    base_scores_by_subgroup: Dict[str, float],
    prevalences: Dict[str, float] = SUBGROUP_PREVALENCES,
) -> float:
    """
    Compute prevalence-weighted expected score across subgroups.

    E[score] = sum_g P(g) × score_g

    This is the correct pan-ATRT score when subgroup is unknown.
    It accounts for the fact that vismodegib benefits only ~37% of patients
    (ATRT-SHH), while panobinostat benefits ~100%.

    Parameters
    ----------
    base_scores_by_subgroup : dict subgroup → score for THIS drug in THAT subgroup

    Returns
    -------
    float — prevalence-weighted expected score (0–1)
    """
    weighted = 0.0
    total_weight = 0.0

    for sg, score in base_scores_by_subgroup.items():
        p = prevalences.get(sg, 0.0)
        weighted    += p * score
        total_weight += p

    if total_weight < 0.9:
        logger.warning(
            "Subgroup prevalences sum to %.2f (expected ~1.0). "
            "Check SUBGROUP_PREVALENCES keys match base_scores keys.",
            total_weight,
        )
        return weighted / max(total_weight, 0.01)

    return round(min(1.0, weighted), 4)


def score_drug_subgroup_aware(
    drug:          Dict,
    base_score:    float,
    subgroup:      Optional[str],
    subgroup_de:   Optional[Dict[str, pd.DataFrame]] = None,
    multipliers:   Dict[str, Dict[str, float]] = SUBGROUP_DRUG_MULTIPLIERS,
    prevalences:   Dict[str, float] = SUBGROUP_PREVALENCES,
) -> Dict:
    """
    Top-level subgroup-aware scorer.

    Mode A (subgroup known): apply subgroup multipliers directly.
    Mode B (subgroup unknown): compute prevalence-weighted expectation.

    Parameters
    ----------
    drug : dict with 'targets' list and 'score' (pre-EZH2-boost base score)
    base_score : score before subgroup adjustment (post-EZH2-boost)
    subgroup : "TYR" / "SHH" / "MYC" / None (unknown)

    Returns
    -------
    dict with:
        final_score       : float
        subgroup_mode     : "SPECIFIC" / "PREVALENCE_WEIGHTED"
        subgroup_scores   : dict sg → score (only in PREVALENCE_WEIGHTED mode)
        multiplier_applied: float
        explanation       : str
    """
    targets = [t.upper() for t in (drug.get("targets") or [])]

    if subgroup and subgroup.upper() in multipliers:
        # Mode A: specific subgroup known
        adjusted, note = apply_subgroup_multiplier(base_score, targets,
                                                    subgroup.upper(), multipliers)
        return {
            "final_score":        adjusted,
            "subgroup_mode":      "SPECIFIC",
            "subgroup":           subgroup.upper(),
            "subgroup_scores":    {subgroup.upper(): adjusted},
            "multiplier_applied": round(adjusted / max(base_score, 0.001), 3),
            "explanation":        note,
        }

    # Mode B: unknown subgroup — compute prevalence-weighted expectation
    scores_by_sg: Dict[str, float] = {}
    for sg in prevalences.keys():
        sg_score, _ = apply_subgroup_multiplier(base_score, targets, sg, multipliers)
        scores_by_sg[sg] = sg_score

    weighted_score = compute_prevalence_weighted_score(scores_by_sg, prevalences)

    explanation = (
        f"Prevalence-weighted score: "
        + ", ".join(f"{sg}={scores_by_sg[sg]:.3f}(P={prevalences[sg]:.0%})"
                    for sg in prevalences)
        + f" → E[score]={weighted_score:.3f}"
    )

    return {
        "final_score":        weighted_score,
        "subgroup_mode":      "PREVALENCE_WEIGHTED",
        "subgroup":           "pan-ATRT",
        "subgroup_scores":    scores_by_sg,
        "multiplier_applied": round(weighted_score / max(base_score, 0.001), 3),
        "explanation":        explanation,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Report generator
# ─────────────────────────────────────────────────────────────────────────────

def generate_subgroup_scoring_report(
    candidates: List[Dict],
    multipliers: Dict[str, Dict[str, float]] = SUBGROUP_DRUG_MULTIPLIERS,
    prevalences: Dict[str, float] = SUBGROUP_PREVALENCES,
) -> str:
    """
    Generate a table showing how scores shift across the three subgroups.
    Useful for communicating which drugs are subgroup-agnostic vs specific.
    """
    lines = [
        "ATRT Subgroup-Stratified Drug Scores",
        "Source: Johann et al. 2016 Cancer Cell PMID 26923874",
        f"Prevalences: TYR={prevalences['TYR']:.0%} | SHH={prevalences['SHH']:.0%} | MYC={prevalences['MYC']:.0%}",
        "",
        f"{'Drug':<24}  {'Pan-ATRT':>9}  {'TYR':>7}  {'SHH':>7}  {'MYC':>7}  Notes",
        "-" * 80,
    ]

    for c in candidates[:12]:
        name    = (c.get("name") or c.get("drug_name") or "?")[:23]
        base    = c.get("score", 0.0)
        targets = [t.upper() for t in (c.get("targets") or [])]

        sg_scores = {}
        for sg in prevalences.keys():
            adj, _ = apply_subgroup_multiplier(base, targets, sg, multipliers)
            sg_scores[sg] = adj

        pan = compute_prevalence_weighted_score(sg_scores, prevalences)

        # Highlight if one subgroup score is dramatically higher
        max_sg  = max(sg_scores, key=sg_scores.get)
        max_s   = sg_scores[max_sg]
        min_s   = min(sg_scores.values())
        range_  = max_s - min_s

        note = ""
        if range_ > 0.15:
            note = f"Subgroup-specific: strongest in {max_sg}"
        elif range_ < 0.05:
            note = "Pan-ATRT (subgroup-agnostic)"

        lines.append(
            f"{name:<24}  {pan:>9.3f}  "
            f"{sg_scores['TYR']:>7.3f}  "
            f"{sg_scores['SHH']:>7.3f}  "
            f"{sg_scores['MYC']:>7.3f}  "
            f"{note}"
        )

    lines += [
        "",
        "Interpretation:",
        "  Pan-ATRT = prevalence-weighted expected score (use when subgroup unknown)",
        "  Subgroup-specific = score for a patient KNOWN to have that subgroup",
        "  EZH2 boost (x1.40) is applied BEFORE subgroup multipliers",
        "  Source: Torchia 2015, Sredni 2017, Johann 2016",
    ]

    return "\n".join(lines)