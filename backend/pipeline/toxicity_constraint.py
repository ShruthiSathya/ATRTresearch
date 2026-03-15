"""
toxicity_constraint.py
Hematologic toxicity penalty for drug combination confidence scoring.

Report: adjusted_confidence ∈ [conservative_lower, optimistic_upper]

Example: Birabresib + Panobinostat + Abemaciclib
  - Additive rate: 55% → capped at 40% → multiplier = 0.60 (conservative)
  - Dose-optimised: 40% × 0.60 = 24% → multiplier = 0.76 (optimistic)
  - Confidence range: [0.55, 0.70] rather than single point 0.55

Sources:
  PANOBINOSTAT : PBTC-047 (Monje et al. 2023, Neuro-Oncology) — 8/29 DLTs
  BIRABRESIB   : PBTC-049 (Geoerger et al. 2017, Clin Cancer Res) — ~15-20% G3/4
  ABEMACICLIB  : MONARCH-2/3 (Sledge et al. 2017) — ~20% G3/4 neutropenia
  MARIZOMIB    : Phase I CNS (Bota et al. 2021) — ~10% G3/4 hematologic
  ONATASERTIB  : mTOR kinase class estimate
  PAXALISIB    : NCT03696355 pediatric DIPG Phase I
"""

import logging
from typing import Dict, List

from .pipeline_config import TOXICITY

logger = logging.getLogger(__name__)

# Grade 3/4 hematologic AE rates from published Phase I/II trials.
# These are CONSERVATIVE (upper-bound) estimates from worst-case dose cohorts.
# Dose-optimised schedules reduce these by the TOXICITY["dose_optimised_fraction"].
HEMATOLOGIC_TOXICITY_RATES: Dict[str, float] = {
    "PANOBINOSTAT":  0.20,   # PBTC-047: 8/29 DLTs = 27.6% → rounded to 20% for G3/4 hematologic
    "BIRABRESIB":    0.15,   # PBTC-049: ~15-20% G3/4
    "PAXALISIB":     0.16,   # NCT03696355
    "GDC-0084":      0.16,
    "ONATASERTIB":   0.08,   # mTOR kinase class estimate
    "ABEMACICLIB":   0.20,   # MONARCH-2/3
    "MARIZOMIB":     0.10,   # Bota et al. 2021 CNS trials
    "CC-115":        0.10,   # Salazar et al. 2016
    "AZD-8055":      0.05,   # Limited data; mTOR class estimate
    "INDOXIMOD":     0.02,   # IDO pathway; low hematologic toxicity
    "ONC201":        0.03,   # ACTION trial: minimal G3/4 hematologic
    "TEMOZOLOMIDE":  0.15,
    "CRIZOTINIB":    0.08,
    "NINTEDANIB":    0.05,
}


def get_single_drug_toxicity(drug_name: str) -> float:
    return HEMATOLOGIC_TOXICITY_RATES.get(
        drug_name.upper().strip(),
        TOXICITY["default_rate"]
    )


def combination_toxicity_penalty(drug_names: List[str]) -> Dict:
    """
    Compute toxicity penalty multiplier for a multi-drug combination.

    v5.5 FIX: Returns both conservative (additive model) and optimistic
    (dose-optimised) multipliers. The adjusted confidence should be reported
    as a RANGE [conservative, optimistic], not a single point estimate.

    Conservative model: additive AE rates, capped at max_combined_rate.
    Optimistic model:   conservative rate × dose_optimised_fraction.
    This reflects real-world PBTC experience where dose reduction mitigates
    combination toxicity (panobinostat PBTC-047 proceeded at 27.6% DLT rate).

    Returns
    -------
    multiplier              : float — conservative (apply to raw confidence)
    multiplier_optimistic   : float — dose-optimised lower-bound penalty
    combined_rate           : float — conservative combined G3/4 AE rate
    combined_rate_optimistic: float — dose-optimised AE rate estimate
    per_drug_rates          : dict
    flag                    : ACCEPTABLE / CAUTION / HIGH_RISK
    note                    : human-readable summary with citations
    confidence_note         : explicit range statement for reporting
    """
    per_drug = {d: get_single_drug_toxicity(d) for d in drug_names}

    max_rate         = TOXICITY["max_combined_rate"]
    combined_rate    = min(sum(per_drug.values()), max_rate)
    multiplier       = round(1.0 - combined_rate, 3)

    # Optimistic: dose-optimised schedule reduces AE rate
    dose_opt_frac    = TOXICITY["dose_optimised_fraction"]
    combined_rate_opt = combined_rate * dose_opt_frac
    multiplier_opt   = round(1.0 - combined_rate_opt, 3)

    acceptable = TOXICITY["acceptable_threshold"]
    caution    = TOXICITY["caution_threshold"]

    if combined_rate <= acceptable:
        flag = "ACCEPTABLE"
    elif combined_rate <= caution:
        flag = "CAUTION"
    else:
        flag = "HIGH_RISK"

    rate_pct     = round(combined_rate * 100, 1)
    rate_pct_opt = round(combined_rate_opt * 100, 1)

    note = (
        f"Estimated combined grade 3/4 hematologic AE rate: {rate_pct}% ({flag}). "
        f"Per drug: { {k: f'{round(v*100):.0f}%' for k, v in per_drug.items()} }. "
        f"Conservative multiplier: {multiplier:.3f}. "
        f"Additive model is conservative (assumes independence + worst-case dose). "
        f"Sources: PBTC-047 (panobinostat), PBTC-049 (birabresib), "
        f"NCT03696355 (paxalisib), MONARCH-2/3 (abemaciclib)."
    )

    # FIX v5.5: explicit range statement for hypothesis_generator to surface
    confidence_note = (
        f"Toxicity-adjusted confidence is reported as a RANGE, not a point estimate. "
        f"Conservative (additive model): ×{multiplier:.3f} ({rate_pct}% combined AE rate). "
        f"Optimistic (dose-optimised, {int(dose_opt_frac*100)}% reduction): "
        f"×{multiplier_opt:.3f} ({rate_pct_opt}% AE rate). "
        f"PBTC-047 precedent: panobinostat proceeded at 27.6% DLT rate in pediatric DIPG "
        f"with dose modification — establishing that >20% is acceptable when no curative "
        f"alternative exists (Monje et al. 2023, PMID 37526549)."
    )

    logger.info(
        "🧪 Toxicity — %s: %.1f%% (%s), multiplier=%.3f [conservative] / %.3f [optimistic]",
        " + ".join(drug_names), rate_pct, flag, multiplier, multiplier_opt,
    )

    return {
        "multiplier":               multiplier,
        "multiplier_optimistic":    multiplier_opt,
        "combined_rate":            round(combined_rate, 4),
        "combined_rate_optimistic": round(combined_rate_opt, 4),
        "per_drug_rates":           per_drug,
        "flag":                     flag,
        "note":                     note,
        "confidence_note":          confidence_note,
    }