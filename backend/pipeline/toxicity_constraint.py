"""
toxicity_constraint.py
======================
Hematologic toxicity penalty for drug combination confidence scoring.

CHANGES FROM v1.0
-----------------
FIX: Added try/except fallback for import of pipeline_config.TOXICITY.
     The original `from .pipeline_config import TOXICITY` raises ImportError
     when this module is imported outside of the package context (e.g. during
     standalone testing or direct execution). All other pipeline modules already
     use this pattern — toxicity_constraint.py was the only module missing it.

Report: adjusted_confidence ∈ [conservative_lower, optimistic_upper]

Example: Tazemetostat + Panobinostat + Alisertib
  - Additive rate: 50% → capped at 40% → multiplier = 0.60 (conservative)
  - Dose-optimised: 40% × 0.60 = 24% → multiplier = 0.76 (optimistic)
  - Confidence range: [conservative, optimistic]

Sources:
  TAZEMETOSTAT : Gounder 2020 JCO PMID 33166238 — 4.8% G3/4 hematologic
  PANOBINOSTAT : PBTC-047 (Monje 2023 Nat Med PMID 37526549) — 8/29 DLTs
  ALISERTIB    : Geller 2015 Cancer PMID 25921089 — ~25% G3/4 (pediatric CNS)
  BIRABRESIB   : PBTC-049 (Geoerger 2017 Clin Cancer Res) — ~15-20% G3/4
  ABEMACICLIB  : MONARCH-2/3 (Sledge 2017 JCO PMID 28580882) — ~20% G3/4
  MARIZOMIB    : Bota 2021 Neuro-Oncology PMID 33300566 — ~10% G3/4
"""

import logging
from typing import Dict, List

# FIX: robust relative/absolute import pattern (was bare `from .pipeline_config import`)
try:
    from .pipeline_config import TOXICITY
except ImportError:
    try:
        from pipeline_config import TOXICITY
    except ImportError:
        # Minimal inline fallback so the module remains functional even when
        # imported outside the package (e.g. in unit tests)
        TOXICITY = {
            "drug_ae_rates": {
                "TAZEMETOSTAT":  0.05,
                "ALISERTIB":     0.25,
                "PANOBINOSTAT":  0.20,
                "BIRABRESIB":    0.15,
                "ABEMACICLIB":   0.20,
                "MARIZOMIB":     0.10,
                "VISMODEGIB":    0.03,
                "SONIDEGIB":     0.04,
                "ONC201":        0.03,
                "PAXALISIB":     0.16,
                "INDOXIMOD":     0.02,
                "ONATASERTIB":   0.08,
                "TEMOZOLOMIDE":  0.15,
            },
            "default_rate":              0.10,
            "max_combined_rate":         0.40,
            "acceptable_threshold":      0.20,
            "caution_threshold":         0.30,
            "dose_optimised_fraction":   0.60,
        }

logger = logging.getLogger(__name__)

# Grade 3/4 hematologic AE rates — these are loaded from TOXICITY config above
# (published Phase I/II trials; conservative upper-bound estimates)
HEMATOLOGIC_TOXICITY_RATES: Dict[str, float] = TOXICITY["drug_ae_rates"]


def get_single_drug_toxicity(drug_name: str) -> float:
    """Return G3/4 hematologic AE rate for a drug. Uses default if unknown."""
    return HEMATOLOGIC_TOXICITY_RATES.get(
        drug_name.upper().strip(),
        TOXICITY.get("default_rate", 0.10),
    )


def combination_toxicity_penalty(drug_names: List[str]) -> Dict:
    """
    Compute toxicity penalty multiplier for a multi-drug combination.

    Returns both conservative (additive model) and optimistic (dose-optimised)
    multipliers. Adjusted confidence should be reported as a RANGE
    [conservative, optimistic], not a single point estimate.

    Conservative model: additive AE rates, capped at max_combined_rate.
    Optimistic model:   conservative rate × dose_optimised_fraction.

    The PBTC-047 precedent (panobinostat 27.6% DLT with dose modification
    in pediatric DIPG) establishes dose optimisation is clinically feasible.
    Source: Monje et al. 2023 Nat Med PMID 37526549.

    Parameters
    ----------
    drug_names : list of drug name strings (case insensitive)

    Returns
    -------
    dict with:
        multiplier               : float — conservative penalty
        multiplier_optimistic    : float — dose-optimised penalty
        combined_rate            : float — conservative combined AE rate
        combined_rate_optimistic : float — dose-optimised AE rate
        per_drug_rates           : dict  — individual drug rates
        flag                     : ACCEPTABLE / CAUTION / HIGH_RISK
        note                     : citation-backed summary
        confidence_note          : explicit range statement for reporting
    """
    per_drug = {d: get_single_drug_toxicity(d) for d in drug_names}

    max_rate      = TOXICITY.get("max_combined_rate", 0.40)
    combined_rate = min(sum(per_drug.values()), max_rate)
    multiplier    = round(1.0 - combined_rate, 3)

    dose_opt_frac     = TOXICITY.get("dose_optimised_fraction", 0.60)
    combined_rate_opt = combined_rate * dose_opt_frac
    multiplier_opt    = round(1.0 - combined_rate_opt, 3)

    acceptable = TOXICITY.get("acceptable_threshold", 0.20)
    caution    = TOXICITY.get("caution_threshold", 0.30)

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
        f"Sources: PBTC-047 (panobinostat PMID 37526549), PBTC-049 (birabresib PMID 28108534), "
        f"Geller 2015 (alisertib PMID 25921089), Gounder 2020 (tazemetostat PMID 33166238)."
    )

    confidence_note = (
        f"Toxicity-adjusted confidence is reported as a RANGE, not a point estimate. "
        f"Conservative (additive model): ×{multiplier:.3f} ({rate_pct}% combined AE rate). "
        f"Optimistic (dose-optimised, {int(dose_opt_frac*100)}% reduction): "
        f"×{multiplier_opt:.3f} ({rate_pct_opt}% AE rate). "
        f"PBTC-047 precedent: panobinostat proceeded at 27.6% DLT rate in pediatric DIPG "
        f"with dose modification — establishing that >20% is acceptable when no curative "
        f"alternative exists (Monje et al. 2023 PMID 37526549)."
    )

    logger.info(
        "Toxicity — %s: %.1f%% (%s), multiplier=%.3f [conservative] / %.3f [optimistic]",
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