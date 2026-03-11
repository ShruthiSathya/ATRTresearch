"""
toxicity_constraint.py — v5.4
Hematologic toxicity penalty for drug combination confidence scoring.

All rates are grade 3/4 hematologic AE frequencies from published trials.
Thresholds and caps are read from pipeline_config.TOXICITY — no magic numbers here.

Sources:
  PANOBINOSTAT : PBTC-047 (Monje et al. 2023, Neuro-Oncology)
  BIRABRESIB   : PBTC-049 (Geoerger et al. 2017, Clin Cancer Res)
  PAXALISIB    : NCT03696355 pediatric DIPG Phase I (PMC7650438)
                 Primary toxicities are rash/hyperglycemia; rate here = neutropenia.
  ONATASERTIB  : mTOR kinase class estimate; no DIPG-specific data published.
  ABEMACICLIB  : MONARCH-2/3 (Sledge et al. 2017/2019)
  MARIZOMIB    : Phase I CNS (Bota et al. 2021, Neuro-Oncology)
  CC-115       : Phase I (Salazar et al. 2016, Lancet Oncol)
"""

import logging
from typing import Dict, List

from .pipeline_config import TOXICITY

logger = logging.getLogger(__name__)

# Grade 3/4 hematologic AE rates from published Phase I/II trials
HEMATOLOGIC_TOXICITY_RATES: Dict[str, float] = {
    "PANOBINOSTAT":  0.20,
    "BIRABRESIB":    0.15,
    "PAXALISIB":     0.16,
    "GDC-0084":      0.16,
    "ONATASERTIB":   0.08,
    "ABEMACICLIB":   0.20,
    "MARIZOMIB":     0.10,
    "CC-115":        0.10,
    "AZD-8055":      0.05,
    "INDOXIMOD":     0.02,
    "ONC201":        0.03,
    "TEMOZOLOMIDE":  0.15,
    "CRIZOTINIB":    0.08,
    "NINTEDANIB":    0.05,
}


def get_single_drug_toxicity(drug_name: str) -> float:
    """Return grade 3/4 hematologic AE rate for a single drug."""
    return HEMATOLOGIC_TOXICITY_RATES.get(
        drug_name.upper().strip(),
        TOXICITY["default_rate"]
    )


def combination_toxicity_penalty(drug_names: List[str]) -> Dict:
    """
    Compute toxicity penalty multiplier for a multi-drug combination.

    Model: additive hematologic AE rates, capped at TOXICITY["max_combined_rate"].
    Multiplier = 1 - combined_rate.

    Apply to raw confidence: adjusted_confidence = raw_confidence * multiplier.

    Returns:
        multiplier      : float — apply to raw confidence
        combined_rate   : float — estimated combined grade 3/4 hematologic rate
        per_drug_rates  : dict  — individual rates per drug
        flag            : str   — ACCEPTABLE / CAUTION / HIGH_RISK
        note            : str   — human-readable summary with citations
    """
    per_drug = {d: get_single_drug_toxicity(d) for d in drug_names}

    max_rate      = TOXICITY["max_combined_rate"]
    combined_rate = min(sum(per_drug.values()), max_rate)
    multiplier    = round(1.0 - combined_rate, 3)

    acceptable = TOXICITY["acceptable_threshold"]
    caution    = TOXICITY["caution_threshold"]

    if combined_rate <= acceptable:
        flag = "ACCEPTABLE"
    elif combined_rate <= caution:
        flag = "CAUTION"
    else:
        flag = "HIGH_RISK"

    rate_pct = round(combined_rate * 100, 1)
    note = (
        f"Estimated combined grade 3/4 hematologic AE rate: {rate_pct}% ({flag}). "
        f"Per drug: { {k: f'{round(v*100):.0f}%' for k, v in per_drug.items()} }. "
        f"Multiplier: {multiplier:.3f}. "
        f"Additive model is conservative. "
        f"Sources: PBTC-047 (panobinostat), PBTC-049 (birabresib), "
        f"NCT03696355 (paxalisib)."
    )

    logger.info(
        "🧪 Toxicity — %s: %.1f%% (%s), multiplier=%.3f",
        " + ".join(drug_names), rate_pct, flag, multiplier,
    )

    return {
        "multiplier":     multiplier,
        "combined_rate":  round(combined_rate, 4),
        "per_drug_rates": per_drug,
        "flag":           flag,
        "note":           note,
    }