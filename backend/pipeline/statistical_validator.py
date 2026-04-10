"""
statistical_validator.py — ATRT Statistical Validation
=======================================================
ATRT differs fundamentally from DIPG in what to test statistically.

DIPG had a Fisher's exact test for H3K27M + CDKN2A-del co-occurrence
(mutual exclusivity), because those are two competing oncogenic mechanisms.

ATRT has NO co-occurrence hypothesis to test. SMARCB1 biallelic loss is
the defining event in ~95% of cases. There is no second alteration to test
co-occurrence against.

Instead, this module:
  1. Validates SMARCB1 loss prevalence against published priors
  2. Tests subgroup distribution vs published priors (Johann 2016)
  3. Reports confidence intervals on cell line coverage

References:
  Hasselblatt 2011 Acta Neuropathologica — SMARCB1 prevalence ~95%
  Johann 2016 Cancer Cell — subgroup prevalence TYR/SHH/MYC
  Frühwald 2020 CNS Oncology — ATRT biology overview
"""

import logging
import math
from typing import Dict, Optional, Tuple

logger = logging.getLogger(__name__)

# Published SMARCB1 prevalence in ATRT
# Source: Hasselblatt M et al. 2011 Acta Neuropathol 122(4):417-424. PMID 20625942
PUBLISHED_SMARCB1_PREVALENCE = 0.95

# Published subgroup prevalences — Johann et al. 2016 Cancer Cell 29(3):379-393
# n=150 ATRT samples with methylation-based subgroup calls
PUBLISHED_SUBGROUP_PREVALENCES = {
    "TYR": 0.36,
    "SHH": 0.37,
    "MYC": 0.27,
}

MIN_SAMPLES_FOR_PROPORTION_CI = 5


class StatisticalValidator:
    """
    ATRT statistical validation module.

    Replaces the DIPG co-occurrence Fisher's exact test with ATRT-appropriate
    statistics: SMARCB1 prevalence validation and subgroup distribution testing.
    """

    def validate_input(self, genomic_stats: dict) -> Tuple[bool, str]:
        if not genomic_stats:
            return False, "genomic_stats dict is empty — genomic data not loaded"

        smarcb1 = genomic_stats.get("smarcb1_loss_count", 0)
        total   = genomic_stats.get("total_samples", 0)

        if total == 0:
            return False, (
                "total_samples is 0. CBTN ATRT files not found. "
                "Pipeline will use GSE70678 fallback or prevalence priors. "
                "This is non-fatal — all scoring continues."
            )

        return True, "ok"

    def calculate_cooccurrence_p_value(
        self,
        genomic_stats: dict,
        total_samples: int = 0,
    ) -> Optional[float]:
        """
        ATRT override: no co-occurrence p-value is computed.

        In ATRT, SMARCB1 loss IS the defining oncogenic event. There is no
        second mutation to test for co-occurrence or mutual exclusivity.

        Returns None to signal "not applicable" to hypothesis_generator,
        which will use the "N/A — SMARCB1 loss is defining event" label.
        """
        return None

    def validate_smarcb1_prevalence(self, genomic_stats: dict) -> Dict:
        """
        Validate observed SMARCB1 loss fraction against published 95% prior.

        Uses Wilson confidence interval for a proportion.
        Source: Wilson EB 1927 J Am Stat Assoc — binomial proportion CI.
        """
        is_valid, reason = self.validate_input(genomic_stats)
        if not is_valid:
            return {
                "validated":    False,
                "reason":       reason,
                "observed_prevalence": None,
                "published_prevalence": PUBLISHED_SMARCB1_PREVALENCE,
                "concordant":   None,
            }

        n        = genomic_stats.get("total_samples", 0)
        smarcb1  = genomic_stats.get("smarcb1_loss_count", 0)
        observed = smarcb1 / max(n, 1)

        # Wilson 95% CI
        z     = 1.96
        denom = 1 + z**2 / n
        centre = (observed + z**2 / (2*n)) / denom
        margin = z * math.sqrt(observed * (1 - observed) / n + z**2 / (4*n**2)) / denom
        ci_lower = max(0.0, centre - margin)
        ci_upper = min(1.0, centre + margin)

        concordant = ci_lower <= PUBLISHED_SMARCB1_PREVALENCE <= ci_upper

        logger.info(
            "SMARCB1 prevalence: observed=%.1f%% [%.1f%%-%.1f%%] vs published=%.0f%% — %s",
            observed * 100, ci_lower * 100, ci_upper * 100,
            PUBLISHED_SMARCB1_PREVALENCE * 100,
            "concordant ✅" if concordant else "discordant ⚠️",
        )

        return {
            "validated":              True,
            "n_samples":              n,
            "smarcb1_loss_count":     smarcb1,
            "observed_prevalence":    round(observed, 4),
            "published_prevalence":   PUBLISHED_SMARCB1_PREVALENCE,
            "ci_lower_95":            round(ci_lower, 4),
            "ci_upper_95":            round(ci_upper, 4),
            "concordant":             concordant,
            "note": (
                f"SMARCB1 biallelic loss: {smarcb1}/{n} ({observed:.0%}) vs "
                f"published {PUBLISHED_SMARCB1_PREVALENCE:.0%} "
                f"(Hasselblatt 2011 Acta Neuropathol). "
                f"95% CI [{ci_lower:.0%}, {ci_upper:.0%}]."
            ),
        }

    def validate_subgroup_distribution(self, subgroup_counts: Dict) -> Dict:
        """
        Test whether observed subgroup distribution matches published priors
        using a chi-squared goodness-of-fit test.

        Published priors: Johann et al. 2016 Cancer Cell, n=150.
        TYR=36%, SHH=37%, MYC=27%.
        """
        observed = {k: subgroup_counts.get(k, 0) for k in PUBLISHED_SUBGROUP_PREVALENCES}
        n_total  = sum(observed.values())

        if n_total < MIN_SAMPLES_FOR_PROPORTION_CI:
            return {
                "tested":  False,
                "reason":  f"n={n_total} too small for chi-squared test (min {MIN_SAMPLES_FOR_PROPORTION_CI})",
                "observed": observed,
                "expected": PUBLISHED_SUBGROUP_PREVALENCES,
            }

        expected_counts = {k: PUBLISHED_SUBGROUP_PREVALENCES[k] * n_total
                          for k in PUBLISHED_SUBGROUP_PREVALENCES}

        # Chi-squared statistic
        chi2 = sum(
            (observed[k] - expected_counts[k])**2 / max(expected_counts[k], 0.01)
            for k in observed
        )
        df = len(observed) - 1

        # Approximate p-value using chi-squared CDF
        # scipy.stats.chi2.sf(chi2, df) — implemented inline to avoid dependency
        p_approx = self._chi2_sf(chi2, df)

        return {
            "tested":            True,
            "n_total":           n_total,
            "chi2_statistic":    round(chi2, 4),
            "degrees_freedom":   df,
            "p_value_approx":    round(p_approx, 4),
            "significant":       p_approx < 0.05,
            "observed":          {k: v / max(n_total, 1) for k, v in observed.items()},
            "expected":          PUBLISHED_SUBGROUP_PREVALENCES,
            "note": (
                f"Subgroup distribution vs Johann 2016 (n={n_total}): "
                f"chi2={chi2:.2f}, df={df}, p={p_approx:.4f}. "
                f"{'Significant deviation from published priors ⚠️' if p_approx < 0.05 else 'Consistent with published priors ✅'}"
            ),
        }

    def format_p_value_for_report(self, p_value: Optional[float]) -> str:
        """Return ATRT-appropriate statistical note (no co-occurrence framing)."""
        if p_value is None:
            return (
                "N/A — SMARCB1 biallelic loss is the defining event in ATRT. "
                "No co-occurrence hypothesis applies. "
                "(Hasselblatt 2011 Acta Neuropathol PMID 20625942)"
            )
        if math.isnan(p_value):
            return "N/A — insufficient sample counts"
        if p_value < 0.001:
            return f"{p_value:.2e} ✅ (significant)"
        if p_value < 0.05:
            return f"{p_value:.4f} ✅ (significant)"
        return f"{p_value:.4f} (not significant)"

    def priority_from_p_value(self, p_value: Optional[float]) -> str:
        """Derive hypothesis priority. ATRT returns COMPUTATIONAL by default."""
        if p_value is None or (isinstance(p_value, float) and math.isnan(p_value)):
            return "COMPUTATIONAL"
        if p_value < 0.05:
            return "HIGH"
        return "MODERATE"

    @staticmethod
    def _chi2_sf(x: float, df: int) -> float:
        """
        Survival function for chi-squared distribution (approximate).
        Uses regularized incomplete gamma function approximation.
        Adequate for df=1,2 (our use case).
        """
        try:
            from scipy.stats import chi2
            return float(chi2.sf(x, df))
        except ImportError:
            # Fallback: simple approximation for df=1,2
            if df == 1:
                return 1.0 - math.erf(math.sqrt(x / 2))
            elif df == 2:
                return math.exp(-x / 2)
            else:
                return 0.5  # unknown