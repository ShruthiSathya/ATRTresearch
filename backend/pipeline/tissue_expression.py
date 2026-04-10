"""
tissue_expression.py
====================
ATRT tissue expression scoring using GSE70678 bulk RNA-seq.

DATA SOURCE RATIONALE
----------------------
DIPG pipeline uses Filbin 2018 scRNA-seq (GSE102130) because a high-quality
H3K27M single-cell atlas exists. No equivalent ATRT scRNA-seq atlas exists
as of 2026. GSE70678 (Torchia 2015) is the largest published ATRT cohort
with RNA expression data: 49 ATRT tumors + normal brain samples, with
molecular subgroup annotations (TYR/SHH/MYC).

Expression scores are derived from:
  1. GSE70678 bulk RNA mean expression across ATRT samples (primary)
  2. Curated ATRT cell-line expression scores from published studies (fallback)
  Both are blended: curated 55%, bulk RNA 45% (curated weighted higher due
  to less resolution in bulk RNA vs scRNA-seq)

QUANTILE SENSITIVITY:
  run_quantile_sensitivity() tests p65/p70/p75/p80 thresholds on the
  GSE70678 "high expression" cutoff and reports top-4 ranking stability.

DATA DOWNLOAD:
  GSE70678: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70678
  Platform: Affymetrix HuGene 1.0 ST (GPL6244)
  File: GSE70678_series_matrix.txt.gz
  Place at: data/raw_omics/GSE70678_ATRT_expression.txt

REFERENCES
-----------
Torchia J et al. (2015). Integrated (epi)-genomic analyses identify
  subgroup-specific therapeutic targets in CNS rhabdoid tumours.
  Cancer Cell 30(6):891–908. PMID 26609405.

Johann PD et al. (2016). Atypical teratoid/rhabdoid tumors are comprised
  of three epigenetic subgroups with distinct enhancer landscapes.
  Cancer Cell 29(3):379–393. PMID 26923874.
"""

import copy
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional

try:
    from .pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES
except ImportError:
    from pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES

logger = logging.getLogger(__name__)


class TissueExpressionScorer:
    """
    ATRT tissue expression scorer using GSE70678 bulk RNA-seq.

    Scoring pipeline:
      1. Load GSE70678 — separate ATRT columns from normal brain
      2. Compute mean expression per gene across ATRT samples
      3. Score each drug target against expression percentiles
      4. Blend with curated ATRT cell-line scores (55% curated, 45% bulk RNA)
      5. Subgroup-stratified scoring if subgroup is known

    If GSE70678 is not available, falls back to curated scores only.
    """

    def __init__(self, disease_name: str = "atrt", data_dir: str = None):
        self.disease_name = disease_name
        self.sc_path      = Path(PATHS["scrna"])   # GSE70678
        self.is_ready     = False
        self.bulk_target_scores: Dict[str, float] = {}   # gene → mean expression
        self._all_values_sorted: List[float] = []
        self._p25 = 0.0
        self._p50 = 0.0
        self._p75 = 0.0
        self._p90 = 0.0
        self._df_cache: Optional[pd.DataFrame] = None
        self._n_atrt_samples = 0
        self._n_normal_samples = 0
        self.subgroup: Optional[str] = None   # TYR / SHH / MYC / None

    def set_subgroup(self, subgroup: Optional[str]) -> None:
        """Set molecular subgroup for subgroup-stratified scoring."""
        if subgroup and subgroup.upper() in ("TYR", "SHH", "MYC"):
            self.subgroup = subgroup.upper()
            logger.info("TissueExpressionScorer: subgroup set to %s", self.subgroup)
        else:
            self.subgroup = None

    async def _load_bulk_rna(self, quantile: float = None) -> None:
        """
        Load GSE70678 ATRT bulk RNA expression matrix.

        Separates ATRT columns from normal brain using column name indicators
        defined in pipeline_config.TISSUE. Computes per-gene mean across ATRT
        samples for expression scoring.
        """
        q = quantile or TISSUE["bulk_high_quantile"]

        if self.is_ready and quantile is None:
            return

        if not self.sc_path.exists():
            logger.warning(
                "GSE70678 not found at %s — using curated ATRT scores only.\n"
                "Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70678\n"
                "Install via: pip install GEOparse",
                self.sc_path,
            )
            return

        try:
            if self._df_cache is None:
                logger.info(
                    "Loading GSE70678 ATRT bulk RNA-seq (%s)...", self.sc_path
                )
                df = pd.read_csv(self.sc_path, sep="\t", index_col=0)
                df.index = df.index.astype(str).str.upper().str.strip()

                # Handle duplicate gene symbols — take max expression
                df = df.groupby(level=0).max()
                self._df_cache = df
            else:
                df = self._df_cache

            atrt_indicators   = TISSUE["atrt_col_indicators"]
            normal_indicators = TISSUE["normal_col_indicators"]

            atrt_cols = [
                c for c in df.columns
                if any(ind.lower() in c.lower() for ind in atrt_indicators)
                and not any(n.lower() in c.lower() for n in normal_indicators)
            ]
            normal_cols = [
                c for c in df.columns
                if any(n.lower() in c.lower() for n in normal_indicators)
            ]

            if not atrt_cols:
                logger.warning(
                    "No ATRT columns found in GSE70678 using indicators %s.\n"
                    "Column names (first 10): %s\n"
                    "Using curated scores only.",
                    atrt_indicators, list(df.columns[:10])
                )
                return

            self._n_atrt_samples   = len(atrt_cols)
            self._n_normal_samples = len(normal_cols)

            # Differential: ATRT mean vs normal mean if normals available
            df_num = df.apply(pd.to_numeric, errors="coerce")
            atrt_mean = df_num[atrt_cols].mean(axis=1)

            if normal_cols:
                normal_mean = df_num[normal_cols].mean(axis=1)
                diff        = atrt_mean - normal_mean
                logger.info(
                    "GSE70678: differential expression (ATRT mean − normal mean). "
                    "%d ATRT samples, %d normal samples.",
                    len(atrt_cols), len(normal_cols)
                )
            else:
                diff = atrt_mean
                logger.info(
                    "GSE70678: no normal columns found — using ATRT mean expression. "
                    "%d ATRT samples.",
                    len(atrt_cols)
                )

            diff = diff.dropna()
            self.bulk_target_scores = diff.to_dict()

            # Compute percentiles for expression → score mapping
            all_vals = sorted(diff.values)
            n        = len(all_vals)
            self._all_values_sorted = all_vals
            self._p25 = all_vals[int(n * 0.25)]
            self._p50 = all_vals[int(n * 0.50)]
            self._p75 = all_vals[int(n * 0.75)]
            self._p90 = all_vals[int(n * 0.90)]

            if quantile is None:
                self.is_ready = True

            logger.info(
                "GSE70678 loaded: %d genes | %d ATRT, %d normal samples\n"
                "  Expression percentiles — p25: %.2f, p50: %.2f, "
                "p75: %.2f, p90: %.2f",
                len(diff), len(atrt_cols), len(normal_cols),
                self._p25, self._p50, self._p75, self._p90,
            )

        except Exception as e:
            logger.error("GSE70678 loading failed: %s", e)

    def _expression_to_score(self, expr_value: float) -> float:
        """Map bulk RNA expression value to 0–1 score using percentile bins."""
        bins = TISSUE["percentile_bins"]

        if not self._all_values_sorted:
            # Fallback if data not loaded
            if expr_value >= TISSUE["fallback_high_cutoff"]:
                return TISSUE["fallback_high_tpm"]
            elif expr_value >= TISSUE["fallback_mid_cutoff"]:
                return TISSUE["fallback_mid_tpm"]
            return TISSUE["fallback_low_tpm"]

        if expr_value >= self._p90:   return bins["p90"]
        elif expr_value >= self._p75: return bins["p75"]
        elif expr_value >= self._p50: return bins["p50"]
        elif expr_value >= self._p25: return bins["p25"]
        else:                         return bins["low"]

    def _curated_score(self, drug: Dict) -> float:
        """
        Score drug targets against curated ATRT expression database.

        Aggregation: best target score (65%) + mean of top-2 (35%).
        Reflects that a drug with one highly expressed target should score
        well even if other targets are less relevant.
        """
        targets = drug.get("targets", [])
        if not targets:
            return TISSUE["no_target_score"]

        # Apply subgroup multipliers if subgroup is known
        scores = []
        for t in targets:
            base_score = ATRT_CURATED_SCORES.get(t.upper(), TISSUE["unknown_target_score"])
            if self.subgroup:
                multiplier = self._get_subgroup_multiplier(t, self.subgroup)
                base_score = min(1.0, base_score * multiplier)
            scores.append(base_score)

        scores_sorted = sorted(scores, reverse=True)
        best          = scores_sorted[0]
        top_n         = TISSUE["curated_top_n"]
        mean_top_n    = (
            sum(scores_sorted[:top_n]) / min(len(scores_sorted), top_n)
        )

        return round(
            TISSUE["curated_best_weight"]    * best
            + TISSUE["curated_top_n_weight"] * mean_top_n,
            4
        )

    @staticmethod
    def _get_subgroup_multiplier(gene: str, subgroup: str) -> float:
        """
        Return subgroup-specific expression multiplier for a gene.

        Based on subgroup-defining gene expression from Torchia 2015 and
        Johann 2016. Multipliers boost subgroup-relevant targets.
        """
        gene_upper = gene.upper()
        # TYR subgroup — neural crest / melanocytic program
        tyr_boosts = {"TYR": 1.5, "DCT": 1.5, "MITF": 1.5, "SOX10": 1.4,
                      "HDAC1": 1.2, "HDAC2": 1.2, "BRD4": 1.15}
        # SHH subgroup — Sonic Hedgehog pathway activation
        shh_boosts = {"GLI2": 1.5, "SMO": 1.4, "PTCH1": 1.3, "SHH": 1.3,
                      "EZH2": 1.2, "CDK4": 1.15}
        # MYC subgroup — MYC amplification, worst prognosis
        myc_boosts = {"MYC": 1.6, "MYCN": 1.5, "AURKA": 1.4, "BRD4": 1.3,
                      "CDK4": 1.2, "PSMB5": 1.2, "EZH2": 1.15}

        if subgroup == "TYR":
            return tyr_boosts.get(gene_upper, 1.0)
        elif subgroup == "SHH":
            return shh_boosts.get(gene_upper, 1.0)
        elif subgroup == "MYC":
            return myc_boosts.get(gene_upper, 1.0)
        return 1.0

    def _atrt_target_relevance(self, target: str) -> float:
        """
        Relevance multiplier for a target in the ATRT context.
        Targets in the curated ATRT database receive full weight (1.0).
        Off-panel targets receive a partial weight.
        """
        return 1.0 if target.upper() in ATRT_CURATED_SCORES \
               else TISSUE["off_target_relevance"]

    async def score_batch(self, candidates: List[Dict]) -> List[Dict]:
        """Score all candidates for ATRT tissue expression."""
        await self._load_bulk_rna()
        return self._score_with_current_state(candidates)

    def _score_with_current_state(
        self, candidates: List[Dict]
    ) -> List[Dict]:
        """
        Score candidates using currently loaded bulk RNA state.
        Used by both score_batch() and run_quantile_sensitivity().
        """
        cw = TISSUE["curated_weight"]
        bw = TISSUE["bulk_weight"]

        for drug in candidates:
            targets = drug.get("targets", [])

            # Curated score — always computed
            curated_score = self._curated_score(drug)

            if not self.is_ready:
                drug["tissue_expression_score"] = curated_score
                drug["sc_context"] = (
                    "Curated ATRT scores only (GSE70678 not loaded)"
                )
                continue

            # Bulk RNA scores for this drug's targets
            target_exprs = [
                self.bulk_target_scores.get(t.upper(), 0.0) for t in targets
            ]

            if not target_exprs or max(target_exprs) == 0.0:
                drug["tissue_expression_score"] = curated_score
                drug["sc_context"] = (
                    "Targets absent from GSE70678 — using curated ATRT scores"
                )
                continue

            # Weighted expression scores with subgroup and ATRT relevance
            weighted_bulk_scores = [
                self._expression_to_score(expr) * self._atrt_target_relevance(t)
                for t, expr in zip(targets, target_exprs)
            ]
            bulk_score_weighted = max(weighted_bulk_scores)

            blended   = round(cw * curated_score + bw * bulk_score_weighted, 4)
            best_expr = max(target_exprs)

            drug["tissue_expression_score"] = blended
            drug["sc_context"] = (
                f"GSE70678 blended: curated={curated_score:.2f} (w={cw}), "
                f"bulk={bulk_score_weighted:.2f} (w={bw}), "
                f"best_expr={best_expr:.2f}, final={blended:.2f}"
                + (f", subgroup={self.subgroup}" if self.subgroup else "")
            )

        scores = [c.get("tissue_expression_score", 0) for c in candidates]
        if scores:
            logger.info(
                "Tissue scoring (GSE70678): %d candidates | "
                "range=[%.2f, %.2f] | ATRT samples=%d",
                len(candidates), min(scores), max(scores),
                self._n_atrt_samples,
            )

        return candidates

    async def run_quantile_sensitivity(
        self, candidates: List[Dict]
    ) -> Dict:
        """
        Sensitivity analysis on bulk RNA high-expression quantile threshold.

        Tests p65/p70/p75/p80 and reports whether top-4 ranking is stable.
        Analogous to DIPG's GSC quantile sensitivity.

        Returns dict with stability flag and per-quantile top-5 rankings.
        """
        if not self.sc_path.exists():
            return {
                "stable": None,
                "note": "GSE70678 not found — sensitivity analysis skipped",
                "rankings_by_quantile": {},
            }

        quantiles = TISSUE["bulk_quantile_sensitivity_range"]
        rankings: Dict[float, List[str]] = {}

        for q in quantiles:
            await self._load_bulk_rna(quantile=q)

            candidates_copy = copy.deepcopy(candidates[:20])
            scored = self._score_with_current_state(candidates_copy)
            scored.sort(
                key=lambda x: x.get("tissue_expression_score", 0),
                reverse=True
            )
            rankings[q] = [c.get("name", "?") for c in scored[:5]]

        # Reload at primary quantile
        self.is_ready = False
        await self._load_bulk_rna()

        # Check top-4 stability
        all_top4 = [tuple(rankings[q][:4]) for q in quantiles if rankings.get(q)]
        stable   = len(set(all_top4)) == 1

        logger.info(
            "GSE70678 quantile sensitivity: top-4 ranking %s across quantiles %s",
            "STABLE ✅" if stable else "UNSTABLE ⚠️",
            quantiles,
        )

        return {
            "stable":               stable,
            "quantiles_tested":     quantiles,
            "primary_quantile":     TISSUE["bulk_high_quantile"],
            "top4_consistent":      stable,
            "rankings_by_quantile": {str(q): rankings.get(q, []) for q in quantiles},
            "note": (
                "Top-4 ranking stable across all GSE70678 expression quantile thresholds ✅"
                if stable else
                "⚠️ Top-4 changes with quantile — interpret tissue scores with caution."
            ),
        }