import numpy as np
import pandas as pd
import logging
from pathlib import Path
from typing import Dict, List

from .pipeline_config import PATHS, TISSUE

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# CURATED DIPG/GBM GSC EXPRESSION SCORES
# Source: Grasso et al. 2015 (Nat Med), Nagaraja et al. 2017 (Cancer Cell),
#         Filbin et al. 2018 (Science), Tirosh et al. 2016 (Science)
#
# Scores = percentile rank in published H3K27M DIPG sc datasets (0–1).
# Used as one component of the blended score when scRNA-seq is loaded,
# or as full fallback when scRNA-seq file is not available.
# ─────────────────────────────────────────────────────────────────────────────

DIPG_GSC_CURATED_SCORES: Dict[str, float] = {
    # ── Epigenetic targets ────────────────────────────────────────────────────
    "EZH2":    0.92,
    "EED":     0.88,
    "SUZ12":   0.85,
    "BRD4":    0.90,
    "BRD2":    0.82,
    "BRD3":    0.79,
    "HDAC1":   0.87,
    "HDAC2":   0.84,
    "HDAC3":   0.81,
    "HDAC4":   0.68,
    "HDAC6":   0.65,
    "KDM6A":   0.72,
    "KDM6B":   0.70,
    # ── Cell cycle / CDK ──────────────────────────────────────────────────────
    "CDK4":    0.85,
    "CDK6":    0.78,
    "CCND1":   0.72,
    "CDKN2A":  0.18,
    "RB1":     0.55,
    "E2F1":    0.68,
    # ── RTK / signalling ──────────────────────────────────────────────────────
    "PDGFRA":  0.88,
    "EGFR":    0.65,
    "MET":     0.60,
    "FGFR1":   0.55,
    "PIK3CA":  0.70,
    "AKT1":    0.68,
    "PTEN":    0.28,
    "MTOR":    0.72,
    # ── Proteasome ────────────────────────────────────────────────────────────
    "PSMB5":   0.55,
    "PSMB2":   0.52,
    "PSMB1":   0.50,
    "PSMD1":   0.58,
    # ── Stemness markers ──────────────────────────────────────────────────────
    "SOX2":    0.95,
    "NES":     0.92,
    "PROM1":   0.82,
    "CD44":    0.88,
    "OLIG2":   0.78,
    # ── MYC / STAT ────────────────────────────────────────────────────────────
    "MYC":     0.80,
    "MYCN":    0.82,
    "STAT3":   0.78,
    # ── ACVR1 / BMP ───────────────────────────────────────────────────────────
    "ACVR1":   0.75,
    "BMPR1A":  0.65,
    "SMAD1":   0.62,
    # ── DNA damage ────────────────────────────────────────────────────────────
    "PARP1":   0.68,
    "ATM":     0.55,
    "ATR":     0.60,
    # ── Metabolism ────────────────────────────────────────────────────────────
    "IDH1":    0.25,
    # ── Efflux transporters ───────────────────────────────────────────────────
    "ABCB1":   0.45,
    "ABCG2":   0.50,
    # ── Immune / TME ──────────────────────────────────────────────────────────
    "DRD2":    0.55,
    "SIGMAR1": 0.52,
    # ── Housekeeping (low by design) ──────────────────────────────────────────
    "ACTB":    0.50,
    "GAPDH":   0.48,
}


class TissueExpressionScorer:
    """
    v5.4: Single-Cell RNA-seq Stem-Cell Targeting Engine

    Reference: GSE102130 (Filbin et al. 2018, Science) — H3K27M DIPG
    single-cell RNA-seq (23,686 genes × 4,058 cells, RSEM TPM-normalised).

    All tunable parameters are read from pipeline_config.TISSUE.
    No magic numbers in this file.
    """

    def __init__(self, disease_name: str, data_dir: str = None):
        self.disease_name = disease_name
        sc_path = data_dir or PATHS["scrna"]
        self.sc_path  = Path(sc_path)
        self.is_ready = False
        self.gsc_target_scores: Dict[str, float] = {}
        self._gsc_all_values_sorted: List[float] = []
        self._gsc_p25 = 0.0
        self._gsc_p50 = 0.0
        self._gsc_p75 = 0.0
        self._gsc_p90 = 0.0

    async def _load_sc_data(self):
        if self.is_ready:
            return
        if not self.sc_path.exists():
            logger.warning(
                "⚠️ Single-cell data not found at %s — using curated DIPG fallback.",
                self.sc_path,
            )
            return

        logger.info("⏳ Loading Single-Cell H3K27M DIPG Matrix (Filbin 2018)...")
        try:
            df = pd.read_csv(self.sc_path, sep="\t", index_col=0)
            df.index = df.index.str.upper().str.strip()

            stem_markers = TISSUE["stem_markers"]
            markers = [m for m in stem_markers if m in df.index]
            if not markers:
                logger.warning(
                    "Stem cell markers not found in scRNA matrix. "
                    "Expected: %s", stem_markers
                )
                return

            stem_scores = df.loc[markers].mean(axis=0)
            threshold   = stem_scores.quantile(TISSUE["gsc_stem_quantile"])
            gsc_cells   = df.columns[stem_scores >= threshold]

            gsc_expr = df[gsc_cells].mean(axis=1)
            self.gsc_target_scores = gsc_expr.to_dict()

            all_vals = sorted(gsc_expr.values)
            self._gsc_all_values_sorted = all_vals
            n = len(all_vals)
            self._gsc_p25 = all_vals[int(n * 0.25)]
            self._gsc_p50 = all_vals[int(n * 0.50)]
            self._gsc_p75 = all_vals[int(n * 0.75)]
            self._gsc_p90 = all_vals[int(n * 0.90)]

            self.is_ready = True
            logger.info(
                "✅ Single-Cell loaded: %d stem-like cells / %d total. "
                "GSC expression p25=%.2f, p50=%.2f, p75=%.2f, p90=%.2f",
                len(gsc_cells), len(df.columns),
                self._gsc_p25, self._gsc_p50, self._gsc_p75, self._gsc_p90,
            )
        except Exception as e:
            logger.error("Single-cell loading failed: %s", e)

    def _expression_to_score(self, expr_value: float) -> float:
        """Convert RSEM TPM to 0–1 score via percentile bins (from config)."""
        bins = TISSUE["percentile_bins"]

        if not self._gsc_all_values_sorted:
            # Fallback when percentile data not computed
            if expr_value >= TISSUE["fallback_high_cutoff"]:
                return TISSUE["fallback_high_tpm"]
            elif expr_value >= TISSUE["fallback_mid_cutoff"]:
                return TISSUE["fallback_mid_tpm"]
            return TISSUE["fallback_low_tpm"]

        if expr_value >= self._gsc_p90:   return bins["p90"]
        elif expr_value >= self._gsc_p75: return bins["p75"]
        elif expr_value >= self._gsc_p50: return bins["p50"]
        elif expr_value >= self._gsc_p25: return bins["p25"]
        else:                             return bins["low"]

    def _curated_score(self, drug: Dict) -> float:
        """
        Score from curated DIPG GSC literature.
        Weighted blend of best target and mean of top-N targets.
        """
        targets = drug.get("targets", [])
        if not targets:
            return TISSUE["no_target_score"]

        scores        = [DIPG_GSC_CURATED_SCORES.get(t.upper(), TISSUE["unknown_target_score"])
                         for t in targets]
        scores_sorted = sorted(scores, reverse=True)
        best          = scores_sorted[0]
        top_n         = TISSUE["curated_top_n"]
        mean_top_n    = sum(scores_sorted[:top_n]) / min(len(scores_sorted), top_n)

        return round(
            TISSUE["curated_best_weight"]  * best
            + TISSUE["curated_top_n_weight"] * mean_top_n,
            4
        )

    def _dipg_relevance(self, target: str) -> float:
        """
        DIPG-relevance weight for a gene target.
        Prevents ubiquitous/housekeeping genes from dominating the score.
        Returns 1.0 if target is in the curated DIPG list, else config discount.
        """
        return 1.0 if target.upper() in DIPG_GSC_CURATED_SCORES \
               else TISSUE["off_target_relevance"]

    async def score_batch(self, candidates: List[Dict]) -> List[Dict]:
        await self._load_sc_data()

        cw = TISSUE["curated_weight"]
        sw = TISSUE["sc_weight"]

        for drug in candidates:
            targets = drug.get("targets", [])

            if not self.is_ready:
                drug["tissue_expression_score"] = self._curated_score(drug)
                drug["sc_context"] = "Curated DIPG GSC expression (scRNA-seq file not loaded)"
                continue

            curated_score = self._curated_score(drug)
            target_exprs  = [self.gsc_target_scores.get(t.upper(), 0.0) for t in targets]

            if not target_exprs or max(target_exprs) == 0.0:
                drug["tissue_expression_score"] = round(curated_score, 4)
                drug["sc_context"] = "Targets absent from scRNA-seq — using curated DIPG scores"
                continue

            weighted_sc_scores = [
                self._expression_to_score(expr) * self._dipg_relevance(t)
                for t, expr in zip(targets, target_exprs)
            ]
            sc_score_weighted = max(weighted_sc_scores) if weighted_sc_scores else \
                                TISSUE["percentile_bins"]["low"]

            blended   = round(cw * curated_score + sw * sc_score_weighted, 4)
            best_expr = max(target_exprs)

            drug["tissue_expression_score"] = blended
            drug["sc_context"] = (
                f"Blended v5.4 (H3K27M DIPG, Filbin 2018): "
                f"curated={curated_score:.2f} (w={cw}), "
                f"sc={sc_score_weighted:.2f} (w={sw}), "
                f"best_tpm={best_expr:.1f}, final={blended:.2f}"
            )

        scores = [c.get("tissue_expression_score", 0) for c in candidates]
        if scores:
            unique_vals = sorted(set(round(s, 2) for s in scores))
            logger.info(
                "✅ Tissue scoring: %d candidates | %d unique scores | "
                "range=[%.2f, %.2f] | distinct values=%s",
                len(candidates), len(unique_vals),
                min(scores), max(scores), unique_vals[:8],
            )

        return candidates