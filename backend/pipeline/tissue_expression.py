"""
tissue_expression.py  (v3.1 — import fixes)
============================================
Fixed: relative imports for probe_mapper and normal_brain_fetcher
moved to module top level so they work correctly as a package.
"""

import copy
import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

try:
    from .pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES
except ImportError:
    from pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES

# ── FIX: top-level imports instead of inside method bodies ──────────────────
try:
    from .probe_mapper import ProbeMapper, parse_geo_series_matrix
    from .normal_brain_fetcher import NormalBrainFetcher
except ImportError:
    try:
        from probe_mapper import ProbeMapper, parse_geo_series_matrix
        from normal_brain_fetcher import NormalBrainFetcher
    except ImportError:
        ProbeMapper = None
        parse_geo_series_matrix = None
        NormalBrainFetcher = None

logger = logging.getLogger(__name__)

GENE_EXPR_PATH   = Path("data/raw_omics/GSE70678_gene_expression.tsv")
SERIES_MTX_PATH  = Path("data/raw_omics/GSE70678_series_matrix.txt.gz")
GTEX_REF_PATH    = Path("data/raw_omics/GTEx_brain_normal_reference.tsv")


class TissueExpressionScorer:
    """ATRT tissue expression scorer using GSE70678 bulk RNA-seq."""

    def __init__(self, disease_name: str = "atrt", data_dir: str = None):
        self.disease_name  = disease_name
        self.is_ready      = False
        self.subgroup: Optional[str] = None
        self._diff_scores: Dict[str, float] = {}
        self._p25 = self._p50 = self._p75 = self._p90 = 0.0
        self._all_diff_sorted: List[float] = []
        self._n_atrt = self._n_normal = 0
        self._raw_gene_df: Optional[pd.DataFrame] = None
        self._gtex_ref: Optional[pd.Series]       = None

    def set_subgroup(self, subgroup: Optional[str]) -> None:
        if subgroup and subgroup.upper() in ("TYR", "SHH", "MYC"):
            self.subgroup = subgroup.upper()
            logger.info("TissueExpressionScorer: subgroup set to %s", self.subgroup)
        else:
            self.subgroup = None

    # ── Loading ───────────────────────────────────────────────────────────────

    async def _load_bulk_rna(self, quantile: float = None) -> None:
        if self.is_ready and quantile is None:
            return

        gene_df = self._load_gene_expression_matrix()
        if gene_df is None:
            logger.warning(
                "No GSE70678 data loaded — tissue scores will use curated values only.\n"
                "Run: python -m backend.pipeline.data_downloader --dataset gse70678 gpl6244\n"
                "     python -m backend.pipeline.data_downloader --process all"
            )
            return

        self._raw_gene_df = gene_df
        self._n_atrt      = len(gene_df.columns)

        gtex_ref = self._load_gtex_reference()
        if gtex_ref is None:
            logger.warning(
                "GTEx reference not loaded — using raw ATRT expression (no subtraction).\n"
                "Run: python -m backend.pipeline.data_downloader --dataset gtex_brain\n"
                "     python -m backend.pipeline.data_downloader --process gtex"
            )
            atrt_mean = gene_df.mean(axis=1)
            atrt_mean.index = atrt_mean.index.str.upper()
            diff = atrt_mean
        else:
            self._gtex_ref = gtex_ref
            gene_df_num = gene_df.apply(pd.to_numeric, errors="coerce")
            atrt_mean   = gene_df_num.mean(axis=1)
            atrt_mean.index = atrt_mean.index.str.upper()
            gtex_ref.index  = gtex_ref.index.str.upper()
            common = atrt_mean.index.intersection(gtex_ref.index)
            diff   = atrt_mean.loc[common] - gtex_ref.loc[common]
            self._n_normal = 4

        diff = diff.dropna()
        diff_sorted = sorted(diff.values)
        n = len(diff_sorted)

        self._diff_scores      = dict(zip(diff.index, diff.values))
        self._all_diff_sorted  = diff_sorted
        q = quantile or TISSUE["bulk_high_quantile"]
        self._p25 = diff_sorted[int(n * 0.25)] if n > 4 else 0.0
        self._p50 = diff_sorted[int(n * 0.50)] if n > 4 else 0.0
        self._p75 = diff_sorted[int(n * q)]    if n > 4 else 0.0
        self._p90 = diff_sorted[int(n * 0.90)] if n > 4 else 0.0

        if quantile is None:
            self.is_ready = True

        n_up   = int((diff > 0).sum())
        n_down = int((diff < 0).sum())
        logger.info(
            "GSE70678 differential loaded: %d genes (%d up, %d down vs GTEx)\n"
            "  Percentiles p25=%.2f p50=%.2f p75=%.2f p90=%.2f | "
            "ATRT n=%d | GTEx tissues n=%d",
            len(diff), n_up, n_down,
            self._p25, self._p50, self._p75, self._p90,
            self._n_atrt, self._n_normal,
        )

    def _load_gene_expression_matrix(self) -> Optional[pd.DataFrame]:
        # Option 1: Pre-processed TSV (fastest)
        if GENE_EXPR_PATH.exists():
            logger.info("Loading pre-processed gene expression: %s", GENE_EXPR_PATH)
            df = pd.read_csv(GENE_EXPR_PATH, sep="\t", index_col=0)
            df.index = df.index.astype(str).str.upper().str.strip()
            logger.info("Loaded %d genes × %d samples", len(df), len(df.columns))
            return df

        # Option 2: Map raw series matrix on the fly
        if SERIES_MTX_PATH.exists() and ProbeMapper is not None:
            logger.info("Mapping probes from: %s", SERIES_MTX_PATH)
            try:
                probe_df, _meta = parse_geo_series_matrix(SERIES_MTX_PATH)
                mapper = ProbeMapper()
                mapper.load()
                coverage = mapper.get_coverage(probe_df.index.tolist())
                logger.info(
                    "Probe coverage: %d/%d (%.1f%%)",
                    coverage["mapped_probes"], coverage["total_probes"],
                    coverage["coverage_pct"],
                )
                gene_df = mapper.map_probes_to_genes(probe_df, aggregation="median")
                gene_df.to_csv(GENE_EXPR_PATH, sep="\t")
                logger.info("Cached gene matrix → %s", GENE_EXPR_PATH)
                return gene_df
            except Exception as e:
                logger.error("On-the-fly probe mapping failed: %s", e)
                return None

        logger.warning(
            "Neither %s nor %s found. "
            "Run: python -m backend.pipeline.data_downloader --dataset gse70678 gpl6244 --process all",
            GENE_EXPR_PATH, SERIES_MTX_PATH,
        )
        return None

    def _load_gtex_reference(self) -> Optional[pd.Series]:
        if not GTEX_REF_PATH.exists():
            return None
        if NormalBrainFetcher is None:
            return None
        try:
            fetcher = NormalBrainFetcher(processed_path=GTEX_REF_PATH)
            fetcher.load()
            return fetcher.get_reference()
        except Exception as e:
            logger.error("Failed to load GTEx reference: %s", e)
            return None

    # ── Scoring ───────────────────────────────────────────────────────────────

    async def score_batch(self, candidates: List[Dict]) -> List[Dict]:
        await self._load_bulk_rna()
        return self._score_with_current_state(candidates)

    def _score_with_current_state(self, candidates: List[Dict]) -> List[Dict]:
        cw = TISSUE["curated_weight"]
        bw = TISSUE["bulk_weight"]

        for drug in candidates:
            curated = self._curated_score(drug)

            if not self.is_ready or not self._diff_scores:
                drug["tissue_expression_score"] = curated
                drug["sc_context"] = "Curated ATRT scores only (GSE70678 not loaded)"
                continue

            targets   = drug.get("targets", [])
            diff_vals = [
                self._diff_scores.get(t.upper())
                for t in targets
                if self._diff_scores.get(t.upper()) is not None
            ]

            if not diff_vals:
                drug["tissue_expression_score"] = curated
                drug["sc_context"] = "Targets not in GSE70678 — curated fallback"
                continue

            bulk_scores = [
                self._diff_to_score(v) * self._atrt_relevance(t)
                for t, v in zip(
                    [t for t in targets if self._diff_scores.get(t.upper()) is not None],
                    diff_vals,
                )
            ]
            bulk_score = max(bulk_scores)

            if self.subgroup:
                for t in targets:
                    mult = self._subgroup_multiplier(t, self.subgroup)
                    bulk_score = min(1.0, bulk_score * max(1.0, mult))

            blended = round(cw * curated + bw * bulk_score, 4)
            best_diff = max(diff_vals)

            drug["tissue_expression_score"] = blended
            drug["sc_context"] = (
                f"GSE70678 vs GTEx: curated={curated:.2f}(w={cw}), "
                f"bulk_diff={bulk_score:.2f}(w={bw}), "
                f"best_diff={best_diff:.2f}, final={blended:.2f}"
                + (f", subgroup={self.subgroup}" if self.subgroup else "")
            )

        scores = [c.get("tissue_expression_score", 0) for c in candidates]
        if scores:
            logger.info(
                "Tissue scoring: %d candidates | range=[%.2f, %.2f]",
                len(candidates), min(scores), max(scores),
            )
        return candidates

    def _diff_to_score(self, diff_value: float) -> float:
        bins = TISSUE["percentile_bins"]
        if not self._all_diff_sorted:
            if diff_value >= 2.0:   return bins["p90"]
            if diff_value >= 1.0:   return bins["p75"]
            if diff_value >= 0.0:   return bins["p50"]
            if diff_value >= -0.5:  return bins["p25"]
            return bins["low"]
        if diff_value >= self._p90:   return bins["p90"]
        if diff_value >= self._p75:   return bins["p75"]
        if diff_value >= self._p50:   return bins["p50"]
        if diff_value >= self._p25:   return bins["p25"]
        return bins["low"]

    def _curated_score(self, drug: Dict) -> float:
        targets = drug.get("targets", [])
        if not targets:
            return TISSUE["no_target_score"]
        scores = []
        for t in targets:
            base = ATRT_CURATED_SCORES.get(t.upper(), TISSUE["unknown_target_score"])
            if self.subgroup:
                base = min(1.0, base * self._subgroup_multiplier(t, self.subgroup))
            scores.append(base)
        scores_sorted = sorted(scores, reverse=True)
        best       = scores_sorted[0]
        top_n      = TISSUE["curated_top_n"]
        mean_top_n = sum(scores_sorted[:top_n]) / min(len(scores_sorted), top_n)
        return round(
            TISSUE["curated_best_weight"]    * best
            + TISSUE["curated_top_n_weight"] * mean_top_n,
            4,
        )

    @staticmethod
    def _atrt_relevance(target: str) -> float:
        return 1.0 if target.upper() in ATRT_CURATED_SCORES \
               else TISSUE["off_target_relevance"]

    @staticmethod
    def _subgroup_multiplier(gene: str, subgroup: str) -> float:
        g = gene.upper()
        tyr = {"TYR": 1.5, "DCT": 1.5, "MITF": 1.5, "SOX10": 1.4,
               "HDAC1": 1.2, "BRD4": 1.15}
        shh = {"GLI2": 1.5, "SMO": 1.4, "PTCH1": 1.3, "SHH": 1.3,
               "EZH2": 1.2, "CDK4": 1.15}
        myc = {"MYC": 1.6, "MYCN": 1.5, "AURKA": 1.4, "BRD4": 1.3,
               "CDK4": 1.2, "PSMB5": 1.2, "EZH2": 1.15}
        if subgroup == "TYR": return tyr.get(g, 1.0)
        if subgroup == "SHH": return shh.get(g, 1.0)
        if subgroup == "MYC": return myc.get(g, 1.0)
        return 1.0

    # ── Sensitivity analysis ──────────────────────────────────────────────────

    async def run_quantile_sensitivity(self, candidates: List[Dict]) -> Dict:
        if not SERIES_MTX_PATH.exists() and not GENE_EXPR_PATH.exists():
            return {
                "stable": None,
                "note":   "GSE70678 not found — sensitivity analysis skipped",
                "rankings_by_quantile": {},
            }
        quantiles = TISSUE["bulk_quantile_sensitivity_range"]
        rankings: Dict[float, List[str]] = {}
        for q in quantiles:
            await self._load_bulk_rna(quantile=q)
            scored = self._score_with_current_state(copy.deepcopy(candidates[:20]))
            scored.sort(key=lambda x: x.get("tissue_expression_score", 0), reverse=True)
            rankings[q] = [c.get("name", "?") for c in scored[:5]]
        self.is_ready = False
        await self._load_bulk_rna()
        all_top4 = [tuple(rankings[q][:4]) for q in quantiles if rankings.get(q)]
        stable   = len(set(all_top4)) == 1
        logger.info(
            "DE quantile sensitivity: top-4 %s across %s",
            "STABLE ✅" if stable else "UNSTABLE ⚠️", quantiles,
        )
        return {
            "stable":               stable,
            "quantiles_tested":     quantiles,
            "primary_quantile":     TISSUE["bulk_high_quantile"],
            "top4_consistent":      stable,
            "rankings_by_quantile": {str(q): rankings.get(q, []) for q in quantiles},
            "note": (
                "Top-4 stable ✅" if stable
                else "⚠️ Top-4 varies — interpret with caution."
            ),
        }