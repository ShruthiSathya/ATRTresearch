"""
tissue_expression.py  (v4.0)
============================
ATRT tissue expression scorer using GSE70678 bulk RNA-seq.

KEY FIXES FROM v3.1
--------------------
1. All file paths now come from pipeline_config.PATHS (absolute).
   The old hardcoded relative paths "data/raw_omics/..." broke when the
   pipeline was launched from any directory other than the repo root.

2. GPL570 probe map loading: the file at PATHS["gpl570_map"] is a TSV
   with columns [probe_id, gene_symbol, entrez_id].  The old code tried
   to instantiate ProbeMapper and call load() on it, which re-downloads
   GPL570.annot.gz even when the processed TSV already exists.
   Fixed: load the TSV directly — if it exists, skip ProbeMapper entirely.

3. GTEx reference loading: reads PATHS["gtex_ref"] directly as a TSV.
   The old NormalBrainFetcher wrapper added complexity without benefit.
   If the file is not present, falls back to curated brain reference.

4. Differential expression: uses simple mean(ATRT) − mean(GTEx_brain).
   Both datasets are log2-scale so the difference is a valid log2FC proxy.
   This is then converted to a 0–1 score via percentile binning.

5. Curated score fallback: always available even if no data files exist.
   Scores come from ATRT_CURATED_SCORES in pipeline_config (literature-based).
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

try:
    from .pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES
except ImportError:
    from pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES

logger = logging.getLogger(__name__)


class TissueExpressionScorer:
    """ATRT tissue expression scorer using GSE70678 bulk RNA-seq."""

    def __init__(self, disease_name: str = "atrt", data_dir: str = None):
        self.disease_name = disease_name
        self.subgroup: Optional[str] = None
        self.is_ready = False

        # Differential expression scores: gene → float (ATRT minus normal)
        self._diff_scores: Dict[str, float] = {}

        # Percentile thresholds computed from the diff score distribution
        self._p25 = self._p50 = self._p75 = self._p90 = 0.0
        self._all_diff_sorted: List[float] = []

        # File paths (absolute, from config)
        self._gene_expr_path = Path(PATHS["scrna"])
        self._gtex_ref_path  = Path(PATHS["gtex_ref"])

    def set_subgroup(self, subgroup: Optional[str]) -> None:
        if subgroup and subgroup.upper() in ("TYR", "SHH", "MYC"):
            self.subgroup = subgroup.upper()
            logger.info("TissueExpressionScorer: subgroup=%s", self.subgroup)
        else:
            self.subgroup = None

    # ── Loading ───────────────────────────────────────────────────────────────

    async def _load_bulk_rna(self) -> None:
        if self.is_ready:
            return

        gene_df  = self._load_gene_expression_matrix()
        gtex_ref = self._load_gtex_reference()

        if gene_df is None:
            logger.warning(
                "GSE70678 gene expression matrix not found at: %s\n"
                "Tissue scores will use curated values only.\n"
                "To add live data: place the file at the path above.",
                self._gene_expr_path,
            )
            self.is_ready = True
            return

        # Compute mean ATRT expression across all samples
        gene_df_num = gene_df.apply(pd.to_numeric, errors="coerce")
        atrt_mean   = gene_df_num.mean(axis=1).dropna()
        atrt_mean.index = atrt_mean.index.str.upper().str.strip()

        if gtex_ref is not None:
            gtex_ref.index = gtex_ref.index.str.upper().str.strip()
            common = atrt_mean.index.intersection(gtex_ref.index)
            if len(common) > 100:
                diff = (atrt_mean.loc[common] - gtex_ref.loc[common]).dropna()
                logger.info(
                    "GSE70678 differential: %d genes in common with GTEx "
                    "(%d ATRT samples vs %d reference tissues)",
                    len(diff), gene_df_num.shape[1], gtex_ref.shape[0]
                    if hasattr(gtex_ref, 'shape') else "?",
                )
            else:
                logger.warning(
                    "Only %d genes in common with GTEx — using raw ATRT expression.",
                    len(common)
                )
                diff = atrt_mean
        else:
            logger.warning("GTEx reference not found — using raw ATRT expression (no subtraction).")
            diff = atrt_mean

        diff_sorted = sorted(diff.values)
        n = len(diff_sorted)
        if n < 10:
            logger.warning("Too few genes in diff (%d) — curated scores only.", n)
            self.is_ready = True
            return

        self._diff_scores     = dict(zip(diff.index, diff.values))
        self._all_diff_sorted = diff_sorted
        q = TISSUE["bulk_high_quantile"]
        self._p25 = diff_sorted[max(0, int(n * 0.25))]
        self._p50 = diff_sorted[max(0, int(n * 0.50))]
        self._p75 = diff_sorted[max(0, int(n * q))]
        self._p90 = diff_sorted[max(0, int(n * 0.90))]

        self.is_ready = True
        n_up   = int((diff > 0).sum())
        n_down = int((diff < 0).sum())
        logger.info(
            "✅ GSE70678 loaded: %d genes | %d up, %d down vs normal brain\n"
            "   Percentiles: p25=%.2f p50=%.2f p75=%.2f p90=%.2f",
            n, n_up, n_down,
            self._p25, self._p50, self._p75, self._p90,
        )

    def _load_gene_expression_matrix(self) -> Optional[pd.DataFrame]:
        """Load pre-processed gene × sample expression matrix."""
        if self._gene_expr_path.exists():
            logger.info("Loading GSE70678 gene expression: %s", self._gene_expr_path)
            try:
                df = pd.read_csv(str(self._gene_expr_path), sep="\t", index_col=0, low_memory=False)
                df.index = df.index.astype(str).str.upper().str.strip()
                logger.info("Loaded %d genes × %d samples", len(df), len(df.columns))
                return df
            except Exception as e:
                logger.error("Failed to load GSE70678: %s", e)
                return None

        # Check if raw series matrix exists (will need probe mapping)
        raw_path = Path(PATHS.get("gse70678_raw", ""))
        map_path = Path(PATHS.get("gpl570_map", ""))

        if raw_path.exists() and map_path.exists():
            logger.info("Found raw GSE70678 — mapping probes to genes...")
            return self._map_probes_to_genes(raw_path, map_path)

        return None

    def _map_probes_to_genes(self, series_matrix: Path, probe_map: Path) -> Optional[pd.DataFrame]:
        """Map Affymetrix probe IDs to gene symbols using GPL570 probe map."""
        try:
            # Load probe map
            pm = pd.read_csv(str(probe_map), sep="\t", dtype=str, na_filter=False)
            pm.columns = pm.columns.str.strip().str.lower()
            probe_col  = next((c for c in pm.columns if "probe" in c), pm.columns[0])
            symbol_col = next((c for c in pm.columns if "gene" in c or "symbol" in c), pm.columns[1])
            probe_to_gene = dict(zip(pm[probe_col], pm[symbol_col].str.upper()))

            # Parse series matrix
            try:
                from .probe_mapper import parse_geo_series_matrix
            except ImportError:
                from probe_mapper import parse_geo_series_matrix

            probe_df, _ = parse_geo_series_matrix(series_matrix)

            # Map probes → genes
            probe_df.index = probe_df.index.astype(str).str.strip()
            genes = probe_df.index.map(lambda p: probe_to_gene.get(p, ""))
            probe_df = probe_df[genes != ""]
            probe_df.index = [probe_to_gene[p] for p in probe_df.index if probe_to_gene.get(p)]

            # Aggregate by gene (median across probes)
            gene_df = probe_df.apply(pd.to_numeric, errors="coerce").groupby(
                level=0, sort=False
            ).median()

            # Cache for next run
            self._gene_expr_path.parent.mkdir(parents=True, exist_ok=True)
            gene_df.to_csv(str(self._gene_expr_path), sep="\t")
            logger.info("Cached gene expression matrix → %s", self._gene_expr_path)
            return gene_df

        except Exception as e:
            logger.error("Probe mapping failed: %s", e, exc_info=True)
            return None

    def _load_gtex_reference(self) -> Optional[pd.Series]:
        """Load GTEx normal brain reference as a gene → mean_log2_expr Series."""
        if not self._gtex_ref_path.exists():
            return None
        try:
            df = pd.read_csv(str(self._gtex_ref_path), sep="\t", index_col=0, low_memory=False)
            df.index = df.index.astype(str).str.upper().str.strip()

            # Select brain-relevant columns if multiple exist
            brain_kws = ["brain", "cerebellum", "cortex", "hippocampus"]
            brain_cols = [c for c in df.columns
                          if any(k.lower() in c.lower() for k in brain_kws)]
            if brain_cols:
                ref = df[brain_cols].mean(axis=1)
            else:
                # Single-column TSV (already a reference vector)
                ref = df.iloc[:, 0]

            logger.info("GTEx reference loaded: %d genes, %d brain tissue columns",
                        len(ref), len(brain_cols) if brain_cols else 1)
            return ref
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
        has_bulk = bool(self._diff_scores)

        for drug in candidates:
            curated = self._curated_score(drug)

            if not has_bulk:
                drug["tissue_expression_score"] = curated
                drug["sc_context"] = "Curated scores only (GSE70678 not loaded)"
                continue

            targets = [t.upper() for t in (drug.get("targets") or [])]
            diff_vals = {
                t: self._diff_scores[t]
                for t in targets if t in self._diff_scores
            }

            if not diff_vals:
                drug["tissue_expression_score"] = curated
                drug["sc_context"] = "Targets not in GSE70678 — curated fallback"
                continue

            # Best + top-2 mean aggregation
            bulk_scores = sorted(
                [self._diff_to_score(v) * self._atrt_relevance(t)
                 for t, v in diff_vals.items()],
                reverse=True,
            )
            bulk_score = 0.65 * bulk_scores[0] + 0.35 * (
                sum(bulk_scores[:2]) / min(len(bulk_scores), 2)
            )

            # Subgroup multiplier
            if self.subgroup:
                for t in targets:
                    mult = self._subgroup_multiplier(t, self.subgroup)
                    if mult > 1.0:
                        bulk_score = min(1.0, bulk_score * mult)

            blended = round(cw * curated + bw * bulk_score, 4)
            drug["tissue_expression_score"] = blended
            drug["sc_context"] = (
                f"GSE70678: curated={curated:.2f}×{cw}, "
                f"bulk={bulk_score:.2f}×{bw}, final={blended:.2f}"
                + (f", subgroup={self.subgroup}" if self.subgroup else "")
            )

        scores = [c.get("tissue_expression_score", 0) for c in candidates]
        if scores:
            logger.info(
                "Tissue scoring: %d candidates | range=[%.2f, %.2f]",
                len(candidates), min(scores), max(scores)
            )
        return candidates

    def _diff_to_score(self, diff_value: float) -> float:
        bins = TISSUE["percentile_bins"]
        if not self._all_diff_sorted:
            # No real data — use magnitude heuristic
            if diff_value >= 2.0:   return bins["p90"]
            if diff_value >= 1.0:   return bins["p75"]
            if diff_value >= 0.0:   return bins["p50"]
            if diff_value >= -0.5:  return bins["p25"]
            return bins["low"]
        if diff_value >= self._p90: return bins["p90"]
        if diff_value >= self._p75: return bins["p75"]
        if diff_value >= self._p50: return bins["p50"]
        if diff_value >= self._p25: return bins["p25"]
        return bins["low"]

    def _curated_score(self, drug: Dict) -> float:
        targets = drug.get("targets") or []
        if not targets:
            return TISSUE["no_target_score"]
        scores = []
        for t in targets:
            base = ATRT_CURATED_SCORES.get(t.upper(), TISSUE["unknown_target_score"])
            if self.subgroup:
                base = min(1.0, base * self._subgroup_multiplier(t, self.subgroup))
            scores.append(base)
        if not scores:
            return TISSUE["no_target_score"]
        scores_sorted = sorted(scores, reverse=True)
        n    = TISSUE["curated_top_n"]
        best = scores_sorted[0]
        mean_top_n = sum(scores_sorted[:n]) / min(len(scores_sorted), n)
        return round(
            TISSUE["curated_best_weight"] * best
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
        tyr_mults = {"TYR": 1.5, "DCT": 1.5, "MITF": 1.5, "SOX10": 1.4,
                     "HDAC1": 1.2, "BRD4": 1.15}
        shh_mults = {"GLI2": 1.5, "SMO": 1.4, "PTCH1": 1.3, "SHH": 1.3,
                     "EZH2": 1.2, "CDK4": 1.15}
        myc_mults = {"MYC": 1.6, "MYCN": 1.5, "AURKA": 1.4, "BRD4": 1.3,
                     "CDK4": 1.2, "PSMB5": 1.2, "EZH2": 1.15}
        if subgroup == "TYR": return tyr_mults.get(g, 1.0)
        if subgroup == "SHH": return shh_mults.get(g, 1.0)
        if subgroup == "MYC": return myc_mults.get(g, 1.0)
        return 1.0