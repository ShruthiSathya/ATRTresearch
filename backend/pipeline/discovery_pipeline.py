import asyncio
import logging
import pandas as pd
from typing import Dict, List, Optional, Set
from pathlib import Path

from .data_fetcher import ProductionDataFetcher
from .cmap_query import CMAPQuery
from .ppi_network import PPINetwork
from .depmap_essentiality import DepMapEssentiality
from .synergy_predictor import SynergyPredictor
from .hypothesis_generator import HypothesisGenerator
from .statistical_validator import StatisticalValidator
from .tissue_expression import TissueExpressionScorer
from .bbb_filter import BBBFilter
from .pipeline_config import (
    COMPOSITE_WEIGHTS,
    SCORE_DEFAULTS,
    ESCAPE,
    GENOMICS,
    PATHS,
)

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# CURATED DIPG RESISTANCE BYPASS NETWORK
# Sources: Grasso 2015, Nagaraja 2017, Filbin 2018, Mackay 2017
# ─────────────────────────────────────────────────────────────────────────────

DIPG_RESISTANCE_BYPASS_MAP: Dict[str, List[str]] = {
    "CDK4":   ["PIK3CA", "AKT1", "MTOR", "CCNE1"],
    "CDK6":   ["PIK3CA", "AKT1", "MTOR"],
    "HDAC1":  ["BCL2", "BCL2L1", "MCL1", "MYC"],
    "HDAC2":  ["BCL2", "BCL2L1", "MCL1"],
    "HDAC3":  ["MYC", "NFKB1", "BCL2"],
    "EZH2":   ["BRD4", "BRD2", "KDM6A", "KDM6B"],
    "BRD4":   ["CDK4", "MYC", "MYCN"],
    "BRD2":   ["MYC", "CDK4"],
    "BRD3":   ["MYC"],
    "EGFR":   ["PIK3CA", "KRAS", "MET", "AXL"],
    "PDGFRA": ["PIK3CA", "MTOR", "AKT1"],
    "MET":    ["EGFR", "AXL", "PIK3CA"],
    "PIK3CA": ["MAPK1", "MTOR", "AKT1"],
    "PIK3CD": ["MAPK1", "MTOR", "AKT1"],
    "PIK3CG": ["MAPK1", "MTOR", "AKT1"],
    "MTOR":   ["PIK3CA", "AKT1"],
    "AKT1":   ["MTOR", "BCL2"],
    "PARP1":  ["RAD51", "BRCA1", "ATM"],
    "PSMB5":  ["ATG5", "BECN1", "SQSTM1"],
    "PSMB2":  ["ATG5", "BECN1"],
    "PSMB1":  ["BECN1"],
    "ACVR1":  ["PDGFRA", "EGFR", "PIK3CA"],
    "MDM2":   ["BCL2L1", "MCL1"],
    "BCL2":   ["MCL1", "BCL2L1"],
    "BCL2L1": ["MCL1", "BCL2"],
    "ATR":    ["RAD51", "CHEK1", "ATM"],
    "ATM":    ["ATR", "CHEK1"],
    "WEE1":   ["CDK4", "CCND1"],
    "H3F3A":  [],
    "STAT3":  ["JAK1", "JAK2"],
    "SOX2":   [],
    "OLIG2":  [],
}


def compute_escape_bypass_score(
    drug_targets: List[str],
    upregulated_genes: Set[str],
) -> float:
    """
    Compute escape bypass score using curated DIPG resistance pathway data.
    Scores per number of active bypass nodes are read from pipeline_config.ESCAPE.
    Constitutive resistance nodes are also from config.
    """
    if not drug_targets:
        return ESCAPE["empty_target_score"]

    constitutive = ESCAPE["constitutive_resistance_nodes"]
    active_bypass_nodes: Set[str] = set()
    covered_targets:     Set[str] = set()

    for target in drug_targets:
        target_upper      = target.upper()
        bypass_candidates = DIPG_RESISTANCE_BYPASS_MAP.get(target_upper, [])
        if bypass_candidates:
            covered_targets.add(target_upper)
            for node in bypass_candidates:
                if node in upregulated_genes or node in constitutive:
                    active_bypass_nodes.add(node)

    if not covered_targets:
        return ESCAPE["no_target_score"]

    n_bypass     = len(active_bypass_nodes)
    bypass_scores = ESCAPE["bypass_scores"]
    # Use the highest key that n_bypass meets or exceeds
    max_key = max(bypass_scores.keys())
    score   = bypass_scores.get(n_bypass, bypass_scores[max_key])
    return score


# ─────────────────────────────────────────────────────────────────────────────
# PedcBioPortalValidator
# ─────────────────────────────────────────────────────────────────────────────

class PedcBioPortalValidator:
    def __init__(self, data_dir: str = None):
        self.data_dir = data_dir or PATHS["genomics"]

    def validate_triple_combo_cohort(self) -> Dict:
        try:
            mut_path = Path(self.data_dir) / "mutations.txt"
            if not mut_path.exists():
                logger.info("Genomic validation files not present — skipping (non-fatal)")
                return {}

            mut = pd.read_csv(mut_path, sep="\t")
            mut.columns = mut.columns.str.lower().str.strip()

            hugo_col   = next((c for c in mut.columns if c in ["hugo_symbol", "gene", "symbol"]), None)
            hgvs_col   = next((c for c in mut.columns if c in ["hgvsp_short", "protein_change",
                                                                 "amino_acid_change", "mutation"]), None)
            sample_col = next((c for c in mut.columns if c in ["tumor_sample_barcode",
                                                                 "sample_id", "sample"]), None)

            if hugo_col and hgvs_col and sample_col:
                return self._parse_standard_format(mut, hugo_col, hgvs_col, sample_col)

            gene_like_cols = [
                c for c in mut.columns
                if c.upper() in GENOMICS["h3_gene_aliases"]
                   | {"DRD2", "HDAC1", "CDK4", "CDKN2A", "EGFR", "PTEN", "EZH2", "ACVR1"}
                or (len(c) <= 8 and c.replace("-", "").replace("_", "").isalnum())
            ]

            if gene_like_cols:
                logger.info(
                    "Detected transposed genomic format (gene columns: %s...). "
                    "Parsing as sample×gene matrix.",
                    ", ".join(gene_like_cols[:5]),
                )
                return self._parse_transposed_format(mut_path)

            logger.warning(
                "Genomic columns not recognised in mutations.txt.\n"
                "  Found columns: %s\n  Genomic validation will be skipped.",
                mut.columns.tolist(),
            )
            return {}

        except Exception as e:
            logger.warning("Genomic validation error: %s", e)
            return {}

    def _parse_standard_format(self, mut, hugo_col, hgvs_col, sample_col) -> Dict:
        h3k27m_values = GENOMICS["h3k27m_values"]
        h3_aliases    = GENOMICS["h3_gene_aliases"]
        h3k27m_samples = set(mut[
            (mut[hugo_col].str.upper().isin(h3_aliases)) &
            (mut[hgvs_col].str.contains("|".join(h3k27m_values), na=False, case=False))
        ][sample_col])
        return self._load_cna_and_rna(h3k27m_samples)

    def _parse_transposed_format(self, mut_path: Path) -> Dict:
        raw       = pd.read_csv(mut_path, sep="\t", index_col=0)
        raw.index = raw.index.astype(str).str.upper().str.strip()

        h3_aliases = GENOMICS["h3_gene_aliases"]
        h3_row     = next((idx for idx in raw.index if idx in h3_aliases), None)

        if h3_row is None:
            logger.warning("Transposed format: H3K27M row not found. Skipping.")
            return {}

        h3_values      = raw.loc[h3_row].astype(str).str.upper().str.strip()
        h3k27m_samples = set(raw.columns[h3_values.isin(GENOMICS["h3k27m_values"])])

        logger.info(
            "Transposed format: found %d H3K27M (K28M) samples out of %d total",
            len(h3k27m_samples), len(raw.columns),
        )
        return self._load_cna_and_rna(h3k27m_samples)

    def _load_cna_and_rna(self, h3k27m_samples: set) -> Dict:
        try:
            data_dir      = Path(self.data_dir)
            cna_path      = data_dir / "cna.txt"
            rna_path      = data_dir / "rna_zscores.txt"
            rna_threshold = GENOMICS["rna_upregulation_zscore"]
            cna_threshold = GENOMICS["cna_deletion_threshold"]

            cdkn2a_del_samples = set()
            upregulated_genes  = set()
            total_samples      = max(len(h3k27m_samples), 1)

            if cna_path.exists():
                cna       = pd.read_csv(cna_path, sep="\t", index_col=0)
                cna.index = cna.index.astype(str).str.upper().str.strip()
                total_samples = len(cna.columns)

                cdkn2a_idx = next((idx for idx in cna.index if idx == "CDKN2A"), None)
                if cdkn2a_idx:
                    cdkn2a_row         = pd.to_numeric(cna.loc[cdkn2a_idx], errors="coerce")
                    cdkn2a_del_samples = set(cna.columns[cdkn2a_row <= cna_threshold])
                    logger.info(
                        "CDKN2A: %d samples with deletion (score ≤ %d), "
                        "%d with deep deletion (score = -2)",
                        len(cdkn2a_del_samples), cna_threshold,
                        len(cna.columns[cdkn2a_row == -2])
                    )

            if rna_path.exists():
                rna = pd.read_csv(rna_path, sep="\t", index_col=0)

                # ── Validate this is a full expression matrix ─────────────────
                metadata_rows     = GENOMICS["rna_metadata_rows"]
                min_genes         = GENOMICS["rna_min_genes_required"]
                rna_threshold     = GENOMICS["rna_upregulation_zscore"]
                h3k27m_indicators = GENOMICS["rna_h3k27m_col_indicators"]
                normal_indicators = GENOMICS["rna_normal_col_indicators"]

                gene_rows = [r for r in rna.index
                             if str(r).upper().strip() not in metadata_rows]

                if len(gene_rows) < min_genes:
                    logger.warning(
                        "⚠️  rna_zscores.txt contains only %d gene rows — "
                        "clinical summary file, not an expression matrix. "
                        "RNA upregulation skipped. Escape bypass uses "
                        "constitutive DIPG resistance nodes only.",
                        len(gene_rows),
                    )
                else:
                    rna_expr    = rna.loc[gene_rows]
                    overlap_rna = list(h3k27m_samples.intersection(rna_expr.columns))

                    if overlap_rna:
                        # ── Matched cohort mode: sample IDs align ─────────────
                        high_exp          = rna_expr[overlap_rna].apply(
                            pd.to_numeric, errors="coerce"
                        ).mean(axis=1)
                        upregulated_genes = set(high_exp[high_exp > rna_threshold].index)
                        logger.info(
                            "RNA (matched cohort): %d genes upregulated "
                            "(value > %.1f) in %d H3K27M samples "
                            "(from %d-gene matrix)",
                            len(upregulated_genes), rna_threshold,
                            len(overlap_rna), len(gene_rows),
                        )
                    else:
                        # ── Reference cohort mode: different dataset (e.g. GSE115397)
                        # Identify H3K27M tumour columns by name pattern
                        h3k27m_ref_cols = [
                            c for c in rna_expr.columns
                            if any(ind.lower() in c.lower()
                                   for ind in h3k27m_indicators)
                            and not any(n.lower() in c.lower()
                                        for n in normal_indicators)
                        ]
                        normal_ref_cols = [
                            c for c in rna_expr.columns
                            if any(n.lower() in c.lower()
                                   for n in normal_indicators)
                        ]

                        if h3k27m_ref_cols:
                            high_exp = rna_expr[h3k27m_ref_cols].apply(
                                pd.to_numeric, errors="coerce"
                            ).mean(axis=1)
                            upregulated_genes = set(
                                high_exp[high_exp > rna_threshold].index
                            )
                            logger.info(
                                "RNA (reference cohort mode): %d genes upregulated "
                                "(value > %.1f) in %d H3K27M reference samples "
                                "(%s...) from %d-gene matrix. "
                                "Normal columns excluded: %s",
                                len(upregulated_genes), rna_threshold,
                                len(h3k27m_ref_cols),
                                ", ".join(h3k27m_ref_cols[:3]),
                                len(gene_rows),
                                normal_ref_cols,
                            )
                        else:
                            logger.warning(
                                "RNA: no sample ID overlap and no H3K27M "
                                "indicator columns found in rna_zscores.txt. "
                                "Columns present: %s. "
                                "RNA upregulation skipped.",
                                rna_expr.columns.tolist()[:5],
                            )

            overlap = h3k27m_samples.intersection(cdkn2a_del_samples)
            logger.info(
                "Co-occurrence: H3K27M=%d, CDKN2A-del=%d, overlap=%d, total=%d",
                len(h3k27m_samples), len(cdkn2a_del_samples), len(overlap), total_samples,
            )

            return {
                "h3k27m_count":      len(h3k27m_samples),
                "cdkn2a_del_count":  len(cdkn2a_del_samples),
                "overlap_count":     len(overlap),
                "prevalence":        len(overlap) / max(len(h3k27m_samples), 1),
                "upregulated_genes": upregulated_genes,
                "total_samples":     total_samples,
            }
        except Exception as e:
            logger.warning("CNA/RNA loading error: %s", e)
            return {}


# ─────────────────────────────────────────────────────────────────────────────
# ProductionPipeline
# ─────────────────────────────────────────────────────────────────────────────

class ProductionPipeline:
    def __init__(self):
        self._data_fetcher      = ProductionDataFetcher()
        self._cmap              = CMAPQuery()
        self._ppi               = PPINetwork()
        self._depmap            = DepMapEssentiality()
        self._tissue            = TissueExpressionScorer("dipg")
        self._synergy           = SynergyPredictor()
        self._hyp_gen           = HypothesisGenerator()
        self._genomic_validator = PedcBioPortalValidator()
        self._stat_validator    = StatisticalValidator()
        self._bbb_filter        = BBBFilter()

    async def initialize(self, disease: str):
        return True

    async def run(self, disease_name: str, top_k: int = 15) -> Dict:
        disease_data = await self._data_fetcher.fetch_disease_data(disease_name)
        candidates   = await self._data_fetcher.fetch_approved_drugs()

        if not candidates:
            logger.warning("⚠️ API returned no drugs. Using fallback safety library.")
            candidates = [
                {"name": "ONC201",       "targets": ["DRD2", "CLPB"]},
                {"name": "Panobinostat", "targets": ["HDAC1", "HDAC2"]},
                {"name": "Abemaciclib",  "targets": ["CDK4", "CDK6"]},
            ]

        # Genomic validation
        genomic_stats = self._genomic_validator.validate_triple_combo_cohort()
        p_val         = self._stat_validator.calculate_cooccurrence_p_value(genomic_stats)
        upregulated   = genomic_stats.get("upregulated_genes", set())

        # ── Escape Route Analysis ─────────────────────────────────────────────
        ew = ESCAPE
        for drug in candidates:
            targets = drug.get("targets", [])

            live_escape_hits = []
            for t in targets:
                neighbors = self._ppi.get_neighbors(t)
                live_escape_hits.extend([n for n in neighbors if n in upregulated])

            drug["resistance_nodes"] = list(set(live_escape_hits))
            curated_score = compute_escape_bypass_score(targets, upregulated)

            if live_escape_hits:
                drug["escape_bypass_score"] = round(
                    ew["string_weight"]  * ew["string_hit_score"]
                    + ew["curated_weight"] * curated_score,
                    4
                )
                drug["escape_note"] = "STRING-DB + curated DIPG bypass"
            else:
                drug["escape_bypass_score"] = round(curated_score, 4)
                drug["escape_note"]         = "Curated DIPG resistance bypass"

        # Multi-omic scoring
        candidates = await self._tissue.score_batch(candidates)
        candidates = await self._depmap.score_batch(candidates, disease_name)
        candidates = await self._ppi.score_batch(candidates, disease_data["genes"])

        # Composite score — weights and defaults from config
        cw = COMPOSITE_WEIGHTS
        sd = SCORE_DEFAULTS
        for c in candidates:
            t_score = c.get("tissue_expression_score", sd["tissue_expression_score"])
            d_score = c.get("depmap_score",            sd["depmap_score"])
            p_score = c.get("ppi_score",               sd["ppi_score"])
            e_score = c.get("escape_bypass_score",     sd["escape_bypass_score"])
            c["score"] = round(
                t_score * cw["tissue"]
                + d_score * cw["depmap"]
                + e_score * cw["escape"]
                + p_score * cw["ppi"],
                4
            )

        # BBB filter
        candidates, _ = self._bbb_filter.filter_and_rank(candidates, apply_penalty=True)
        sorted_candidates = sorted(candidates, key=lambda x: x.get("score", 0), reverse=True)

        # Generate hypotheses
        # NOTE: toxicity penalty is computed INSIDE hypothesis_generator only.
        # Do NOT call combination_toxicity_penalty here — it would duplicate the log.
        hypotheses = self._hyp_gen.generate(
            candidates        = sorted_candidates[:top_k],
            cmap_results      = [],
            synergy_combos    = [],
            differential_cmap = [],
            genomic_stats     = genomic_stats,
            p_value           = p_val,
        )

        # Extract toxicity stats from the generated hypothesis for the stats dict
        tox_breakdown = {}
        if hypotheses and "confidence_breakdown" in hypotheses[0]:
            bd = hypotheses[0]["confidence_breakdown"]
            tox_breakdown = {
                "toxicity_flag":       bd.get("toxicity_flag", "UNKNOWN"),
                "toxicity_multiplier": bd.get("toxicity_multiplier", 1.0),
                "toxicity_note":       bd.get("toxicity_note", ""),
            }

        return {
            "hypotheses":     hypotheses,
            "top_candidates": sorted_candidates[:top_k],
            "stats": {
                "p_value":             p_val,
                "p_value_label":       self._stat_validator.format_p_value_for_report(p_val),
                "genomic_data_loaded": bool(genomic_stats),
                "h3k27m_count":        genomic_stats.get("h3k27m_count", 0),
                "cdkn2a_del_count":    genomic_stats.get("cdkn2a_del_count", 0),
                "overlap_count":       genomic_stats.get("overlap_count", 0),
                "total_samples":       genomic_stats.get("total_samples", 0),
                "n_screened":          len(sorted_candidates),
                "composite_weights":   cw,
                **tox_breakdown,
            },
        }