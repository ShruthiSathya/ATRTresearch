import asyncio
import logging
import math
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

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# CURATED DIPG RESISTANCE BYPASS NETWORK
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

ALL_RESISTANCE_NODES: Set[str] = set()
for bypass_targets in DIPG_RESISTANCE_BYPASS_MAP.values():
    ALL_RESISTANCE_NODES.update(bypass_targets)


def compute_escape_bypass_score(
    drug_targets: List[str],
    upregulated_genes: Set[str],
) -> float:
    """
    Compute escape bypass score using curated DIPG resistance pathway data.
    1.00 = no active bypass | 0.40 = 4+ active bypass routes
    """
    if not drug_targets:
        return 0.70

    active_bypass_nodes: Set[str] = set()
    covered_targets: Set[str] = set()

    for target in drug_targets:
        target_upper = target.upper()
        bypass_candidates = DIPG_RESISTANCE_BYPASS_MAP.get(target_upper, [])
        if bypass_candidates:
            covered_targets.add(target_upper)
            for node in bypass_candidates:
                is_tumor_upregulated = node in upregulated_genes
                is_constitutive_dipg_resistance = node in {
                    "PIK3CA", "MTOR", "MYC", "MYCN", "BCL2L1", "MCL1", "CDK4", "BRD4",
                }
                if is_tumor_upregulated or is_constitutive_dipg_resistance:
                    active_bypass_nodes.add(node)

    n_bypass  = len(active_bypass_nodes)
    n_covered = len(covered_targets)

    if n_covered == 0:
        return 0.75

    if n_bypass == 0:
        return 1.00
    elif n_bypass == 1:
        return 0.85
    elif n_bypass == 2:
        return 0.72
    elif n_bypass == 3:
        return 0.58
    else:
        return 0.40


# ─────────────────────────────────────────────────────────────────────────────
# PedcBioPortalValidator
# ─────────────────────────────────────────────────────────────────────────────

class PedcBioPortalValidator:
    def __init__(self, data_dir: str = "data/validation/cbtn_genomics/"):
        self.data_dir = data_dir

    def validate_triple_combo_cohort(self) -> Dict:
        try:
            mut_path = Path(f"{self.data_dir}mutations.txt")
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
                if c.upper() in {"H3-3A", "H3F3A", "DRD2", "HDAC1", "CDK4",
                                  "CDKN2A", "EGFR", "PTEN", "EZH2", "ACVR1"}
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
                "  Found columns: %s\n"
                "  Genomic validation will be skipped.",
                mut.columns.tolist(),
            )
            return {}

        except Exception as e:
            logger.warning("Genomic validation error: %s", e)
            return {}

    def _parse_standard_format(self, mut, hugo_col, hgvs_col, sample_col) -> Dict:
        h3k27m_samples = set(mut[
            (mut[hugo_col].str.upper().isin(["H3-3A", "H3F3A"])) &
            (mut[hgvs_col].str.contains("K28M|K27M", na=False, case=False))
        ][sample_col])
        return self._load_cna_and_rna(h3k27m_samples)

    def _parse_transposed_format(self, mut_path: Path) -> Dict:
        raw = pd.read_csv(mut_path, sep="\t", index_col=0)
        raw.index = raw.index.astype(str).str.upper().str.strip()

        H3_ALIASES = {"H3-3A", "H3F3A", "HIST1H3B", "H33A", "H3.3A"}
        h3_row = next((idx for idx in raw.index if idx in H3_ALIASES), None)

        if h3_row is None:
            logger.warning("Transposed format: H3K27M row not found. Skipping.")
            return {}

        h3_values      = raw.loc[h3_row].astype(str).str.upper().str.strip()
        H3K27M_VALUES  = {"K28M", "K27M"}
        h3k27m_samples = set(raw.columns[h3_values.isin(H3K27M_VALUES)])

        logger.info(
            "Transposed format: found %d H3K27M (K28M) samples out of %d total",
            len(h3k27m_samples), len(raw.columns),
        )
        return self._load_cna_and_rna(h3k27m_samples)

    def _load_cna_and_rna(self, h3k27m_samples: set) -> Dict:
        try:
            cna_path = Path(f"{self.data_dir}cna.txt")
            rna_path = Path(f"{self.data_dir}rna_zscores.txt")

            cdkn2a_del_samples = set()
            upregulated_genes  = set()
            total_samples      = max(len(h3k27m_samples), 1)

            if cna_path.exists():
                cna = pd.read_csv(cna_path, sep="\t", index_col=0)
                cna.index     = cna.index.astype(str).str.upper().str.strip()
                total_samples = len(cna.columns)

                cdkn2a_idx = next((idx for idx in cna.index if idx == "CDKN2A"), None)
                if cdkn2a_idx:
                    cdkn2a_row         = pd.to_numeric(cna.loc[cdkn2a_idx], errors="coerce")
                    cdkn2a_del_samples = set(cna.columns[cdkn2a_row <= -1])
                    logger.info(
                        "CDKN2A: %d samples with deletion (score ≤ -1), %d with deep deletion (score = -2)",
                        len(cdkn2a_del_samples), len(cna.columns[cdkn2a_row == -2])
                    )

            if rna_path.exists():
                rna         = pd.read_csv(rna_path, sep="\t", index_col=0)
                overlap_rna = list(h3k27m_samples.intersection(rna.columns))
                if overlap_rna:
                    high_exp          = rna[overlap_rna].apply(pd.to_numeric, errors="coerce").mean(axis=1)
                    upregulated_genes = set(high_exp[high_exp > 2.0].index)
                    logger.info(
                        "RNA: %d genes upregulated (z-score > 2.0) in H3K27M samples",
                        len(upregulated_genes),
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

        # p-value
        p_val = self._stat_validator.calculate_cooccurrence_p_value(genomic_stats)

        upregulated = genomic_stats.get("upregulated_genes", set())

        # ── Escape Route Analysis — curated DIPG bypass map ──────────────────
        for drug in candidates:
            targets = drug.get("targets", [])

            live_escape_hits = []
            for t in targets:
                neighbors = self._ppi.get_neighbors(t)
                live_escape_hits.extend([n for n in neighbors if n in upregulated])

            drug["resistance_nodes"] = list(set(live_escape_hits))

            curated_score = compute_escape_bypass_score(targets, upregulated)

            if live_escape_hits:
                string_score = 1.0 if not live_escape_hits else 0.4
                drug["escape_bypass_score"] = round(0.6 * string_score + 0.4 * curated_score, 4)
                drug["escape_note"] = "STRING-DB + curated DIPG bypass"
            else:
                drug["escape_bypass_score"] = round(curated_score, 4)
                drug["escape_note"] = "Curated DIPG resistance bypass (STRING-DB not available)"

        # Multi-omic scoring streams
        candidates = await self._tissue.score_batch(candidates)
        candidates = await self._depmap.score_batch(candidates, disease_name)
        candidates = await self._ppi.score_batch(candidates, disease_data["genes"])

        # Composite score
        for c in candidates:
            t_score = c.get("tissue_expression_score", 0.1)
            d_score = c.get("depmap_score", 0.1)
            p_score = c.get("ppi_score", 0.1)
            e_score = c.get("escape_bypass_score", 0.4)
            c["score"] = (t_score * 0.40) + (d_score * 0.30) + (e_score * 0.20) + (p_score * 0.10)

        # BBB filter
        candidates, _ = self._bbb_filter.filter_and_rank(candidates, apply_penalty=True)

        sorted_candidates = sorted(candidates, key=lambda x: x.get("score", 0), reverse=True)

        hypotheses = self._hyp_gen.generate(
            candidates        = sorted_candidates[:top_k],
            cmap_results      = [],
            synergy_combos    = [],
            differential_cmap = [],
            genomic_stats     = genomic_stats,
            p_value           = p_val,
        )

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
                # ── FIX v5.3 ─────────────────────────────────────────────────
                # len(sorted_candidates) = full OpenTargets count BEFORE [:top_k]
                # slice. save_results.py reads this as stats["n_screened"].
                # Without this fix, n_drugs_screened was always == top_k (20).
                "n_screened":          len(sorted_candidates),
            },
        }