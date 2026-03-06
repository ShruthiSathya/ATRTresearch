
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


class PedcBioPortalValidator:
    """
    Validates drug combinations against PedcBioPortal DIPG/GBM genomic data.

    Supports two file formats:
      1. Standard cBioPortal format: rows=samples, columns include hugo_symbol,
         hgvsp_short, tumor_sample_barcode (most TCGA/CBTN datasets)
      2. Transposed study format: index=gene names (H3-3A, CDKN2A etc.),
         columns=sample IDs (PNOC/PBTA format from cBioPortal download)
    """

    def __init__(self, data_dir: str = "data/validation/cbtn_genomics/"):
        self.data_dir = data_dir

    def validate_triple_combo_cohort(self) -> Dict:
        try:
            mut_path = Path(f"{self.data_dir}mutations.txt")
            if not mut_path.exists():
                logger.info("Genomic validation files not present — skipping (non-fatal)")
                return {}

            mut = pd.read_csv(mut_path, sep='\t')
            mut.columns = mut.columns.str.lower().str.strip()

            # ── Format 1: Standard cBioPortal mutation format ──────────────────
            hugo_col   = next((c for c in mut.columns if c in ['hugo_symbol', 'gene', 'symbol']), None)
            hgvs_col   = next((c for c in mut.columns if c in ['hgvsp_short', 'protein_change',
                                                                 'amino_acid_change', 'mutation']), None)
            sample_col = next((c for c in mut.columns if c in ['tumor_sample_barcode',
                                                                 'sample_id', 'sample']), None)

            if hugo_col and hgvs_col and sample_col:
                return self._parse_standard_format(mut, hugo_col, hgvs_col, sample_col)

            # ── Format 2: Transposed study format (genes as row index) ─────────
            # Detect: if column names look like gene symbols (short, uppercase-able)
            gene_like_cols = [
                c for c in mut.columns
                if c.upper() in {'H3-3A', 'H3F3A', 'DRD2', 'HDAC1', 'CDK4',
                                  'CDKN2A', 'EGFR', 'PTEN', 'EZH2', 'ACVR1'}
                or (len(c) <= 8 and c.replace('-', '').replace('_', '').isalnum())
            ]

            if gene_like_cols:
                logger.info(
                    "Detected transposed genomic format (gene columns: %s...). "
                    "Parsing as sample×gene matrix.",
                    ', '.join(gene_like_cols[:5])
                )
                return self._parse_transposed_format(mut_path)

            # ── Neither format matched ─────────────────────────────────────────
            logger.warning(
                "Genomic columns not recognised in mutations.txt.\n"
                "  Found columns: %s\n"
                "  Expected either:\n"
                "    Standard format: hugo_symbol, hgvsp_short, tumor_sample_barcode\n"
                "    Transposed format: gene names as row index (e.g. H3-3A, CDKN2A)\n"
                "  Genomic validation will be skipped. This is non-fatal — the pipeline\n"
                "  will run without CBTN co-occurrence statistics.",
                mut.columns.tolist()
            )
            return {}

        except Exception as e:
            logger.warning("Genomic validation error: %s", e)
            return {}

    def _parse_standard_format(self, mut, hugo_col, hgvs_col, sample_col) -> Dict:
        """Parse standard cBioPortal mutation format."""
        h3k27m_samples = set(mut[
            (mut[hugo_col].str.upper().isin(['H3-3A', 'H3F3A'])) &
            (mut[hgvs_col].str.contains('K28M|K27M', na=False, case=False))
        ][sample_col])

        return self._load_cna_and_rna(h3k27m_samples)

    def _parse_transposed_format(self, mut_path: Path) -> Dict:
        """
        Parse transposed format (PNOC/PBTA cBioPortal export) where:
          - First column = gene names (H3-3A, CDKN2A, etc.) → becomes row index
          - Remaining columns = sample IDs
          - Cell values = mutation calls (K28M, WT, NP, G35R, etc.)

        This is the format produced by cBioPortal "Transposed Matrix" download.
        Row 0 = STUDY_ID, Row 1 = SAMPLE_ID are metadata rows that get skipped
        when index_col=0 is used — pandas treats first column as index automatically.
        """
        # Re-read with index_col=0 so gene names are the row index
        raw = pd.read_csv(mut_path, sep='\t', index_col=0)

        # Normalize index to uppercase for robust matching
        raw.index = raw.index.astype(str).str.upper().str.strip()

        # Find the H3-3A row (gene encoding histone H3.3, mutated in DIPG)
        H3_ALIASES = {'H3-3A', 'H3F3A', 'HIST1H3B', 'H33A', 'H3.3A'}
        h3_row = next(
            (idx for idx in raw.index if idx in H3_ALIASES),
            None
        )

        if h3_row is None:
            logger.warning(
                "Transposed format detected but could not find H3K27M column or sample ID column. "
                "Skipping genomic validation.\n"
                "  Gene rows found: %s",
                raw.index.tolist()[:10]
            )
            return {}

        # Get mutation values for H3-3A across all samples
        h3_values = raw.loc[h3_row].astype(str).str.upper().str.strip()

        # K28M = H3K27M in updated HGVS nomenclature (1-indexed vs 0-indexed)
        # K27M = older nomenclature (also valid)
        # EXCLUDE G35R, G35V — these are H3.3 mutations but NOT H3K27M
        H3K27M_VALUES = {'K28M', 'K27M'}
        h3k27m_samples = set(raw.columns[h3_values.isin(H3K27M_VALUES)])

        logger.info(
            "Transposed format: found %d H3K27M (K28M) samples out of %d total",
            len(h3k27m_samples), len(raw.columns)
        )

        return self._load_cna_and_rna(h3k27m_samples)

    def _load_cna_and_rna(self, h3k27m_samples: Set) -> Dict:
        """Load CNA and RNA data and compute co-occurrence statistics."""
        try:
            cna_path = Path(f"{self.data_dir}cna.txt")
            rna_path = Path(f"{self.data_dir}rna_zscores.txt")

            cdkn2a_del_samples = set()
            upregulated_genes  = set()
            total_samples      = max(len(h3k27m_samples), 1)

            if cna_path.exists():
                cna = pd.read_csv(cna_path, sep='\t', index_col=0)
                # Normalize index
                cna.index = cna.index.astype(str).str.upper().str.strip()
                total_samples = len(cna.columns)

                cdkn2a_idx = next(
                    (idx for idx in cna.index if idx == 'CDKN2A'), None
                )
                if cdkn2a_idx:
                    # Deep deletion = -2, hemizygous = -1
                    # Use -2 for strict CDKN2A homozygous deletion
                    cdkn2a_row = pd.to_numeric(cna.loc[cdkn2a_idx], errors='coerce')
                    cdkn2a_del_samples = set(cna.columns[cdkn2a_row <= -1])
                    logger.info(
                        "CDKN2A: %d samples with deletion (score ≤ -1), "
                        "%d with deep deletion (score = -2)",
                        len(cdkn2a_del_samples),
                        len(cna.columns[cdkn2a_row == -2])
                    )

            if rna_path.exists():
                rna = pd.read_csv(rna_path, sep='\t', index_col=0)
                overlap_rna = list(h3k27m_samples.intersection(rna.columns))
                if overlap_rna:
                    high_exp = rna[overlap_rna].apply(
                        pd.to_numeric, errors='coerce'
                    ).mean(axis=1)
                    upregulated_genes = set(high_exp[high_exp > 2.0].index)
                    logger.info(
                        "RNA: %d genes upregulated (z-score > 2.0) in H3K27M samples",
                        len(upregulated_genes)
                    )

            overlap = h3k27m_samples.intersection(cdkn2a_del_samples)

            logger.info(
                "Co-occurrence: H3K27M=%d, CDKN2A-del=%d, overlap=%d, total=%d",
                len(h3k27m_samples), len(cdkn2a_del_samples),
                len(overlap), total_samples
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

        # Genomic validation (handles missing files gracefully)
        genomic_stats = self._genomic_validator.validate_triple_combo_cohort()

        # p-value (returns None/nan instead of 1.0 when data absent)
        p_val = self._stat_validator.calculate_cooccurrence_p_value(genomic_stats)

        upregulated = genomic_stats.get("upregulated_genes", set())

        # Escape Route Analysis
        for drug in candidates:
            targets     = drug.get("targets", [])
            escape_hits = []
            for t in targets:
                neighbors = self._ppi.get_neighbors(t)
                escape_hits.extend([n for n in neighbors if n in upregulated])
            drug["resistance_nodes"]     = list(set(escape_hits))
            drug["escape_bypass_score"]  = 1.0 if not escape_hits else 0.4

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

        # BBB filter — penalises known GBM clinical failures (e.g. Cilengitide)
        # to 10% of their score BEFORE sorting, so they never appear in top results
        candidates, _ = self._bbb_filter.filter_and_rank(candidates, apply_penalty=True)

        sorted_candidates = sorted(candidates, key=lambda x: x.get("score", 0), reverse=True)

        hypotheses = self._hyp_gen.generate(
            candidates      = sorted_candidates[:top_k],
            cmap_results    = [],
            synergy_combos  = [],
            differential_cmap = [],
            genomic_stats   = genomic_stats,
            p_value         = p_val,
        )

        return {
            "hypotheses": hypotheses,
            "stats": {
                "p_value":            p_val,
                "p_value_label":      self._stat_validator.format_p_value_for_report(p_val),
                "genomic_data_loaded": bool(genomic_stats),
            },
        }