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
    Compute escape bypass score using curated DIPG resistance pathway data
    combined with RNA-derived upregulated gene set.

    FIX v5.5: The curated bypass map is now a FALLBACK signal only.
    When RNA data is available (upregulated_genes is populated), it drives
    the bypass score directly. The curated map supplements it for targets
    not covered by the RNA data.

    Score interpretation:
        1.00 = no active bypass routes found (drug should work)
        0.40 = 4+ active bypass routes (high resistance risk)
    """
    if not drug_targets:
        return ESCAPE["empty_target_score"]

    constitutive = ESCAPE["constitutive_resistance_nodes"]
    active_bypass_nodes: Set[str] = set()
    covered_targets:     Set[str] = set()

    has_rna_data = len(upregulated_genes) >= ESCAPE["min_rna_genes_for_string_weight"]

    for target in drug_targets:
        target_upper      = target.upper()
        bypass_candidates = DIPG_RESISTANCE_BYPASS_MAP.get(target_upper, [])
        if bypass_candidates:
            covered_targets.add(target_upper)
            for node in bypass_candidates:
                node_upper = node.upper()
                # FIX: when RNA data is present, only count bypass nodes that
                # are CONFIRMED upregulated in actual patient data, plus
                # constitutive nodes. Without RNA data, all curated nodes count.
                if node_upper in constitutive:
                    active_bypass_nodes.add(node_upper)
                elif has_rna_data and node_upper in upregulated_genes:
                    active_bypass_nodes.add(node_upper)
                elif not has_rna_data:
                    # No RNA data — fall back to curated-only scoring
                    active_bypass_nodes.add(node_upper)

    if not covered_targets:
        return ESCAPE["no_target_score"]

    n_bypass      = len(active_bypass_nodes)
    bypass_scores = ESCAPE["bypass_scores"]
    max_key       = max(bypass_scores.keys())
    score         = bypass_scores.get(n_bypass, bypass_scores[max_key])

    if has_rna_data and n_bypass > 0:
        logger.debug(
            "Escape bypass (RNA-confirmed): %d active bypass nodes for targets %s: %s",
            n_bypass, drug_targets[:3], sorted(active_bypass_nodes)[:5]
        )

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
                    "Detected transposed genomic format. Parsing as sample×gene matrix.",
                )
                return self._parse_transposed_format(mut_path)

            logger.warning(
                "Genomic columns not recognised in mutations.txt. "
                "Found columns: %s. Genomic validation will be skipped.",
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

            # ── Also detect ACVR1-mutant samples for subgroup stratification ──
            acvr1_samples = set()

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

                metadata_rows     = GENOMICS["rna_metadata_rows"]
                min_genes         = GENOMICS["rna_min_genes_required"]
                h3k27m_indicators = GENOMICS["rna_h3k27m_col_indicators"]
                normal_indicators = GENOMICS["rna_normal_col_indicators"]

                gene_rows = [r for r in rna.index
                             if str(r).upper().strip() not in metadata_rows]

                if len(gene_rows) < min_genes:
                    logger.warning(
                        "⚠️  rna_zscores.txt contains only %d gene rows — "
                        "clinical summary file, not an expression matrix. "
                        "RNA upregulation skipped. Escape bypass will use "
                        "constitutive DIPG resistance nodes only "
                        "(curated fallback mode — less precise).",
                        len(gene_rows),
                    )
                else:
                    rna_expr    = rna.loc[gene_rows]
                    overlap_rna = list(h3k27m_samples.intersection(rna_expr.columns))

                    if overlap_rna:
                        high_exp          = rna_expr[overlap_rna].apply(
                            pd.to_numeric, errors="coerce"
                        ).mean(axis=1)
                        upregulated_genes = set(high_exp[high_exp > rna_threshold].index)
                        logger.info(
                            "RNA (matched cohort): %d genes upregulated "
                            "(value > %.1f) in %d H3K27M samples",
                            len(upregulated_genes), rna_threshold, len(overlap_rna),
                        )
                    else:
                        h3k27m_ref_cols = [
                            c for c in rna_expr.columns
                            if any(ind.lower() in c.lower() for ind in h3k27m_indicators)
                            and not any(n.lower() in c.lower() for n in normal_indicators)
                        ]
                        normal_ref_cols = [
                            c for c in rna_expr.columns
                            if any(n.lower() in c.lower() for n in normal_indicators)
                        ]

                        if h3k27m_ref_cols:
                            high_exp = rna_expr[h3k27m_ref_cols].apply(
                                pd.to_numeric, errors="coerce"
                            ).mean(axis=1)
                            upregulated_genes = set(
                                high_exp[high_exp > rna_threshold].index
                            )
                            logger.info(
                                "RNA (reference cohort, n=%d H3K27M samples): "
                                "%d genes upregulated. Normal excluded: %s. "
                                "NOTE: n=5 reference samples is small — "
                                "escape bypass scores carry higher uncertainty.",
                                len(h3k27m_ref_cols),
                                len(upregulated_genes),
                                normal_ref_cols,
                            )
                        else:
                            logger.warning(
                                "RNA: no sample overlap and no H3K27M indicator columns. "
                                "Escape bypass uses constitutive nodes only (curated fallback)."
                            )

            overlap = h3k27m_samples.intersection(cdkn2a_del_samples)
            logger.info(
                "Co-occurrence: H3K27M=%d, CDKN2A-del=%d, overlap=%d, total=%d",
                len(h3k27m_samples), len(cdkn2a_del_samples), len(overlap), total_samples,
            )

            # ── ACVR1 subgroup: DIRECT mutation calling ────────────────────────
            # Known activating ACVR1 mutations in DIPG (Taylor 2014 Nat Genet)
            # R206H (~50% of ACVR1-mutant), G328V, G328E, G356D, R258G
            acvr1_mutations = {"R206H", "G328V", "G328E", "G356D", "R258G",
                               "G328W", "R206C", "G325A"}
            acvr1_samples   = set()

            # Try to call from mutations.txt if columns are available
            mut_path2 = data_dir / "mutations.txt"
            if mut_path2.exists():
                try:
                    mut2 = pd.read_csv(mut_path2, sep="\t")
                    mut2.columns = mut2.columns.str.lower().str.strip()
                    hugo_col   = next((c for c in mut2.columns
                                       if c in ["hugo_symbol", "gene", "symbol"]), None)
                    hgvs_col   = next((c for c in mut2.columns
                                       if c in ["hgvsp_short", "protein_change",
                                                "amino_acid_change", "mutation"]), None)
                    sample_col = next((c for c in mut2.columns
                                       if c in ["tumor_sample_barcode", "sample_id", "sample"]), None)

                    if hugo_col and hgvs_col and sample_col:
                        acvr1_rows = mut2[mut2[hugo_col].str.upper().str.strip() == "ACVR1"]
                        for _, row in acvr1_rows.iterrows():
                            aa_change = str(row[hgvs_col]).upper().strip()
                            for mut in acvr1_mutations:
                                if mut in aa_change:
                                    acvr1_samples.add(row[sample_col])
                                    break
                        logger.info(
                            "ACVR1 direct mutation calling: %d samples with "
                            "known activating mutations (R206H/G328V/G328E/G356D/R258G). "
                            "Prevalence: %.1f%% of H3K27M+ samples.",
                            len(acvr1_samples),
                            len(acvr1_samples) / max(len(h3k27m_samples), 1) * 100,
                        )
                    else:
                        # Transposed format — scan for ACVR1 row
                        raw2 = pd.read_csv(mut_path2, sep="\t", index_col=0)
                        raw2.index = raw2.index.astype(str).str.upper().str.strip()
                        if "ACVR1" in raw2.index:
                            acvr1_vals = raw2.loc["ACVR1"].astype(str).str.upper().str.strip()
                            for col in raw2.columns:
                                for mut in acvr1_mutations:
                                    if mut in acvr1_vals[col]:
                                        acvr1_samples.add(col)
                                        break
                            logger.info(
                                "ACVR1 (transposed format): %d samples called directly",
                                len(acvr1_samples),
                            )
                except Exception as e:
                    logger.debug("ACVR1 direct calling failed: %s", e)

            # Fall back to prevalence prior if direct calling found nothing
            if not acvr1_samples:
                acvr1_estimated_n = round(len(h3k27m_samples) * 0.25)
                acvr1_calling_method = "prevalence_prior_25pct"
                logger.info(
                    "ACVR1: direct calling found 0 samples — using 25%% prevalence prior "
                    "(estimated n=%d). Check mutations.txt for ACVR1 amino acid changes.",
                    acvr1_estimated_n,
                )
            else:
                acvr1_estimated_n    = len(acvr1_samples)
                acvr1_calling_method = "direct_mutation_calling"

            return {
                "h3k27m_count":         len(h3k27m_samples),
                "cdkn2a_del_count":     len(cdkn2a_del_samples),
                "overlap_count":        len(overlap),
                "prevalence":           len(overlap) / max(len(h3k27m_samples), 1),
                "upregulated_genes":    upregulated_genes,
                "total_samples":        total_samples,
                "rna_gene_count":       len(upregulated_genes),
                "has_rna_data":         len(upregulated_genes) >= ESCAPE["min_rna_genes_for_string_weight"],
                "acvr1_estimated_n":    acvr1_estimated_n,
                "acvr1_samples":        acvr1_samples,
                "acvr1_calling_method": acvr1_calling_method,
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
        has_rna_data  = genomic_stats.get("has_rna_data", False)

        logger.info(
            "RNA data status: %s (%d upregulated genes). "
            "Escape bypass mode: %s",
            "loaded" if has_rna_data else "not loaded / too sparse",
            len(upregulated),
            "RNA-confirmed" if has_rna_data else "curated fallback (less precise)",
        )

        # ── Escape Route Analysis ─────────────────────────────────────────────
        ew = ESCAPE
        for drug in candidates:
            targets = drug.get("targets", [])

            live_escape_hits = []
            for t in targets:
                neighbors = self._ppi.get_neighbors(t)
                live_escape_hits.extend([n for n in neighbors if n in upregulated])

            drug["resistance_nodes"] = list(set(live_escape_hits))

            # FIX: compute_escape_bypass_score now uses RNA-confirmed nodes
            # when has_rna_data=True, curated-only when False
            curated_score = compute_escape_bypass_score(targets, upregulated)

            if live_escape_hits and has_rna_data:
                # STRING-DB + RNA-confirmed escape hits
                drug["escape_bypass_score"] = round(
                    ew["string_weight"]  * ew["string_hit_score"]
                    + ew["curated_weight"] * curated_score,
                    4
                )
                drug["escape_note"] = (
                    f"STRING-DB + RNA-confirmed DIPG bypass "
                    f"({len(live_escape_hits)} live hits)"
                )
            elif live_escape_hits and not has_rna_data:
                # STRING-DB hits but no RNA confirmation — use equal weighting
                drug["escape_bypass_score"] = round(
                    0.50 * ew["string_hit_score"]
                    + 0.50 * curated_score,
                    4
                )
                drug["escape_note"] = (
                    "STRING-DB bypass (no RNA confirmation — curated fallback)"
                )
            else:
                drug["escape_bypass_score"] = round(curated_score, 4)
                drug["escape_note"]         = (
                    "Curated DIPG resistance bypass"
                    if not has_rna_data
                    else "RNA-informed curated bypass (no STRING-DB hits)"
                )

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

        # DIPG-specific BBB penalty
        # Brainstem BBB is tighter than cortex. LOW BBB drugs are non-viable
        # for DIPG regardless of target essentiality. Apply score penalty to
        # push them below CNS-penetrant alternatives in the ranking.
        # LOW=0.50x  (effectively excluded from top candidates)
        # UNKNOWN=0.85x (mild penalty — lack of data, not confirmed poor)
        # Sources: Pardridge 2003; Fischer 1998; DIPG brainstem PK literature
        dipg_bbb_penalties = {"LOW": 0.50, "UNKNOWN": 0.85}
        n_penalised = 0
        for c in candidates:
            bbb = c.get("bbb_penetrance", "UNKNOWN")
            if bbb in dipg_bbb_penalties:
                c["score_before_dipg_bbb"] = c["score"]
                c["score"] = round(c["score"] * dipg_bbb_penalties[bbb], 4)
                c["dipg_bbb_penalty"] = dipg_bbb_penalties[bbb]
                n_penalised += 1
        if n_penalised:
            logger.info(
                "DIPG BBB penalty applied: %d candidates penalised "
                "(LOW x0.50, UNKNOWN x0.85) — brainstem penetrance requirement",
                n_penalised,
            )

        sorted_candidates = sorted(candidates, key=lambda x: x.get("score", 0), reverse=True)

        # IC50 validation — annotate top candidates with published cell-line data
        # Non-fatal: if module missing or data absent, candidates are unchanged
        try:
            from .published_ic50_validation import annotate_candidates_with_ic50
            sorted_candidates = annotate_candidates_with_ic50(sorted_candidates)
        except Exception as e:
            logger.debug("IC50 validation skipped: %s", e)

        # Generate hypotheses
        hypotheses = self._hyp_gen.generate(
            candidates        = sorted_candidates[:top_k],
            cmap_results      = [],
            synergy_combos    = [],
            differential_cmap = [],
            genomic_stats     = genomic_stats,
            p_value           = p_val,
        )

        # Extract toxicity stats from hypothesis
        tox_breakdown = {}
        if hypotheses and "confidence_breakdown" in hypotheses[0]:
            bd = hypotheses[0]["confidence_breakdown"]
            tox_breakdown = {
                "toxicity_flag":              bd.get("toxicity_flag", "UNKNOWN"),
                "toxicity_multiplier":        bd.get("toxicity_multiplier", 1.0),
                "toxicity_multiplier_optimistic": bd.get("toxicity_multiplier_optimistic", 1.0),
                "toxicity_note":              bd.get("toxicity_note", ""),
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
                "escape_bypass_mode":  "RNA-confirmed" if has_rna_data else "curated fallback",
                "rna_upregulated_genes": len(upregulated),
                "acvr1_estimated_n":   genomic_stats.get("acvr1_estimated_n", 0),
                **tox_breakdown,
            },
        }