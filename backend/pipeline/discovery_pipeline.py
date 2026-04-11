"""
discovery_pipeline.py  (v1.1 — generic drug filter + import fixes)
====================================================================
Changes from v1.0:
  - Added generic_only parameter to run() — filters to drugs with
    FDA-approved generic formulations (more accessible for paediatric use)
  - Fixed _ATRT_FALLBACK_CANDIDATES to include generic status flags
  - Fixed import order to avoid circular dependencies
"""

import asyncio
import logging
import pandas as pd
from typing import Dict, List, Optional, Set
from pathlib import Path

try:
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
        COMPOSITE_WEIGHTS, SCORE_DEFAULTS, ESCAPE, GENOMICS, PATHS,
        EZH2_INHIBITOR, AURKA_INHIBITOR, SMO_INHIBITOR, OPENTARGETS,
        BBB as BBB_CONFIG,
    )
except ImportError:
    from data_fetcher import ProductionDataFetcher
    from cmap_query import CMAPQuery
    from ppi_network import PPINetwork
    from depmap_essentiality import DepMapEssentiality
    from synergy_predictor import SynergyPredictor
    from hypothesis_generator import HypothesisGenerator
    from statistical_validator import StatisticalValidator
    from tissue_expression import TissueExpressionScorer
    from bbb_filter import BBBFilter
    from pipeline_config import (
        COMPOSITE_WEIGHTS, SCORE_DEFAULTS, ESCAPE, GENOMICS, PATHS,
        EZH2_INHIBITOR, AURKA_INHIBITOR, SMO_INHIBITOR, OPENTARGETS,
        BBB as BBB_CONFIG,
    )

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Generic drug flag — drugs confirmed to have FDA-approved generics (April 2026)
# Source: FDA Orange Book; verified generic availability
# Used when generic_only=True is passed to run()
# ─────────────────────────────────────────────────────────────────────────────

KNOWN_GENERIC_DRUGS: set = {
    "valproic acid", "metformin", "metformin hcl",
    "chloroquine", "hydroxychloroquine", "sirolimus",
    "bortezomib",    # US patent expired 2022
    "imatinib",      # US patent expired 2016
    "temozolomide",  # patent expired 2014
    "dexamethasone", "lomustine", "carmustine", "itraconazole",
    "tretinoin",     # ATRA (all-trans retinoic acid)
    "arsenic trioxide", "vorinostat",
}

# Drugs with generic status that are particularly relevant to ATRT
ATRT_GENERIC_RELEVANT: set = {
    "valproic acid",      # HDAC inhibitor — fully generic, low cost
    "sirolimus",          # mTOR inhibitor — generic rapamycin
    "chloroquine",        # autophagy inhibitor — generic
    "hydroxychloroquine", # autophagy inhibitor — generic
    "itraconazole",       # SMO inhibitor (repurposed antifungal) — generic
    "arsenic trioxide",   # GLI inhibitor (repurposed APL drug) — generic
    "metformin",          # AMPK activator — generic, pennies/pill
    "vorinostat",         # HDAC inhibitor — generic filed
    "temozolomide",       # alkylating — generic
    "dexamethasone",      # steroid — generic (use sparingly — worsens immunity)
}


def _is_generic(drug: Dict) -> bool:
    """Return True if drug has a verified generic formulation."""
    # First check pipeline annotation from GenericDrugFetcher
    if drug.get("has_generic"):
        return True
    # Fallback: check known generic set
    name = (drug.get("name") or drug.get("drug_name") or "").lower().strip()
    return name in KNOWN_GENERIC_DRUGS


# ─────────────────────────────────────────────────────────────────────────────
# ATRT Resistance Bypass Network
# ─────────────────────────────────────────────────────────────────────────────

ATRT_RESISTANCE_BYPASS_MAP: Dict[str, List[str]] = {
    "EZH2":   ["BRD4", "CDK4", "MYC", "MYCN"],
    "EED":    ["BRD4", "CDK4", "MYC"],
    "SUZ12":  ["BRD4", "MYC"],
    "BRD4":   ["CDK4", "MYC", "MYCN", "MTOR"],
    "BRD2":   ["MYC", "CDK4", "EZH2"],
    "BRD3":   ["MYC", "EZH2"],
    "HDAC1":  ["BCL2", "BCL2L1", "MCL1", "MYC"],
    "HDAC2":  ["BCL2", "MCL1", "MYC"],
    "HDAC3":  ["BCL2", "NFKB1", "MYC"],
    "CDK4":   ["PIK3CA", "MTOR", "AKT1", "CCNE1"],
    "CDK6":   ["PIK3CA", "MTOR", "AKT1"],
    "AURKA":  ["CDK4", "MYC", "BCL2L1", "MYCN"],
    "AURKB":  ["CDK4", "MYC"],
    "MTOR":   ["PIK3CA", "AKT1", "EIF4E"],
    "PIK3CA": ["MTOR", "AKT1", "MAPK1"],
    "PSMB5":  ["ATG5", "BECN1", "SQSTM1", "MCL1"],
    "PSMB2":  ["ATG5", "BECN1", "MCL1"],
    "MYC":    ["AURKA", "CDK4", "BCL2L1", "BRD4"],
    "MYCN":   ["AURKA", "CDK4", "BRD4"],
    "BCL2":   ["MCL1", "BCL2L1"],
    "BCL2L1": ["MCL1", "BCL2"],
    "MCL1":   ["BCL2", "BCL2L1"],
    "SMO":    ["GLI2", "PIK3CA", "MTOR"],
    "GLI2":   ["MYC", "CDK4", "PIK3CA"],
    "PARP1":  ["RAD51", "BRCA1", "ATM"],
}

ATRT_CONSTITUTIVE_RESISTANCE: Set[str] = {
    "MYC", "MYCN", "BCL2L1", "MCL1", "CDK4", "BRD4", "PIK3CA", "EZH2",
}


def compute_escape_bypass_score(
    drug_targets: List[str],
    upregulated_genes: Set[str],
) -> float:
    if not drug_targets:
        return ESCAPE["empty_target_score"]

    constitutive  = ESCAPE["constitutive_resistance_nodes"]
    active_bypass: Set[str] = set()
    covered:       Set[str] = set()
    has_rna_data = len(upregulated_genes) >= ESCAPE["min_rna_genes_for_string_weight"]

    for target in drug_targets:
        target_upper = target.upper()
        bypass_candidates = ATRT_RESISTANCE_BYPASS_MAP.get(target_upper, [])
        if bypass_candidates:
            covered.add(target_upper)
            for node in bypass_candidates:
                node_upper = node.upper()
                if node_upper in constitutive:
                    active_bypass.add(node_upper)
                elif has_rna_data and node_upper in upregulated_genes:
                    active_bypass.add(node_upper)
                elif not has_rna_data:
                    active_bypass.add(node_upper)

    if not covered:
        return ESCAPE["no_target_score"]

    n_bypass = len(active_bypass)
    scores   = ESCAPE["bypass_scores"]
    max_key  = max(scores.keys())
    return scores.get(n_bypass, scores[max_key])


# ─────────────────────────────────────────────────────────────────────────────
# ATRT Genomic Validator
# ─────────────────────────────────────────────────────────────────────────────

class ATRTGenomicValidator:
    def __init__(self, data_dir: str = None):
        self.data_dir = Path(data_dir or PATHS["genomics"])

    def validate_atrt_cohort(self) -> Dict:
        result = {
            "smarcb1_loss_count":  0,
            "smarca4_loss_count":  0,
            "upregulated_genes":   set(),
            "downregulated_genes": set(),
            "total_samples":       0,
            "has_rna_data":        False,
            "subgroup_counts":     {"TYR": 0, "SHH": 0, "MYC": 0},
            "calling_method":      "none",
        }

        cna_path = self.data_dir / "cna.txt"
        rna_path = self.data_dir / "rna_zscores.txt"

        if not cna_path.exists() and not rna_path.exists():
            logger.info(
                "CBTN ATRT files not found at %s — using GSE70678 fallback.",
                self.data_dir,
            )
            return self._load_gse70678_fallback(result)

        if cna_path.exists():
            result = self._load_cna(cna_path, result)
        if rna_path.exists():
            result = self._load_rna(rna_path, result)

        return result

    def _load_cna(self, cna_path: Path, result: Dict) -> Dict:
        try:
            cna = pd.read_csv(cna_path, sep="\t", index_col=0)
            cna.index = cna.index.astype(str).str.upper().str.strip()
            result["total_samples"] = len(cna.columns)
            threshold = GENOMICS["smarcb1_del_threshold"]
            for alias in GENOMICS["smarcb1_aliases"]:
                if alias in cna.index:
                    row = pd.to_numeric(cna.loc[alias], errors="coerce")
                    result["smarcb1_loss_count"] = int((row <= threshold).sum())
                    logger.info("SMARCB1 loss: %d/%d samples",
                                result["smarcb1_loss_count"], result["total_samples"])
                    break
            for alias in GENOMICS["smarca4_aliases"]:
                if alias in cna.index:
                    row = pd.to_numeric(cna.loc[alias], errors="coerce")
                    result["smarca4_loss_count"] = int((row <= threshold).sum())
                    break
        except Exception as e:
            logger.warning("CNA loading failed: %s", e)
        return result

    def _load_rna(self, rna_path: Path, result: Dict) -> Dict:
        try:
            rna = pd.read_csv(rna_path, sep="\t", index_col=0)
            threshold = GENOMICS["rna_upregulation_zscore"]
            down_threshold = GENOMICS["rna_downregulation_zscore"]
            atrt_indicators = GENOMICS["rna_h3k27m_col_indicators"]
            normal_indicators = GENOMICS["rna_normal_col_indicators"]
            metadata_rows = GENOMICS["rna_metadata_rows"]

            gene_rows = [r for r in rna.index
                         if str(r).upper().strip() not in metadata_rows]
            if len(gene_rows) < GENOMICS["rna_min_genes_required"]:
                logger.warning("Too few gene rows in rna_zscores.txt (%d)", len(gene_rows))
                return result

            rna_expr = rna.loc[gene_rows]
            atrt_cols = [c for c in rna_expr.columns
                         if any(ind.lower() in c.lower() for ind in atrt_indicators)
                         and not any(n.lower() in c.lower() for n in normal_indicators)]

            if not atrt_cols:
                return result

            expr_numeric = rna_expr[atrt_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
            result["upregulated_genes"]   = set(expr_numeric[expr_numeric > threshold].index)
            result["downregulated_genes"] = set(expr_numeric[expr_numeric < down_threshold].index)
            result["has_rna_data"]        = (
                len(result["upregulated_genes"]) >= ESCAPE["min_rna_genes_for_string_weight"]
            )
            result["calling_method"] = "cbtn_rna"
        except Exception as e:
            logger.warning("RNA loading failed: %s", e)
        return result

    def _load_gse70678_fallback(self, result: Dict) -> Dict:
        gse_path = Path(PATHS["scrna"])
        if not gse_path.exists():
            result["subgroup_counts"] = {"TYR": 54, "SHH": 55, "MYC": 41}
            result["calling_method"] = "prevalence_prior"
            return result

        try:
            df = pd.read_csv(gse_path, sep="\t", index_col=0)
            df.index = df.index.astype(str).str.upper().str.strip()
            atrt_indicators   = GENOMICS["rna_h3k27m_col_indicators"]
            normal_indicators = GENOMICS["rna_normal_col_indicators"]

            atrt_cols = [c for c in df.columns
                         if any(ind.lower() in c.lower() for ind in atrt_indicators)
                         and not any(n.lower() in c.lower() for n in normal_indicators)]
            normal_cols = [c for c in df.columns
                           if any(n.lower() in c.lower() for n in normal_indicators)]

            if not atrt_cols:
                logger.warning("No ATRT columns in GSE70678.")
                return result

            result["total_samples"]      = len(atrt_cols)
            result["smarcb1_loss_count"] = len(atrt_cols)
            result["calling_method"]     = "gse70678_fallback"

            df_num    = df.apply(pd.to_numeric, errors="coerce")
            atrt_mean = df_num[atrt_cols].mean(axis=1)

            if normal_cols:
                diff = atrt_mean - df_num[normal_cols].mean(axis=1)
            else:
                diff = atrt_mean

            threshold = GENOMICS["rna_upregulation_zscore"]
            down_threshold = GENOMICS["rna_downregulation_zscore"]
            result["upregulated_genes"]   = set(diff[diff > threshold].index)
            result["downregulated_genes"] = set(diff[diff < down_threshold].index)
            result["has_rna_data"]        = (
                len(result["upregulated_genes"]) >= ESCAPE["min_rna_genes_for_string_weight"]
            )
            logger.info(
                "GSE70678 fallback: %d ATRT samples, %d upregulated genes",
                len(atrt_cols), len(result["upregulated_genes"]),
            )
        except Exception as e:
            logger.warning("GSE70678 fallback failed: %s", e)

        return result

    def calculate_smarcb1_statistics(self, genomic_stats: Dict) -> Dict:
        total   = genomic_stats.get("total_samples", 0)
        smarcb1 = genomic_stats.get("smarcb1_loss_count", 0)
        smarca4 = genomic_stats.get("smarca4_loss_count", 0)
        return {
            "smarcb1_loss_count":   smarcb1,
            "smarca4_loss_count":   smarca4,
            "total_loss_count":     smarcb1 + smarca4,
            "total_samples":        total,
            "smarcb1_prevalence":   smarcb1 / max(total, 1),
            "combined_prevalence":  (smarcb1 + smarca4) / max(total, 1),
            "p_value":              None,
            "p_value_label": (
                "N/A — SMARCB1 biallelic loss is the defining event in ATRT. "
                "(Hasselblatt 2011 Acta Neuropathol PMID 20625942)"
            ),
            "statistical_note": (
                f"SMARCB1 biallelic loss: {smarcb1}/{total} samples "
                f"({smarcb1 / max(total, 1):.0%}). "
                "All SMARCB1-null samples are candidates for EZH2 synthetic lethality "
                "(Knutson 2013 PNAS 110(19):7922, PMID 23620515)."
            ),
        }


# ─────────────────────────────────────────────────────────────────────────────
# ATRT Drug Score Boosts
# ─────────────────────────────────────────────────────────────────────────────

def _apply_atrt_drug_boosts(
    candidates: List[Dict],
    subgroup: Optional[str] = None,
) -> List[Dict]:
    n_ezh2 = n_aurka = n_smo = 0

    for c in candidates:
        name_lower = (c.get("name") or c.get("drug_name") or "").lower().strip()
        mech_lower = (c.get("mechanism") or "").lower()
        targets    = [t.upper() for t in (c.get("targets") or [])]

        for suffix in (" hydrochloride", " hcl", " sodium", " mesylate", " phosphate"):
            name_lower = name_lower.replace(suffix, "")

        base_score = c.get("score", 0.0)
        c["atrt_boosts_applied"] = []

        # EZH2 inhibitor BOOST (not penalty — SMARCB1 synthetic lethality)
        is_ezh2 = (
            name_lower in EZH2_INHIBITOR["known_inhibitors"]
            or any(kw in mech_lower for kw in EZH2_INHIBITOR["mechanism_keywords"])
            or (targets == ["EZH2"])
        )
        if is_ezh2:
            mult = EZH2_INHIBITOR["composite_boost"]
            c["score"] = round(min(1.0, base_score * mult), 4)
            c["atrt_boosts_applied"].append(
                f"EZH2-inhibitor ×{mult} (SMARCB1 synthetic lethality — Knutson 2013 PNAS)"
            )
            c["ezh2_boosted"] = True
            n_ezh2 += 1
            continue

        c["ezh2_boosted"] = False

        # AURKA inhibitor BOOST
        is_aurka = (
            name_lower in AURKA_INHIBITOR["known_inhibitors"]
            or any(kw in mech_lower for kw in AURKA_INHIBITOR["mechanism_keywords"])
            or "AURKA" in targets
        )
        if is_aurka:
            mult = (
                AURKA_INHIBITOR["composite_boost_myc_subgroup"]
                if subgroup and subgroup.upper() == "MYC"
                else AURKA_INHIBITOR["composite_boost_other"]
            )
            c["score"] = round(min(1.0, base_score * mult), 4)
            c["atrt_boosts_applied"].append(
                f"AURKA-inhibitor ×{mult} (MYCN stabilisation — Sredni 2017)"
            )
            c["aurka_boosted"] = True
            n_aurka += 1
        else:
            c["aurka_boosted"] = False

        # SMO/GLI inhibitor BOOST (SHH subgroup)
        is_smo = (
            name_lower in SMO_INHIBITOR["known_inhibitors"]
            or any(kw in mech_lower for kw in SMO_INHIBITOR["mechanism_keywords"])
            or bool({"SMO", "GLI2", "GLI1", "PTCH1"} & set(targets))
        )
        if is_smo:
            mult = (
                SMO_INHIBITOR["composite_boost_shh_subgroup"]
                if subgroup and subgroup.upper() == "SHH"
                else SMO_INHIBITOR["composite_boost_other"]
            )
            if mult > 1.0:
                c["score"] = round(min(1.0, c.get("score", base_score) * mult), 4)
                c["atrt_boosts_applied"].append(
                    f"SMO-inhibitor ×{mult} (SHH subgroup — Torchia 2015)"
                )
                n_smo += 1

    logger.info("ATRT drug boosts: EZH2×%d | AURKA×%d | SMO×%d", n_ezh2, n_aurka, n_smo)
    return candidates


# ─────────────────────────────────────────────────────────────────────────────
# Production Pipeline
# ─────────────────────────────────────────────────────────────────────────────

class ProductionPipeline:
    """
    ATRT Drug Repurposing Pipeline — Production Orchestrator.

    Parameters
    ----------
    generic_only : bool  (can also be set per-run in run())
        If True, filters candidates to those with FDA-approved generic
        formulations. More accessible/affordable for paediatric use.
        Generic drugs in ATRT context: valproic acid (HDAC-i), sirolimus
        (mTOR-i), itraconazole (SMO-i), arsenic trioxide (GLI-i), metformin.
    """

    def __init__(self, generic_only: bool = False):
        self.generic_only       = generic_only
        self._data_fetcher      = ProductionDataFetcher()
        self._cmap              = CMAPQuery()
        self._ppi               = PPINetwork()
        self._depmap            = DepMapEssentiality()
        self._tissue            = TissueExpressionScorer("atrt")
        self._synergy           = SynergyPredictor()
        self._hyp_gen           = HypothesisGenerator()
        self._genomic_validator = ATRTGenomicValidator()
        self._stat_validator    = StatisticalValidator()
        self._bbb_filter        = BBBFilter()

    async def initialize(self, disease: str = "atrt"):
        return True

    async def analyze_disease(
        self,
        disease_name: str = "atrt",
        min_score: float = 0.2,
        max_results: int = 20,
        generic_only: bool = False,
    ) -> Dict:
        result = await self.run(
            disease_name=disease_name,
            top_k=max_results,
            generic_only=generic_only or self.generic_only,
        )
        candidates = result.get("top_candidates", [])
        for c in candidates:
            if "indication" not in c:
                c["indication"] = "ATRT / rhabdoid tumor"
            if "mechanism" not in c:
                c["mechanism"] = ""
            c["drug_name"]       = c.get("name", "")
            c["composite_score"] = c.get("score", 0.0)
        return {
            "success":    True,
            "candidates": [c for c in candidates if c.get("score", 0) >= min_score],
            "stats":      result.get("stats", {}),
        }

    async def run(
        self,
        disease_name:  str  = "atrt",
        top_k:         int  = 20,
        subgroup:      Optional[str] = None,
        location:      str  = "unknown_location",
        generic_only:  bool = False,
    ) -> Dict:
        """
        Run the full ATRT pipeline.

        Parameters
        ----------
        generic_only : bool
            Filter output to drugs with FDA-approved generic formulations.
            Useful for resource-limited settings or cost-sensitive analysis.
        """
        if subgroup:
            self._tissue.set_subgroup(subgroup)

        effective_generic = generic_only or self.generic_only

        # ── Step 1: Fetch candidates ──────────────────────────────────────────
        disease_data = await self._data_fetcher.fetch_disease_data(disease_name)
        candidates   = await self._data_fetcher.fetch_approved_drugs(
            annotate_generics=True
        )

        if not candidates:
            logger.warning("OpenTargets returned 0 drugs. Using ATRT fallback library.")
            candidates = _ATRT_FALLBACK_CANDIDATES()

        # ── Generic filter (applied early to reduce compute) ──────────────────
        if effective_generic:
            original_n  = len(candidates)
            candidates  = [c for c in candidates if _is_generic(c)]
            logger.info(
                "Generic filter: %d → %d candidates (kept drugs with FDA generics).",
                original_n, len(candidates),
            )
            if len(candidates) < 5:
                logger.warning(
                    "Only %d generic candidates found. "
                    "Adding ATRT-relevant generic drugs from curated list.",
                    len(candidates),
                )
                candidates = _ATRT_GENERIC_FALLBACK_CANDIDATES()

        # ── Step 2: Genomic data ──────────────────────────────────────────────
        genomic_stats = self._genomic_validator.validate_atrt_cohort()
        smarcb1_stats = self._genomic_validator.calculate_smarcb1_statistics(genomic_stats)
        upregulated   = genomic_stats.get("upregulated_genes", set())
        has_rna_data  = genomic_stats.get("has_rna_data", False)

        logger.info(
            "ATRT genomics: %s | %d upregulated genes | SMARCB1 loss: %d | "
            "escape bypass: %s",
            genomic_stats.get("calling_method", "none"),
            len(upregulated),
            genomic_stats.get("smarcb1_loss_count", 0),
            "RNA-confirmed" if has_rna_data else "curated fallback",
        )

        # ── Step 3: Escape bypass scoring ────────────────────────────────────
        ew = ESCAPE
        for drug in candidates:
            targets = drug.get("targets", [])
            live_escape_hits = []
            for t in targets:
                neighbors = self._ppi.get_neighbors(t)
                live_escape_hits.extend([n for n in neighbors if n in upregulated])
            drug["resistance_nodes"] = list(set(live_escape_hits))

            curated_score = compute_escape_bypass_score(targets, upregulated)
            if live_escape_hits and has_rna_data:
                drug["escape_bypass_score"] = round(
                    ew["string_weight"]  * ew["string_hit_score"]
                    + ew["curated_weight"] * curated_score, 4,
                )
                drug["escape_note"] = (
                    f"STRING-DB + RNA-confirmed bypass ({len(live_escape_hits)} hits)"
                )
            elif live_escape_hits:
                drug["escape_bypass_score"] = round(
                    0.50 * ew["string_hit_score"] + 0.50 * curated_score, 4
                )
                drug["escape_note"] = "STRING-DB bypass (no RNA confirmation)"
            else:
                drug["escape_bypass_score"] = round(curated_score, 4)
                drug["escape_note"] = (
                    "RNA-informed curated bypass" if has_rna_data
                    else "Curated ATRT resistance bypass"
                )

        # ── Step 4: Multi-omic scoring ────────────────────────────────────────
        candidates = await self._tissue.score_batch(candidates)
        candidates = await self._depmap.score_batch(candidates, disease_name)
        candidates = await self._ppi.score_batch(candidates, disease_data["genes"])

        # ── Step 5: Composite score ───────────────────────────────────────────
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
                4,
            )

        # ── Step 6: ATRT drug boosts ──────────────────────────────────────────
        candidates = _apply_atrt_drug_boosts(candidates, subgroup=subgroup)

        # ── Step 7: BBB filter + location penalty ─────────────────────────────
        candidates, _ = self._bbb_filter.filter_and_rank(candidates, apply_penalty=True)
        bbb_penalties = BBB_CONFIG["dipg_bbb_penalties"].get(
            location, BBB_CONFIG["dipg_bbb_penalties"]["unknown_location"]
        )
        for c in candidates:
            bbb = c.get("bbb_penetrance", "UNKNOWN")
            if bbb in bbb_penalties:
                c["score_before_bbb"] = c["score"]
                c["score"] = round(c["score"] * bbb_penalties[bbb], 4)
                c["bbb_penalty"] = bbb_penalties[bbb]

        # ── Step 8: Sort + IC50 annotation ───────────────────────────────────
        sorted_candidates = sorted(
            candidates, key=lambda x: x.get("score", 0), reverse=True
        )

        try:
            from .published_ic50_atrt_validation import annotate_atrt_candidates_with_ic50
            sorted_candidates = annotate_atrt_candidates_with_ic50(sorted_candidates)
        except ImportError:
            try:
                from published_ic50_atrt_validation import annotate_atrt_candidates_with_ic50
                sorted_candidates = annotate_atrt_candidates_with_ic50(sorted_candidates)
            except Exception as e:
                logger.debug("IC50 annotation skipped: %s", e)

        # ── Step 9: Hypothesis generation ────────────────────────────────────
        hypotheses = self._hyp_gen.generate(
            candidates     = sorted_candidates[:top_k],
            cmap_results   = [],
            synergy_combos = [],
            differential_cmap = [],
            genomic_stats  = {
                **smarcb1_stats,
                "has_rna_data":    has_rna_data,
                "subgroup_counts": genomic_stats.get("subgroup_counts", {}),
                "total_samples":   genomic_stats.get("total_samples", 0),
            },
            p_value = smarcb1_stats.get("p_value"),
        )

        # ── Step 10: CMap integration ────────────────────────────────────────
        if self._cmap.has_precomputed_scores():
            for c in sorted_candidates:
                name = c.get("name", "")
                cmap_score = self._cmap.get_precomputed_score(name)
                if cmap_score is not None:
                    c["cmap_score"] = cmap_score
                    c["is_reverser"] = cmap_score > 0.75

        return {
            "hypotheses":     hypotheses,
            "top_candidates": sorted_candidates[:top_k],
            "stats": {
                "smarcb1_loss_count":   genomic_stats.get("smarcb1_loss_count", 0),
                "smarca4_loss_count":   genomic_stats.get("smarca4_loss_count", 0),
                "total_samples":        genomic_stats.get("total_samples", 0),
                "upregulated_genes":    len(upregulated),
                "p_value":              smarcb1_stats.get("p_value"),
                "p_value_label":        smarcb1_stats.get("p_value_label", "N/A"),
                "statistical_note":     smarcb1_stats.get("statistical_note", ""),
                "n_screened":           len(sorted_candidates),
                "subgroup":             subgroup or "pan-ATRT",
                "location":             location,
                "composite_weights":    cw,
                "escape_bypass_mode":   "RNA-confirmed" if has_rna_data else "curated fallback",
                "rna_upregulated_genes": len(upregulated),
                "n_ezh2_boosted":       sum(1 for c in sorted_candidates if c.get("ezh2_boosted")),
                "n_aurka_boosted":      sum(1 for c in sorted_candidates if c.get("aurka_boosted")),
                "calling_method":       genomic_stats.get("calling_method", "none"),
                "generic_only":         effective_generic,
            },
        }


def _ATRT_FALLBACK_CANDIDATES() -> List[Dict]:
    """Minimal fallback when OpenTargets API is unavailable."""
    return [
        {"name": "Tazemetostat",  "targets": ["EZH2"], "has_generic": False, "cost_category": "BRAND"},
        {"name": "Panobinostat",  "targets": ["HDAC1", "HDAC2", "HDAC3", "HDAC6"], "has_generic": False, "cost_category": "BRAND"},
        {"name": "Alisertib",     "targets": ["AURKA"], "has_generic": False, "cost_category": "INVESTIGATIONAL"},
        {"name": "Birabresib",    "targets": ["BRD4", "BRD2", "BRD3"], "has_generic": False, "cost_category": "INVESTIGATIONAL"},
        {"name": "Abemaciclib",   "targets": ["CDK4", "CDK6"], "has_generic": False, "cost_category": "BRAND"},
        {"name": "Marizomib",     "targets": ["PSMB5", "PSMB2", "PSMB1"], "has_generic": False, "cost_category": "INVESTIGATIONAL"},
        {"name": "ONC201",        "targets": ["DRD2", "CLPB"], "has_generic": False, "cost_category": "BRAND"},
        {"name": "Paxalisib",     "targets": ["PIK3CA", "PIK3CD", "PIK3CG", "MTOR"], "has_generic": False, "cost_category": "INVESTIGATIONAL"},
        {"name": "Vismodegib",    "targets": ["SMO"], "has_generic": False, "cost_category": "BRAND"},
        {"name": "Vorinostat",    "targets": ["HDAC1", "HDAC2", "HDAC3"], "has_generic": False, "cost_category": "BRAND"},
        # Generic drugs
        {"name": "Valproic acid", "targets": ["HDAC1", "HDAC2", "HDAC3"], "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Sirolimus",     "targets": ["MTOR", "RPTOR"], "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Itraconazole",  "targets": ["SMO", "PTCH1"], "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Arsenic trioxide", "targets": ["GLI1", "GLI2"], "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Metformin",     "targets": ["PRKAB1", "PRKAB2"], "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Chloroquine",   "targets": ["ATP6V0A1"], "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Hydroxychloroquine", "targets": ["ATP6V0A1"], "has_generic": True, "cost_category": "GENERIC"},
    ]


def _ATRT_GENERIC_FALLBACK_CANDIDATES() -> List[Dict]:
    """Generic-drug-only fallback for resource-limited settings."""
    return [c for c in _ATRT_FALLBACK_CANDIDATES() if c.get("has_generic")]