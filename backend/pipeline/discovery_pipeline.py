"""
discovery_pipeline.py  (v2.0)
==============================
ATRT Drug Repurposing Pipeline — Main Orchestrator

KEY FIXES FROM v1.1
--------------------
1. Generic drug identification is now purely data-driven:
   - Checks pipeline_config.KNOWN_GENERICS (FDA Orange Book seed)
   - Augments with has_generic flag from GenericDrugFetcher (live FDA/RxNorm API)
   - No hardcoded composite scores anywhere

2. All file I/O uses absolute paths from pipeline_config.PATHS.

3. Import guards: every cross-module import uses try/except relative/absolute
   pattern so the module works both as a package and standalone.

4. ATRTGenomicValidator now reads GSE70678 correctly from config path
   and handles the case where neither CBTN nor GSE70678 is present.

5. _apply_atrt_drug_boosts: EZH2 boost only applied once (no double-boosting
   if a drug appears in both known_inhibitors and mechanism_keywords).

6. Escape bypass scoring: curated ATRT resistance map is the primary source;
   STRING-DB PPI neighbors supplement it when available.

GENERIC DRUG PRIORITY
---------------------
When generic_only=True:
  - Only drugs in KNOWN_GENERICS or with has_generic=True are kept
  - Confidence gets a small GENERIC_CONFIDENCE_BONUS (from config)
  - This makes affordable, accessible drugs rank higher in confidence

Top generic candidates for ATRT based on published biology:
  1. Valproic acid   — pan-HDAC; ~0.5 mM IC50 in ATRT; pennies/pill
  2. Vorinostat      — class I/II HDAC; ANDA filed; some generic supply
  3. Sirolimus       — mTOR inhibitor; rapamycin; generic; some CNS data
  4. Itraconazole    — SMO inhibitor (repurposed); generic; MODERATE BBB
  5. Arsenic trioxide— GLI1/2 inhibitor; generic (APL therapy repurposed)
  6. Chloroquine     — autophagy/lysosome; generic; HIGH BBB
  7. Metformin       — AMPK/mTOR; fully generic; very low toxicity
  8. Bortezomib      — proteasome; patent expired 2022; but LOW BBB

References for generic ATRT rationale:
  Valproic acid HDAC: Balasubramanian 2009 Cancer Res
  Itraconazole SMO: Kim 2010 Cancer Cell
  Arsenic trioxide GLI: Beauchamp 2011 Nat Med
  Metformin AMPK/mTOR: Cerezo 2013 Cancer Res
"""

import asyncio
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

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
        BBB as BBB_CONFIG, KNOWN_GENERICS, GENERIC_CONFIDENCE_BONUS,
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
        BBB as BBB_CONFIG, KNOWN_GENERICS, GENERIC_CONFIDENCE_BONUS,
    )

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Generic drug detection
# ─────────────────────────────────────────────────────────────────────────────

def _is_generic(drug: Dict) -> bool:
    """
    Return True if this drug has a generic formulation available.

    Priority:
      1. has_generic flag set by GenericDrugFetcher (live FDA/RxNorm query)
      2. KNOWN_GENERICS seed set (FDA Orange Book, April 2026)
    """
    if drug.get("has_generic"):
        return True
    name = (drug.get("name") or drug.get("drug_name") or "").lower().strip()
    # Strip salt suffixes
    for suffix in (" hydrochloride", " hcl", " sodium", " phosphate",
                   " sulfate", " mesylate", " acetate"):
        name = name.replace(suffix, "")
    return name.strip() in KNOWN_GENERICS


# ─────────────────────────────────────────────────────────────────────────────
# ATRT Resistance Bypass Map (SMARCB1-loss context)
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
    "HDAC3":  ["BCL2", "MYC"],
    "CDK4":   ["PIK3CA", "MTOR", "AKT1"],
    "CDK6":   ["PIK3CA", "MTOR"],
    "AURKA":  ["CDK4", "MYC", "BCL2L1", "MYCN"],
    "MTOR":   ["PIK3CA", "AKT1"],
    "PIK3CA": ["MTOR", "AKT1"],
    "PSMB5":  ["ATG5", "BECN1", "MCL1"],
    "MYC":    ["AURKA", "CDK4", "BCL2L1"],
    "BCL2":   ["MCL1", "BCL2L1"],
    "PARP1":  ["RAD51", "BRCA1", "ATM"],
    "SMO":    ["GLI2", "PIK3CA"],
    "GLI2":   ["MYC", "CDK4"],
    # Generic drug targets
    "PRKAB1": ["MTOR", "PIK3CA"],          # metformin via AMPK
    "ATP6V0A1":["BECN1", "MCL1"],          # chloroquine
    "HDAC1":  ["BCL2", "BCL2L1", "MCL1"], # valproic acid / vorinostat
    "IDO1":   ["MYC", "PIK3CA"],           # indoximod
}


def compute_escape_bypass_score(
    drug_targets: List[str],
    upregulated_genes: Set[str],
) -> float:
    """Score resistance bypass potential for a drug's target set."""
    if not drug_targets:
        return ESCAPE["empty_target_score"]

    constitutive = ESCAPE["constitutive_resistance_nodes"]
    active_bypass: Set[str] = set()
    covered: Set[str] = set()
    has_rna = len(upregulated_genes) >= ESCAPE["min_rna_genes_for_string_weight"]

    for target in drug_targets:
        t_upper = target.upper()
        bypass_candidates = ATRT_RESISTANCE_BYPASS_MAP.get(t_upper, [])
        if bypass_candidates:
            covered.add(t_upper)
            for node in bypass_candidates:
                n_upper = node.upper()
                if n_upper in constitutive:
                    active_bypass.add(n_upper)
                elif has_rna and n_upper in upregulated_genes:
                    active_bypass.add(n_upper)
                elif not has_rna:
                    active_bypass.add(n_upper)

    if not covered:
        return ESCAPE["no_target_score"]

    n_bypass = len(active_bypass)
    scores   = ESCAPE["bypass_scores"]
    max_key  = max(scores.keys())
    return scores.get(n_bypass, scores[max_key])


# ─────────────────────────────────────────────────────────────────────────────
# Genomic Validator (CBTN → GSE70678 fallback → prevalence prior)
# ─────────────────────────────────────────────────────────────────────────────

class ATRTGenomicValidator:
    def __init__(self):
        self.data_dir = Path(PATHS["genomics"])
        self._gse_path = Path(PATHS["scrna"])

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

        if cna_path.exists() or rna_path.exists():
            if cna_path.exists():
                result = self._load_cna(cna_path, result)
            if rna_path.exists():
                result = self._load_rna(rna_path, result)
            return result

        if self._gse_path.exists():
            return self._load_gse70678_fallback(result)

        # No data files — use prevalence priors from Hasselblatt 2011
        logger.info(
            "No genomic data files found. Using published prevalence priors.\n"
            "GSE70678 path: %s (exists=%s)", self._gse_path, self._gse_path.exists()
        )
        result["smarcb1_loss_count"] = 47   # ~95% of 49 samples (Torchia 2015)
        result["total_samples"]      = 49
        result["subgroup_counts"]    = {"TYR": 18, "SHH": 18, "MYC": 13}
        result["calling_method"]     = "prevalence_prior"
        return result

    def _load_cna(self, cna_path, result):
        try:
            cna = pd.read_csv(str(cna_path), sep="\t", index_col=0)
            cna.index = cna.index.astype(str).str.upper().str.strip()
            result["total_samples"] = len(cna.columns)
            threshold = GENOMICS["smarcb1_del_threshold"]
            for alias in GENOMICS["smarcb1_aliases"]:
                if alias in cna.index:
                    row = pd.to_numeric(cna.loc[alias], errors="coerce")
                    result["smarcb1_loss_count"] = int((row <= threshold).sum())
                    break
            for alias in GENOMICS["smarca4_aliases"]:
                if alias in cna.index:
                    row = pd.to_numeric(cna.loc[alias], errors="coerce")
                    result["smarca4_loss_count"] = int((row <= threshold).sum())
                    break
        except Exception as e:
            logger.warning("CNA loading failed: %s", e)
        return result

    def _load_rna(self, rna_path, result):
        try:
            import pandas as pd
            rna = pd.read_csv(str(rna_path), sep="\t", index_col=0)
            up_thr   = GENOMICS["rna_upregulation_zscore"]
            down_thr = GENOMICS["rna_downregulation_zscore"]
            atrt_ind   = GENOMICS["rna_h3k27m_col_indicators"]
            normal_ind = GENOMICS["rna_normal_col_indicators"]
            meta_rows  = GENOMICS["rna_metadata_rows"]

            gene_rows = [r for r in rna.index
                         if str(r).upper().strip() not in meta_rows]
            if len(gene_rows) < GENOMICS["rna_min_genes_required"]:
                return result

            rna_expr  = rna.loc[gene_rows]
            atrt_cols = [c for c in rna_expr.columns
                         if any(ind.lower() in c.lower() for ind in atrt_ind)
                         and not any(n.lower() in c.lower() for n in normal_ind)]
            if not atrt_cols:
                return result

            expr = rna_expr[atrt_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
            result["upregulated_genes"]   = set(expr[expr > up_thr].index)
            result["downregulated_genes"] = set(expr[expr < down_thr].index)
            result["has_rna_data"]        = (
                len(result["upregulated_genes"]) >= ESCAPE["min_rna_genes_for_string_weight"]
            )
            result["calling_method"] = "cbtn_rna"
        except Exception as e:
            logger.warning("RNA loading failed: %s", e)
        return result

    def _load_gse70678_fallback(self, result):
        try:
            import pandas as pd
            df = pd.read_csv(str(self._gse_path), sep="\t", index_col=0, low_memory=False)
            df.index = df.index.astype(str).str.upper().str.strip()

            atrt_ind   = GENOMICS["rna_h3k27m_col_indicators"]
            normal_ind = GENOMICS["rna_normal_col_indicators"]

            atrt_cols = [c for c in df.columns
                         if any(ind.lower() in c.lower() for ind in atrt_ind)
                         and not any(n.lower() in c.lower() for n in normal_ind)]
            normal_cols = [c for c in df.columns
                           if any(n.lower() in c.lower() for n in normal_ind)]

            if not atrt_cols:
                # All columns are ATRT (no normal label) — treat all as tumor
                atrt_cols = list(df.columns)

            result["total_samples"]      = len(atrt_cols)
            result["smarcb1_loss_count"] = len(atrt_cols)   # ~100% in cohort
            result["calling_method"]     = "gse70678_fallback"

            df_num    = df.apply(pd.to_numeric, errors="coerce")
            atrt_mean = df_num[atrt_cols].mean(axis=1)

            if normal_cols:
                diff = atrt_mean - df_num[normal_cols].mean(axis=1)
            else:
                diff = atrt_mean

            up_thr   = GENOMICS["rna_upregulation_zscore"]
            down_thr = GENOMICS["rna_downregulation_zscore"]
            result["upregulated_genes"]   = set(diff[diff > up_thr].index)
            result["downregulated_genes"] = set(diff[diff < down_thr].index)
            result["has_rna_data"]        = (
                len(result["upregulated_genes"]) >= ESCAPE["min_rna_genes_for_string_weight"]
            )
            logger.info(
                "GSE70678 fallback: %d ATRT samples, %d upregulated genes",
                len(atrt_cols), len(result["upregulated_genes"])
            )
        except Exception as e:
            logger.warning("GSE70678 fallback loading failed: %s", e)
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
            "p_value":              None,
            "p_value_label": (
                "N/A — SMARCB1 biallelic loss is the defining event in ATRT "
                "(Hasselblatt 2011 PMID 20625942)"
            ),
            "statistical_note": (
                f"SMARCB1 biallelic loss: {smarcb1}/{total} samples "
                f"({smarcb1 / max(total, 1):.0%}). All are candidates for "
                "EZH2 synthetic lethality (Knutson 2013 PNAS PMID 23620515)."
            ),
        }


# ─────────────────────────────────────────────────────────────────────────────
# ATRT Drug Score Boosts
# ─────────────────────────────────────────────────────────────────────────────

def _apply_atrt_drug_boosts(
    candidates: List[Dict],
    subgroup: Optional[str] = None,
) -> List[Dict]:
    """
    Apply ATRT-specific multiplicative boosts to composite scores.

    Each drug gets AT MOST ONE boost (whichever is largest).
    This prevents double-boosting if a drug matches multiple criteria.
    """
    n_ezh2 = n_aurka = n_smo = 0

    for c in candidates:
        name_lower = (c.get("name") or c.get("drug_name") or "").lower().strip()
        mech_lower = (c.get("mechanism") or "").lower()
        targets    = [t.upper() for t in (c.get("targets") or [])]

        for suffix in (" hydrochloride", " hcl", " sodium", " mesylate", " phosphate"):
            name_lower = name_lower.replace(suffix, "")

        base_score = c.get("score", 0.0)
        c["atrt_boosts_applied"] = []
        c["ezh2_boosted"]  = False
        c["aurka_boosted"] = False

        # Check each inhibitor type — apply the first match only
        is_ezh2 = (
            name_lower.strip() in EZH2_INHIBITOR["known_inhibitors"]
            or any(kw in mech_lower for kw in EZH2_INHIBITOR["mechanism_keywords"])
            or (len(targets) == 1 and targets[0] == "EZH2")
        )
        if is_ezh2:
            mult = EZH2_INHIBITOR["composite_boost"]
            c["score"] = round(min(1.0, base_score * mult), 4)
            c["atrt_boosts_applied"].append(
                f"EZH2-inhibitor ×{mult} (SMARCB1 synthetic lethality — Knutson 2013)"
            )
            c["ezh2_boosted"] = True
            n_ezh2 += 1
            continue   # skip other boosts for this drug

        is_aurka = (
            name_lower.strip() in AURKA_INHIBITOR["known_inhibitors"]
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
            continue

        is_smo = (
            name_lower.strip() in SMO_INHIBITOR["known_inhibitors"]
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
                c["score"] = round(min(1.0, base_score * mult), 4)
                c["atrt_boosts_applied"].append(
                    f"SMO-inhibitor ×{mult} (SHH subgroup — Torchia 2015)"
                )
                n_smo += 1

    logger.info(
        "ATRT boosts applied: EZH2×%d | AURKA×%d | SMO×%d",
        n_ezh2, n_aurka, n_smo
    )
    return candidates


# ─────────────────────────────────────────────────────────────────────────────
# Fallback candidate lists
# ─────────────────────────────────────────────────────────────────────────────

def _atrt_fallback_candidates() -> List[Dict]:
    """Complete curated candidate list used when OpenTargets is unavailable."""
    return [
        # --- EZH2 synthetic lethality ---
        {"name": "Tazemetostat",   "targets": ["EZH2"],
         "mechanism": "EZH2 inhibitor", "has_generic": False, "cost_category": "BRAND"},
        # --- Pan-HDAC ---
        {"name": "Panobinostat",   "targets": ["HDAC1", "HDAC2", "HDAC3", "HDAC6"],
         "mechanism": "pan-HDAC inhibitor", "has_generic": False, "cost_category": "BRAND"},
        # --- AURKA ---
        {"name": "Alisertib",      "targets": ["AURKA"],
         "mechanism": "aurora kinase A inhibitor", "has_generic": False, "cost_category": "INVESTIGATIONAL"},
        # --- BET ---
        {"name": "Birabresib",     "targets": ["BRD4", "BRD2", "BRD3"],
         "mechanism": "BET bromodomain inhibitor", "has_generic": False, "cost_category": "INVESTIGATIONAL"},
        # --- CDK4/6 ---
        {"name": "Abemaciclib",    "targets": ["CDK4", "CDK6"],
         "mechanism": "CDK4/CDK6 inhibitor", "has_generic": False, "cost_category": "BRAND"},
        # --- Proteasome ---
        {"name": "Marizomib",      "targets": ["PSMB5", "PSMB2", "PSMB1"],
         "mechanism": "proteasome inhibitor", "has_generic": False, "cost_category": "INVESTIGATIONAL"},
        # --- SHH ---
        {"name": "Vismodegib",     "targets": ["SMO"],
         "mechanism": "smoothened inhibitor", "has_generic": False, "cost_category": "BRAND"},
        # --- DRD2/TRAIL ---
        {"name": "ONC201",         "targets": ["DRD2", "CLPB"],
         "mechanism": "DRD2 antagonist / TRAIL inducer", "has_generic": False, "cost_category": "BRAND"},
        # --- PI3K/mTOR ---
        {"name": "Paxalisib",      "targets": ["PIK3CA", "PIK3CD", "MTOR"],
         "mechanism": "PI3K inhibitor (CNS-penetrant)", "has_generic": False, "cost_category": "INVESTIGATIONAL"},
        # =========================================================
        # GENERIC DRUGS — ATRT biological rationale
        # =========================================================
        {"name": "Valproic acid",  "targets": ["HDAC1", "HDAC2", "HDAC3"],
         "mechanism": "HDAC inhibitor (repurposed anticonvulsant)",
         "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Vorinostat",     "targets": ["HDAC1", "HDAC2", "HDAC3", "HDAC6"],
         "mechanism": "pan-HDAC inhibitor",
         "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Sirolimus",      "targets": ["MTOR", "RPTOR"],
         "mechanism": "mTOR inhibitor (rapamycin)",
         "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Itraconazole",   "targets": ["SMO", "PTCH1"],
         "mechanism": "SMO inhibitor (repurposed antifungal)",
         "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Arsenic trioxide","targets": ["GLI1", "GLI2"],
         "mechanism": "GLI inhibitor (repurposed APL therapy)",
         "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Chloroquine",    "targets": ["ATP6V0A1", "BECN1"],
         "mechanism": "autophagy inhibitor (repurposed antimalarial)",
         "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Hydroxychloroquine","targets": ["ATP6V0A1", "BECN1"],
         "mechanism": "autophagy inhibitor (repurposed antimalarial)",
         "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Metformin",      "targets": ["PRKAB1", "PRKAB2"],
         "mechanism": "AMPK activator / mTOR inhibitor (repurposed antidiabetic)",
         "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Bortezomib",     "targets": ["PSMB5", "PSMB1"],
         "mechanism": "proteasome inhibitor (boronic acid, IV)",
         "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Tretinoin",      "targets": ["RARA", "RARB", "RARG"],
         "mechanism": "retinoid receptor agonist / differentiation therapy",
         "has_generic": True, "cost_category": "GENERIC"},
        {"name": "Temozolomide",   "targets": ["MGMT"],
         "mechanism": "alkylating agent (standard of care for CNS tumors)",
         "has_generic": True, "cost_category": "GENERIC"},
    ]


def _generic_only_fallback() -> List[Dict]:
    return [c for c in _atrt_fallback_candidates() if c.get("has_generic")]


# ─────────────────────────────────────────────────────────────────────────────
# Production Pipeline
# ─────────────────────────────────────────────────────────────────────────────

import pandas as pd  # noqa: E402  (after fallback defs which use pd in methods)


class ProductionPipeline:
    """
    ATRT Drug Repurposing Pipeline — Production Orchestrator v2.0

    Parameters
    ----------
    generic_only : bool
        If True, restrict output to drugs with FDA-approved generic formulations.
        These are more accessible/affordable for pediatric patients in all settings.
    """

    def __init__(self, generic_only: bool = False):
        self.generic_only  = generic_only
        self._data_fetcher = ProductionDataFetcher()
        self._cmap         = CMAPQuery()
        self._ppi          = PPINetwork()
        self._depmap       = DepMapEssentiality()
        self._tissue       = TissueExpressionScorer("atrt")
        self._synergy      = SynergyPredictor()
        self._hyp_gen      = HypothesisGenerator()
        self._genomic_val  = ATRTGenomicValidator()
        self._stat_val     = StatisticalValidator()
        self._bbb_filter   = BBBFilter()
        logger.info("ProductionPipeline v2.0 (generic_only=%s)", generic_only)

    async def initialize(self, disease: str = "atrt"):
        return True

    async def analyze_disease(
        self,
        disease_name: str = "atrt",
        min_score:    float = 0.2,
        max_results:  int   = 20,
        generic_only: bool  = False,
    ) -> Dict:
        result = await self.run(
            disease_name=disease_name,
            top_k=max_results,
            generic_only=generic_only or self.generic_only,
        )
        candidates = result.get("top_candidates", [])
        for c in candidates:
            c.setdefault("indication", "ATRT / rhabdoid tumor")
            c.setdefault("mechanism", "")
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
        effective_generic = generic_only or self.generic_only
        if subgroup:
            self._tissue.set_subgroup(subgroup)

        logger.info(
            "Pipeline run: disease=%s subgroup=%s location=%s generic_only=%s top_k=%d",
            disease_name, subgroup or "pan-ATRT", location, effective_generic, top_k
        )

        # ── Step 1: Fetch candidates ──────────────────────────────────────────
        disease_data = await self._data_fetcher.fetch_disease_data(disease_name)
        try:
            candidates = await self._data_fetcher.fetch_approved_drugs(
                annotate_generics=True
            )
        except Exception as e:
            logger.warning("OpenTargets fetch failed (%s) — using curated fallback.", e)
            candidates = []

        if not candidates:
            logger.info("Using curated ATRT fallback candidate list.")
            candidates = _atrt_fallback_candidates()

        # Merge curated list to ensure key drugs are always present
        existing_names = {c.get("name", "").lower() for c in candidates}
        for fb in _atrt_fallback_candidates():
            if fb["name"].lower() not in existing_names:
                candidates.append(fb)

        # Generic filter (applied early)
        if effective_generic:
            all_n = len(candidates)
            candidates = [c for c in candidates if _is_generic(c)]
            logger.info(
                "Generic filter: %d → %d candidates", all_n, len(candidates)
            )
            if len(candidates) < 5:
                logger.warning(
                    "Too few generic candidates (%d) — adding curated generics.",
                    len(candidates)
                )
                gen_fallback = _generic_only_fallback()
                gen_names = {c.get("name", "").lower() for c in candidates}
                for fb in gen_fallback:
                    if fb["name"].lower() not in gen_names:
                        candidates.append(fb)

        logger.info("Scoring %d candidates...", len(candidates))

        # ── Step 2: Genomic data ──────────────────────────────────────────────
        genomic_stats = self._genomic_val.validate_atrt_cohort()
        smarcb1_stats = self._genomic_val.calculate_smarcb1_statistics(genomic_stats)
        upregulated   = genomic_stats.get("upregulated_genes", set())
        has_rna_data  = genomic_stats.get("has_rna_data", False)

        logger.info(
            "Genomics: method=%s | SMARCB1_loss=%d/%d | upregulated_genes=%d",
            genomic_stats.get("calling_method"),
            genomic_stats.get("smarcb1_loss_count", 0),
            genomic_stats.get("total_samples", 0),
            len(upregulated),
        )

        # ── Step 3: Escape bypass ─────────────────────────────────────────────
        ew = ESCAPE
        for drug in candidates:
            targets = drug.get("targets") or []
            curated_score = compute_escape_bypass_score(targets, upregulated)

            # Live STRING-DB neighbors (via PPI module cache)
            live_hits = []
            for t in targets:
                for nb in self._ppi.get_neighbors(t.upper()):
                    if nb.upper() in upregulated:
                        live_hits.append(nb)
            live_hits = list(set(live_hits))
            drug["resistance_nodes"] = live_hits

            if live_hits and has_rna_data:
                drug["escape_bypass_score"] = round(
                    ew["string_weight"] * ew["string_hit_score"]
                    + ew["curated_weight"] * curated_score, 4
                )
                drug["escape_note"] = f"STRING+RNA ({len(live_hits)} bypass hits)"
            else:
                drug["escape_bypass_score"] = round(curated_score, 4)
                drug["escape_note"] = "Curated ATRT resistance bypass"

        # ── Step 4: Multi-omic scoring ────────────────────────────────────────
        candidates = await self._tissue.score_batch(candidates)
        candidates = await self._depmap.score_batch(candidates, disease_name)
        candidates = await self._ppi.score_batch(candidates, disease_data.get("genes", []))

        # ── Step 5: Composite score ───────────────────────────────────────────
        cw = COMPOSITE_WEIGHTS
        sd = SCORE_DEFAULTS
        for c in candidates:
            t = c.get("tissue_expression_score", sd["tissue_expression_score"])
            d = c.get("depmap_score",            sd["depmap_score"])
            p = c.get("ppi_score",               sd["ppi_score"])
            e = c.get("escape_bypass_score",     sd["escape_bypass_score"])
            c["score"] = round(
                t * cw["tissue"] + d * cw["depmap"]
                + e * cw["escape"] + p * cw["ppi"], 4
            )

        # ── Step 6: ATRT-specific boosts ──────────────────────────────────────
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

        # ── Step 8: Generic confidence bonus ──────────────────────────────────
        if effective_generic:
            for c in candidates:
                if _is_generic(c):
                    c["generic_confidence_bonus"] = GENERIC_CONFIDENCE_BONUS

        # ── Step 9: Sort and annotate IC50 ───────────────────────────────────
        sorted_candidates = sorted(
            candidates, key=lambda x: x.get("score", 0), reverse=True
        )

        try:
            try:
                from .published_ic50_atrt_validation import annotate_atrt_candidates_with_ic50
            except ImportError:
                from published_ic50_atrt_validation import annotate_atrt_candidates_with_ic50
            sorted_candidates = annotate_atrt_candidates_with_ic50(sorted_candidates)
        except Exception as e:
            logger.debug("IC50 annotation skipped: %s", e)

        # ── Step 10: CMap integration ─────────────────────────────────────────
        if self._cmap.has_precomputed_scores():
            for c in sorted_candidates:
                cmap_score = self._cmap.get_precomputed_score(c.get("name", ""))
                if cmap_score is not None:
                    c["cmap_score"]  = cmap_score
                    c["is_reverser"] = cmap_score > 0.75

        # ── Step 11: Hypothesis generation ───────────────────────────────────
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

        n_generic_out = sum(1 for c in sorted_candidates[:top_k] if _is_generic(c))

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
                "n_generic_in_top_k":  n_generic_out,
                "subgroup":             subgroup or "pan-ATRT",
                "location":             location,
                "composite_weights":    cw,
                "escape_bypass_mode":   "RNA-confirmed" if has_rna_data else "curated fallback",
                "rna_upregulated_genes": len(upregulated),
                "n_ezh2_boosted":       sum(1 for c in sorted_candidates if c.get("ezh2_boosted")),
                "n_aurka_boosted":      sum(1 for c in sorted_candidates if c.get("aurka_boosted")),
                "calling_method":       genomic_stats.get("calling_method", "none"),
                "generic_only":         effective_generic,
                "depmap_source":        (
                    "live CSV" if not self._depmap.using_fallback
                    else "verified fallback (Knutson 2013 + others)"
                ),
                "tissue_source":        (
                    "GSE70678 + GTEx" if self._tissue._diff_scores
                    else "curated ATRT scores only"
                ),
            },
        }