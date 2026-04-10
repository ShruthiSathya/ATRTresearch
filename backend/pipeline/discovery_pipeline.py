"""
discovery_pipeline.py
======================
ATRT Drug Repurposing Pipeline — Main Orchestrator

ATRT SCORING ARCHITECTURE
--------------------------
Same weight structure as DIPG v5.5 but with ATRT-specific data sources:

  composite_score = 0.40 × tissue_expression      (GSE70678 bulk RNA, Torchia 2015)
                  + 0.35 × depmap_essentiality     (DepMap, ATRT/rhabdoid lines)
                  + 0.20 × escape_bypass           (SMARCB1-loss resistance map)
                  + 0.05 × ppi_network             (STRING-DB proximity)

KEY DIFFERENCES FROM DIPG PIPELINE
-------------------------------------
1. EZH2 inhibitors BOOSTED (×1.40) — synthetic lethality with SMARCB1 loss
   (vs DIPG where they are penalised ×0.50 due to H3K27M suppression)
2. AURKA inhibitors boosted — MYCN stabilisation in ATRT-MYC subgroup
3. BBB penalty less severe — ATRT is NOT exclusively brainstem
4. Tissue data: GSE70678 bulk RNA (not Filbin 2018 scRNA)
5. DepMap cell lines: BT16/BT37/G401/A204 (not GBM lines)
6. Genomics: SMARCB1 loss calling (not H3K27M mutation calling)

REFERENCES
-----------
Torchia 2015 : Cancer Cell 30:891. [GSE70678]
Knutson 2013 : PNAS 110:7922. [EZH2 synthetic lethality]
Johann 2016  : Cancer Cell 29:379. [subgroups]
Frühwald 2020: CNS Oncol 9:CNS56. [ATRT overview]
Behan 2019   : Nature 568:511. [DepMap]
"""

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
    EZH2_INHIBITOR,
    AURKA_INHIBITOR,
    SMO_INHIBITOR,
    OPENTARGETS,
    BBB as BBB_CONFIG,
)

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# ATRT RESISTANCE BYPASS NETWORK
# Based on SMARCB1-loss biology and published resistance mechanisms.
# Key difference from DIPG: bypass routes run through PRC2/BET/AURKA axis.
#
# Sources:
#   Wilson & Roberts 2011 Nat Rev Cancer — SMARCB1 biology
#   Knutson 2013 PNAS — EZH2 resistance bypass (CDK4, BRD4)
#   Frühwald 2020 CNS Oncology — resistance overview
#   Lin 2019 Sci Transl Med — proteasome-autophagy crosstalk
# ─────────────────────────────────────────────────────────────────────────────

ATRT_RESISTANCE_BYPASS_MAP: Dict[str, List[str]] = {
    # EZH2/PRC2 — resistance via BET bromodomain compensation
    "EZH2":   ["BRD4", "CDK4", "MYC", "MYCN"],
    "EED":    ["BRD4", "CDK4", "MYC"],
    "SUZ12":  ["BRD4", "MYC"],

    # BET bromodomain — resistance via CDK-mediated transcription
    "BRD4":   ["CDK4", "MYC", "MYCN", "MTOR"],
    "BRD2":   ["MYC", "CDK4", "EZH2"],
    "BRD3":   ["MYC", "EZH2"],

    # HDAC — resistance via anti-apoptotic upregulation
    "HDAC1":  ["BCL2", "BCL2L1", "MCL1", "MYC"],
    "HDAC2":  ["BCL2", "MCL1", "MYC"],
    "HDAC3":  ["BCL2", "NFKB1", "MYC"],

    # CDK4/6 — resistance via PI3K/mTOR activation
    "CDK4":   ["PIK3CA", "MTOR", "AKT1", "CCNE1"],
    "CDK6":   ["PIK3CA", "MTOR", "AKT1"],

    # Aurora kinase — resistance via MYC/CDK4 backup
    "AURKA":  ["CDK4", "MYC", "BCL2L1", "MYCN"],
    "AURKB":  ["CDK4", "MYC"],

    # mTOR / PI3K
    "MTOR":   ["PIK3CA", "AKT1", "EIF4E"],
    "PIK3CA": ["MTOR", "AKT1", "MAPK1"],

    # Proteasome — resistance via autophagy induction
    "PSMB5":  ["ATG5", "BECN1", "SQSTM1", "MCL1"],
    "PSMB2":  ["ATG5", "BECN1", "MCL1"],

    # MYC/MYCN axis
    "MYC":    ["AURKA", "CDK4", "BCL2L1", "BRD4"],
    "MYCN":   ["AURKA", "CDK4", "BRD4"],

    # Anti-apoptotic
    "BCL2":   ["MCL1", "BCL2L1"],
    "BCL2L1": ["MCL1", "BCL2"],
    "MCL1":   ["BCL2", "BCL2L1"],

    # SHH subgroup — resistance via PI3K when SMO blocked
    "SMO":    ["GLI2", "PIK3CA", "MTOR"],
    "GLI2":   ["MYC", "CDK4", "PIK3CA"],

    # PARP
    "PARP1":  ["RAD51", "BRCA1", "ATM"],
}

# Constitutive resistance nodes in ATRT — always active
ATRT_CONSTITUTIVE_RESISTANCE: Set[str] = {
    "MYC", "MYCN", "BCL2L1", "MCL1", "CDK4", "BRD4", "PIK3CA", "EZH2",
}


def compute_escape_bypass_score(
    drug_targets: List[str],
    upregulated_genes: Set[str],
) -> float:
    """
    Compute SMARCB1-loss escape bypass score for an ATRT drug candidate.

    Logic mirrors DIPG v5.5 escape bypass but uses ATRT_RESISTANCE_BYPASS_MAP
    and ATRT_CONSTITUTIVE_RESISTANCE instead of DIPG-specific maps.

    Score interpretation:
        1.00 = no active bypass routes (drug should work)
        0.40 = 4+ active bypass routes (high resistance risk)
    """
    if not drug_targets:
        return ESCAPE["empty_target_score"]

    constitutive  = ESCAPE["constitutive_resistance_nodes"]
    active_bypass: Set[str] = set()
    covered:       Set[str] = set()

    has_rna_data = len(upregulated_genes) >= ESCAPE["min_rna_genes_for_string_weight"]

    for target in drug_targets:
        target_upper      = target.upper()
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

    n_bypass  = len(active_bypass)
    scores    = ESCAPE["bypass_scores"]
    max_key   = max(scores.keys())
    score     = scores.get(n_bypass, scores[max_key])

    return score


# ─────────────────────────────────────────────────────────────────────────────
# ATRT Genomic Validator
# Loads CBTN ATRT data and computes SMARCB1 loss statistics.
# Analogous to PedcBioPortalValidator in DIPG pipeline.
# ─────────────────────────────────────────────────────────────────────────────

class ATRTGenomicValidator:
    """
    Validates ATRT genomic data from CBTN Cavatica or GSE70678 fallback.

    Outputs:
      smarcb1_loss_count  : samples with SMARCB1 homozygous deletion
      smarca4_loss_count  : SMARCA4 deletion (~5% of ATRT)
      upregulated_genes   : set of RNA-upregulated genes vs normal
      total_samples       : cohort size
      has_rna_data        : bool — enough genes for RNA-confirmed bypass

    Data: CBTN ATRT from Cavatica (controlled access)
    URL: https://portal.kidsfirstdrc.org
    Fallback: GSE70678 expression data when CBTN unavailable
    """

    def __init__(self, data_dir: str = None):
        self.data_dir = Path(data_dir or PATHS["genomics"])

    def validate_atrt_cohort(self) -> Dict:
        result = {
            "smarcb1_loss_count": 0,
            "smarca4_loss_count": 0,
            "upregulated_genes":  set(),
            "downregulated_genes": set(),
            "total_samples":      0,
            "has_rna_data":       False,
            "subgroup_counts":    {"TYR": 0, "SHH": 0, "MYC": 0},
            "calling_method":     "none",
        }

        cna_path = self.data_dir / "cna.txt"
        rna_path = self.data_dir / "rna_zscores.txt"

        if not cna_path.exists() and not rna_path.exists():
            logger.info(
                "CBTN ATRT genomic files not found at %s.\n"
                "Falling back to GSE70678 for RNA data.\n"
                "To enable CBTN data: https://portal.kidsfirstdrc.org",
                self.data_dir,
            )
            return self._load_gse70678_fallback(result)

        # Load CNA for SMARCB1 deletion
        if cna_path.exists():
            result = self._load_cna(cna_path, result)

        # Load RNA for upregulated genes
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
                    result["smarcb1_loss_count"] = int(
                        (row <= threshold).sum()
                    )
                    logger.info(
                        "SMARCB1 loss: %d/%d samples (CNA ≤ %d)",
                        result["smarcb1_loss_count"],
                        result["total_samples"],
                        threshold,
                    )
                    break

            for alias in GENOMICS["smarca4_aliases"]:
                if alias in cna.index:
                    row = pd.to_numeric(cna.loc[alias], errors="coerce")
                    result["smarca4_loss_count"] = int(
                        (row <= threshold).sum()
                    )
                    break

        except Exception as e:
            logger.warning("CNA loading failed: %s", e)

        return result

    def _load_rna(self, rna_path: Path, result: Dict) -> Dict:
        try:
            rna = pd.read_csv(rna_path, sep="\t", index_col=0)

            threshold       = GENOMICS["rna_upregulation_zscore"]
            down_threshold  = GENOMICS["rna_downregulation_zscore"]
            atrt_indicators = GENOMICS["rna_h3k27m_col_indicators"]  # "ATRT", "MRT"...
            normal_indicators = GENOMICS["rna_normal_col_indicators"]
            metadata_rows   = GENOMICS["rna_metadata_rows"]
            min_genes       = GENOMICS["rna_min_genes_required"]

            gene_rows = [
                r for r in rna.index
                if str(r).upper().strip() not in metadata_rows
            ]

            if len(gene_rows) < min_genes:
                logger.warning(
                    "rna_zscores.txt: only %d gene rows — likely a summary file. "
                    "RNA-confirmed bypass disabled. Need a full expression matrix.",
                    len(gene_rows)
                )
                return result

            rna_expr = rna.loc[gene_rows]

            atrt_cols = [
                c for c in rna_expr.columns
                if any(ind.lower() in c.lower() for ind in atrt_indicators)
                and not any(n.lower() in c.lower() for n in normal_indicators)
            ]

            if not atrt_cols:
                logger.warning(
                    "No ATRT indicator columns in rna_zscores.txt. "
                    "Column sample: %s", list(rna_expr.columns[:5])
                )
                return result

            expr_numeric = rna_expr[atrt_cols].apply(
                pd.to_numeric, errors="coerce"
            ).mean(axis=1)

            result["upregulated_genes"]   = set(expr_numeric[expr_numeric > threshold].index)
            result["downregulated_genes"] = set(expr_numeric[expr_numeric < down_threshold].index)
            result["has_rna_data"]        = (
                len(result["upregulated_genes"]) >= ESCAPE["min_rna_genes_for_string_weight"]
            )
            result["calling_method"]      = "cbtn_rna"

            logger.info(
                "CBTN ATRT RNA: %d upregulated, %d downregulated genes "
                "across %d ATRT samples",
                len(result["upregulated_genes"]),
                len(result["downregulated_genes"]),
                len(atrt_cols),
            )

        except Exception as e:
            logger.warning("RNA loading failed: %s", e)

        return result

    def _load_gse70678_fallback(self, result: Dict) -> Dict:
        """Load GSE70678 as RNA data source when CBTN unavailable."""
        gse_path = Path(PATHS["scrna"])
        if not gse_path.exists():
            logger.info(
                "GSE70678 not found. Using prevalence priors only.\n"
                "Download: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70678"
            )
            # Apply published prevalence priors — Johann 2016, n=150
            result["subgroup_counts"] = {
                "TYR": 54,   # 36% × 150
                "SHH": 55,   # 37% × 150
                "MYC": 41,   # 27% × 150
            }
            result["calling_method"] = "prevalence_prior"
            return result

        try:
            df = pd.read_csv(gse_path, sep="\t", index_col=0)
            df.index = df.index.astype(str).str.upper().str.strip()

            atrt_indicators   = GENOMICS["rna_h3k27m_col_indicators"]
            normal_indicators = GENOMICS["rna_normal_col_indicators"]

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
                    "No ATRT columns in GSE70678. Check column indicators. "
                    "Sample columns: %s", list(df.columns[:10])
                )
                return result

            result["total_samples"]    = len(atrt_cols)
            result["smarcb1_loss_count"] = len(atrt_cols)   # all ATRT are SMARCB1-null
            result["calling_method"]   = "gse70678_fallback"

            threshold = GENOMICS["rna_upregulation_zscore"]
            down_threshold = GENOMICS["rna_downregulation_zscore"]

            df_num    = df.apply(pd.to_numeric, errors="coerce")
            atrt_mean = df_num[atrt_cols].mean(axis=1)

            if normal_cols:
                normal_mean = df_num[normal_cols].mean(axis=1)
                diff        = atrt_mean - normal_mean
            else:
                diff = atrt_mean

            result["upregulated_genes"]   = set(diff[diff > threshold].index)
            result["downregulated_genes"] = set(diff[diff < down_threshold].index)
            result["has_rna_data"]        = (
                len(result["upregulated_genes"]) >= ESCAPE["min_rna_genes_for_string_weight"]
            )

            logger.info(
                "GSE70678 fallback: %d ATRT samples, %d upregulated genes, "
                "%d downregulated genes",
                len(atrt_cols),
                len(result["upregulated_genes"]),
                len(result["downregulated_genes"]),
            )

        except Exception as e:
            logger.warning("GSE70678 fallback loading failed: %s", e)

        return result

    def calculate_smarcb1_statistics(self, genomic_stats: Dict) -> Dict:
        """
        Compute summary statistics for SMARCB1 loss in the cohort.
        Used for hypothesis report generation.
        """
        total    = genomic_stats.get("total_samples", 0)
        smarcb1  = genomic_stats.get("smarcb1_loss_count", 0)
        smarca4  = genomic_stats.get("smarca4_loss_count", 0)

        return {
            "smarcb1_loss_count":     smarcb1,
            "smarca4_loss_count":     smarca4,
            "total_loss_count":       smarcb1 + smarca4,
            "total_samples":          total,
            "smarcb1_prevalence":     smarcb1 / max(total, 1),
            "combined_prevalence":    (smarcb1 + smarca4) / max(total, 1),
            "p_value":                None,   # No Fisher's test for ATRT (no 2×2 table)
            "p_value_label":          "N/A — SMARCB1 loss is the defining event in ATRT",
            "statistical_note": (
                f"SMARCB1 biallelic loss in {smarcb1}/{total} samples "
                f"({smarcb1/max(total,1):.0%}). "
                "No co-occurrence test — SMARCB1 loss is the defining ATRT event, "
                "not a co-occurrence hypothesis."
            ),
        }


# ─────────────────────────────────────────────────────────────────────────────
# EZH2 / AURKA / SMO Boost Application
# These are applied after composite scoring as disease-specific adjustments.
# ─────────────────────────────────────────────────────────────────────────────

def _apply_atrt_drug_boosts(
    candidates: List[Dict],
    subgroup: Optional[str] = None,
) -> List[Dict]:
    """
    Apply ATRT-specific score boosts/penalties based on disease biology.

    1. EZH2 inhibitors: BOOST ×1.40 (SMARCB1 synthetic lethality)
    2. AURKA inhibitors: BOOST ×1.30 in MYC subgroup, ×1.15 otherwise
    3. SMO/GLI inhibitors: BOOST ×1.25 in SHH subgroup, neutral otherwise

    All boosts are capped at 1.0 final score.
    Applied BEFORE BBB penalty so BBB correctly down-weights the boosted score.
    """
    n_ezh2 = n_aurka = n_smo = 0

    for c in candidates:
        name_lower = (c.get("name") or c.get("drug_name") or "").lower().strip()
        mech_lower = (c.get("mechanism") or "").lower()
        targets    = [t.upper() for t in (c.get("targets") or [])]

        # Remove common salt suffixes for matching
        for suffix in (" hydrochloride", " hcl", " sodium", " mesylate", " phosphate"):
            name_lower = name_lower.replace(suffix, "")

        base_score = c.get("score", 0.0)
        c["atrt_boosts_applied"] = []

        # ── EZH2 inhibitor BOOST ────────────────────────────────────────────────
        is_ezh2 = (
            name_lower in EZH2_INHIBITOR["known_inhibitors"]
            or any(kw in mech_lower for kw in EZH2_INHIBITOR["mechanism_keywords"])
            or (targets == ["EZH2"])
        )
        if is_ezh2:
            mult = EZH2_INHIBITOR["composite_boost"]
            c["score"] = round(min(1.0, base_score * mult), 4)
            c["atrt_boosts_applied"].append(
                f"EZH2-inhibitor ×{mult} (SMARCB1 synthetic lethality)"
            )
            c["ezh2_boosted"] = True
            n_ezh2 += 1
            continue   # EZH2 boost is exclusive — skip other boosts for this drug

        c["ezh2_boosted"] = False

        # ── AURKA inhibitor BOOST ───────────────────────────────────────────────
        is_aurka = (
            name_lower in AURKA_INHIBITOR["known_inhibitors"]
            or any(kw in mech_lower for kw in AURKA_INHIBITOR["mechanism_keywords"])
            or "AURKA" in targets
        )
        if is_aurka:
            if subgroup and subgroup.upper() == "MYC":
                mult = AURKA_INHIBITOR["composite_boost_myc_subgroup"]
            else:
                mult = AURKA_INHIBITOR["composite_boost_other"]
            c["score"] = round(min(1.0, base_score * mult), 4)
            c["atrt_boosts_applied"].append(
                f"AURKA-inhibitor ×{mult} (MYCN stabilisation"
                + (", MYC subgroup" if subgroup == "MYC" else "") + ")"
            )
            c["aurka_boosted"] = True
            n_aurka += 1
        else:
            c["aurka_boosted"] = False

        # ── SMO/GLI inhibitor BOOST (SHH subgroup) ─────────────────────────────
        is_smo = (
            name_lower in SMO_INHIBITOR["known_inhibitors"]
            or any(kw in mech_lower for kw in SMO_INHIBITOR["mechanism_keywords"])
            or bool({"SMO", "GLI2", "GLI1", "PTCH1"} & set(targets))
        )
        if is_smo:
            if subgroup and subgroup.upper() == "SHH":
                mult = SMO_INHIBITOR["composite_boost_shh_subgroup"]
            else:
                mult = SMO_INHIBITOR["composite_boost_other"]
            if mult > 1.0:
                c["score"] = round(min(1.0, c.get("score", base_score) * mult), 4)
                c["atrt_boosts_applied"].append(
                    f"SMO-inhibitor ×{mult} (SHH subgroup)"
                )
                n_smo += 1

    logger.info(
        "ATRT drug boosts applied: EZH2×%d | AURKA×%d | SMO×%d",
        n_ezh2, n_aurka, n_smo,
    )
    return candidates


# ─────────────────────────────────────────────────────────────────────────────
# Production Pipeline
# ─────────────────────────────────────────────────────────────────────────────

class ProductionPipeline:
    """
    ATRT Drug Repurposing Pipeline — Production Orchestrator.

    Scoring pipeline:
      1. Fetch 557 CNS/Oncology drugs from OpenTargets API
      2. Load ATRT genomic data (CBTN or GSE70678 fallback)
      3. Score tissue expression using GSE70678 bulk RNA
      4. Score DepMap CRISPR essentiality (ATRT/rhabdoid lines)
      5. Score PPI network proximity (STRING-DB)
      6. Compute escape bypass score (SMARCB1-loss resistance map)
      7. Compute composite score
      8. Apply ATRT-specific boosts (EZH2, AURKA, SMO)
      9. Apply location-aware BBB penalty
      10. Annotate with published ATRT IC50 data
      11. Generate triple combination hypothesis
    """

    def __init__(self):
        self._data_fetcher    = ProductionDataFetcher()
        self._cmap            = CMAPQuery()
        self._ppi             = PPINetwork()
        self._depmap          = DepMapEssentiality()
        self._tissue          = TissueExpressionScorer("atrt")
        self._synergy         = SynergyPredictor()
        self._hyp_gen         = HypothesisGenerator()
        self._genomic_validator = ATRTGenomicValidator()
        self._stat_validator  = StatisticalValidator()
        self._bbb_filter      = BBBFilter()

    async def initialize(self, disease: str = "atrt"):
        return True

    async def run(
        self,
        disease_name: str = "atrt",
        top_k: int = 20,
        subgroup: Optional[str] = None,
        location: str = "unknown_location",
    ) -> Dict:
        """
        Run the full ATRT pipeline.

        Parameters
        ----------
        disease_name : str
        top_k        : int — number of top candidates to return
        subgroup     : "TYR", "SHH", "MYC", or None (pan-ATRT)
        location     : "infratentorial", "supratentorial", "unknown_location"
        """
        # Set subgroup for tissue expression scorer
        if subgroup:
            self._tissue.set_subgroup(subgroup)

        # ── Step 1: Fetch disease data and drug candidates ───────────────────
        disease_data = await self._data_fetcher.fetch_disease_data(disease_name)
        candidates   = await self._data_fetcher.fetch_approved_drugs()

        if not candidates:
            logger.warning("OpenTargets returned 0 drugs. Using fallback library.")
            candidates = _ATRT_FALLBACK_CANDIDATES()

        # ── Step 2: Load ATRT genomic data ───────────────────────────────────
        genomic_stats   = self._genomic_validator.validate_atrt_cohort()
        smarcb1_stats   = self._genomic_validator.calculate_smarcb1_statistics(
            genomic_stats
        )
        upregulated     = genomic_stats.get("upregulated_genes", set())
        downregulated   = genomic_stats.get("downregulated_genes", set())
        has_rna_data    = genomic_stats.get("has_rna_data", False)

        logger.info(
            "ATRT RNA data: %s | %d upregulated, %d downregulated genes\n"
            "SMARCB1 loss: %d samples | Escape bypass mode: %s",
            "loaded" if has_rna_data else "not loaded (using curated fallback)",
            len(upregulated), len(downregulated),
            genomic_stats.get("smarcb1_loss_count", 0),
            "RNA-confirmed" if has_rna_data else "curated fallback",
        )

        # ── Step 3: Escape bypass scoring ─────────────────────────────────────
        ew = ESCAPE
        for drug in candidates:
            targets = drug.get("targets", [])

            # STRING-DB live escape hits
            live_escape_hits = []
            for t in targets:
                neighbors = self._ppi.get_neighbors(t)
                live_escape_hits.extend(
                    [n for n in neighbors if n in upregulated]
                )
            drug["resistance_nodes"] = list(set(live_escape_hits))

            curated_score = compute_escape_bypass_score(targets, upregulated)

            if live_escape_hits and has_rna_data:
                drug["escape_bypass_score"] = round(
                    ew["string_weight"]  * ew["string_hit_score"]
                    + ew["curated_weight"] * curated_score,
                    4,
                )
                drug["escape_note"] = (
                    f"STRING-DB + RNA-confirmed ATRT bypass "
                    f"({len(live_escape_hits)} live hits)"
                )
            elif live_escape_hits and not has_rna_data:
                drug["escape_bypass_score"] = round(
                    0.50 * ew["string_hit_score"]
                    + 0.50 * curated_score,
                    4,
                )
                drug["escape_note"] = "STRING-DB bypass (curated fallback, no RNA)"
            else:
                drug["escape_bypass_score"] = round(curated_score, 4)
                drug["escape_note"] = (
                    "RNA-informed curated bypass" if has_rna_data
                    else "Curated ATRT resistance bypass"
                )

        # ── Step 4: Multi-omic scoring ─────────────────────────────────────────
        candidates = await self._tissue.score_batch(candidates)
        candidates = await self._depmap.score_batch(candidates, disease_name)
        candidates = await self._ppi.score_batch(candidates, disease_data["genes"])

        # ── Step 5: Composite score ────────────────────────────────────────────
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

        # ── Step 6: ATRT-specific drug boosts ──────────────────────────────────
        # EZH2 × 1.40, AURKA × 1.30 (MYC) or × 1.15, SMO × 1.25 (SHH)
        candidates = _apply_atrt_drug_boosts(candidates, subgroup=subgroup)

        # ── Step 7: BBB filter + ATRT location-aware penalty ──────────────────
        candidates, _ = self._bbb_filter.filter_and_rank(
            candidates, apply_penalty=True
        )

        # ATRT BBB penalty — less severe than DIPG brainstem
        bbb_penalties = BBB_CONFIG["dipg_bbb_penalties"].get(
            location, BBB_CONFIG["dipg_bbb_penalties"]["unknown_location"]
        )
        n_penalised = 0
        for c in candidates:
            bbb = c.get("bbb_penetrance", "UNKNOWN")
            if bbb in bbb_penalties:
                c["score_before_bbb"] = c["score"]
                c["score"] = round(c["score"] * bbb_penalties[bbb], 4)
                c["bbb_penalty"] = bbb_penalties[bbb]
                n_penalised += 1

        if n_penalised:
            logger.info(
                "ATRT BBB penalty (%s location): %d candidates penalised",
                location, n_penalised,
            )

        # ── Step 8: Sort and annotate ─────────────────────────────────────────
        sorted_candidates = sorted(
            candidates, key=lambda x: x.get("score", 0), reverse=True
        )

        # Annotate with published ATRT IC50 data
        try:
            from .published_ic50_atrt_validation import (
                annotate_atrt_candidates_with_ic50,
            )
            sorted_candidates = annotate_atrt_candidates_with_ic50(sorted_candidates)
        except Exception as e:
            logger.debug("ATRT IC50 annotation skipped: %s", e)

        # ── Step 9: Generate hypothesis ────────────────────────────────────────
        hypotheses = self._hyp_gen.generate(
            candidates     = sorted_candidates[:top_k],
            cmap_results   = [],
            synergy_combos = [],
            differential_cmap = [],
            genomic_stats  = {
                **smarcb1_stats,
                "has_rna_data": has_rna_data,
            },
            p_value        = smarcb1_stats.get("p_value"),
        )

        # ── Step 10: CMap integration (if pre-computed scores available) ───────
        if self._cmap.has_precomputed_scores():
            for c in sorted_candidates:
                name = c.get("name", "")
                cmap_score = self._cmap.get_precomputed_score(name)
                if cmap_score is not None:
                    c["cmap_score"] = cmap_score
                    c["is_reverser"] = cmap_score > 0.75   # norm_cs < -0.9 → score > 0.75

        return {
            "hypotheses":     hypotheses,
            "top_candidates": sorted_candidates[:top_k],
            "stats": {
                "smarcb1_loss_count":     genomic_stats.get("smarcb1_loss_count", 0),
                "smarca4_loss_count":     genomic_stats.get("smarca4_loss_count", 0),
                "total_samples":          genomic_stats.get("total_samples", 0),
                "upregulated_genes":      len(upregulated),
                "downregulated_genes":    len(downregulated),
                "p_value":                smarcb1_stats.get("p_value"),
                "p_value_label":          smarcb1_stats.get("p_value_label", "N/A"),
                "statistical_note":       smarcb1_stats.get("statistical_note", ""),
                "n_screened":             len(sorted_candidates),
                "subgroup":               subgroup or "pan-ATRT",
                "location":               location,
                "composite_weights":      cw,
                "escape_bypass_mode":     "RNA-confirmed" if has_rna_data else "curated fallback",
                "rna_upregulated_genes":  len(upregulated),
                "n_ezh2_boosted":         sum(
                    1 for c in sorted_candidates if c.get("ezh2_boosted")
                ),
                "n_aurka_boosted":        sum(
                    1 for c in sorted_candidates if c.get("aurka_boosted")
                ),
                "calling_method":         genomic_stats.get("calling_method", "none"),
            },
        }


def _ATRT_FALLBACK_CANDIDATES() -> List[Dict]:
    """
    Minimal fallback drug library when OpenTargets API is unavailable.
    Covers key ATRT therapeutic targets from published biology.
    """
    return [
        {"name": "Tazemetostat",  "targets": ["EZH2"]},
        {"name": "Panobinostat",  "targets": ["HDAC1", "HDAC2", "HDAC3", "HDAC6"]},
        {"name": "Alisertib",     "targets": ["AURKA"]},
        {"name": "Birabresib",    "targets": ["BRD4", "BRD2", "BRD3"]},
        {"name": "Abemaciclib",   "targets": ["CDK4", "CDK6"]},
        {"name": "Marizomib",     "targets": ["PSMB5", "PSMB2", "PSMB1"]},
        {"name": "ONC201",        "targets": ["DRD2", "CLPB"]},
        {"name": "Paxalisib",     "targets": ["PIK3CA", "PIK3CD", "PIK3CG", "MTOR"]},
        {"name": "Vismodegib",    "targets": ["SMO"]},
        {"name": "Vorinostat",    "targets": ["HDAC1", "HDAC2", "HDAC3"]},
    ]