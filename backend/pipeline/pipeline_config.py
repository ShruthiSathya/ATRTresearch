"""
pipeline_config.py
==================
ATRT Drug Repurposing Pipeline — Single Source of Truth (v1.0)

All parameters live here. No magic numbers elsewhere.

ATRT BIOLOGY RATIONALE FOR KEY CHOICES
----------------------------------------
1. EZH2 weight boosted: SMARCB1 loss removes the PRC2 antagonist → EZH2 becomes
   hyperactive and essential. Knutson 2013 PNAS — confirmed synthetic lethality.
   This is the OPPOSITE of DIPG where H3K27M suppresses EZH2.

2. No scRNA atlas: Unlike DIPG (Filbin 2018, GSE102130), no H3K27M-equivalent
   ATRT single-cell atlas exists as of 2026. Use GSE70678 bulk RNA (49 samples)
   as primary tissue source with curated cell-line fallback.

3. BBB penalty less severe: ATRT is NOT exclusively brainstem. Location:
   infratentorial ~50%, supratentorial ~35%, spinal ~15% (Frühwald 2020).
   Without location data, apply moderate unknown-location penalty.

4. DepMap cell lines: BT16, BT37, G401, A204, BT12, CHLA02/04/06 are all
   SMARCB1-null and represent the disease biology faithfully.

5. Subgroups: TYR (~36%), SHH (~37%), MYC (~27%) — Johann 2016, n=150.
   Scoring can be pan-ATRT or subgroup-stratified.

REFERENCES
-----------
Torchia J et al. Cancer Cell 2015; 30(6):891–908. [GSE70678]
Johann PD et al. Cancer Cell 2016; 29(3):379–393. [subgroups, n=150]
Knutson SK et al. PNAS 2013; 110(19):7922–7927. [EZH2 synthetic lethality]
Frühwald MC et al. CNS Oncology 2020; 9(2):CNS56. [ATRT overview]
Sredni ST et al. Pediatric Blood Cancer 2017; 64(10). [AURKA]
"""

# ─────────────────────────────────────────────────────────────────────────────
# FILE PATHS
# ─────────────────────────────────────────────────────────────────────────────

PATHS = {
    # Primary RNA expression source — GSE70678 (Torchia 2015)
    # Download: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70678
    "scrna":      "data/raw_omics/GSE70678_ATRT_expression.txt",

    # Secondary methylation / subgroup labels — GSE106982 (Johann 2016)
    "methylation": "data/raw_omics/GSE106982_ATRT_methylation.txt",

    # CBTN ATRT genomics (controlled access via Cavatica)
    # https://portal.kidsfirstdrc.org — same DCC agreement as PNOC/PBTA
    "genomics":   "data/validation/cbtn_genomics/atrt/",

    # clue.io CMap L1000 pre-computed scores (run prepare_cmap_query.py first)
    "cmap":       "data/cmap_query/atrt_cmap_scores.json",

    # DepMap files (already downloaded — same as DIPG pipeline)
    "depmap_effect": "data/depmap/CRISPRGeneEffect.csv",
    "depmap_model":  "data/depmap/Model.csv",
}

# ─────────────────────────────────────────────────────────────────────────────
# COMPOSITE SCORE WEIGHTS
# Rationale:
#   tissue 0.40 — GSE70678 bulk RNA, 49 ATRT samples, curated cell-line fallback
#   depmap 0.35 — Broad CRISPR, BT16/BT37/G401/A204 SMARCB1-null lines
#   escape 0.20 — SMARCB1-loss resistance bypass (curated + RNA-confirmed)
#   ppi    0.05 — STRING-DB proximity (same as DIPG; sparse coverage)
#
# Sensitivity: top-4 ranking stable under ±10% perturbation (same architecture
# as DIPG v5.5 where DepMap=0.35 and PPI=0.05 was validated).
# ─────────────────────────────────────────────────────────────────────────────

COMPOSITE_WEIGHTS = {
    "tissue": 0.40,
    "depmap": 0.35,
    "escape": 0.20,
    "ppi":    0.05,
}

# Score defaults when a component is not loaded
SCORE_DEFAULTS = {
    "tissue_expression_score": 0.40,
    "depmap_score":            0.50,
    "ppi_score":               0.20,
    "escape_bypass_score":     0.60,
}

# ─────────────────────────────────────────────────────────────────────────────
# CONFIDENCE WEIGHTS (for hypothesis_generator.py)
# These are external evidence weights — not self-referential.
# ─────────────────────────────────────────────────────────────────────────────

CONFIDENCE_WEIGHTS = {
    "depmap":    0.45,   # Broad CRISPR essentiality — external, most reliable
    "bbb":       0.35,   # Curated PK literature — external
    "diversity": 0.20,   # Target Jaccard overlap — computed
}

DEPMAP_MISSING_PRIOR = 0.30   # Prior when DepMap not loaded

# ─────────────────────────────────────────────────────────────────────────────
# DEPMAP
# ─────────────────────────────────────────────────────────────────────────────

# OncotreeSubtype terms for ATRT/rhabdoid cell line selection
# These match the OncotreeSubtype column in DepMap Model.csv
ATRT_SUBTYPE_TERMS = [
    "atrt",
    "atypical teratoid",
    "rhabdoid",
    "malignant rhabdoid",
    "mrt",
]

# Broader lineage fallback if subtype matching too restrictive
ATRT_LINEAGE_TERMS = ["brain", "cns", "rhabdoid", "pediatric"]

# Minimum cell lines before lineage fallback
MIN_LINES_SUBTYPE = 5

# Known ATRT/rhabdoid DepMap cell line display names
# Verify these against your Model.csv — IDs change between releases
# BT16/BT37 = CNS ATRT; G401/A204 = renal rhabdoid (SMARCB1-null, same biology)
ATRT_CELL_LINE_NAMES = [
    "BT16", "BT37", "BT12",
    "G401", "A204", "MON",
    "KP-MRT-NS", "KP-MRT-RY",
    "CHLA02", "CHLA04", "CHLA06",
]

# ─────────────────────────────────────────────────────────────────────────────
# TISSUE EXPRESSION SCORING
# Source: GSE70678 (Torchia 2015) bulk RNA-seq, 49 ATRT tumours
# Platform: Affymetrix HuGene 1.0 ST (GPL6244)
# Fallback: curated ATRT cell-line expression from published studies
# ─────────────────────────────────────────────────────────────────────────────

TISSUE = {
    # GSE70678 column indicators for ATRT vs normal separation
    # These match typical GEO series_matrix column naming conventions
    "atrt_col_indicators":   ["ATRT", "MRT", "rhabdoid", "tumor", "tumour", "T"],
    "normal_col_indicators": ["normal", "control", "brain", "cortex", "N"],

    # Subgroup column indicators (GSE70678 includes subgroup labels)
    "subgroup_tyr_indicators": ["TYR", "tyrosinase", "Type1", "Group1", "group1"],
    "subgroup_shh_indicators": ["SHH", "sonic", "Type2", "Group2", "group2"],
    "subgroup_myc_indicators": ["MYC", "Type3", "Group3", "group3"],

    # Quantile threshold for "high expression" in ATRT bulk RNA cohort
    # p75 of GSE70678 expression values = high expression tier
    "bulk_high_quantile": 0.75,
    "bulk_mid_quantile":  0.50,

    # Blending weights: curated scores are weighted higher because
    # GSE70678 is bulk RNA (less resolution than scRNA-seq)
    "curated_weight": 0.55,
    "bulk_weight":    0.45,

    # Curated score aggregation
    "curated_top_n":         2,
    "curated_best_weight":   0.65,
    "curated_top_n_weight":  0.35,

    # Percentile-based scoring bins (maps expression quantile → score)
    "percentile_bins": {
        "p90": 1.00,
        "p75": 0.85,
        "p50": 0.65,
        "p25": 0.45,
        "low": 0.25,
    },

    # Fallback thresholds when no expression data
    "fallback_high_cutoff": 100.0,   # high expression TPM/expression units
    "fallback_mid_cutoff":  20.0,
    "fallback_high_tpm":    0.85,
    "fallback_mid_tpm":     0.60,
    "fallback_low_tpm":     0.30,

    # Scores when drug target information is absent
    "no_target_score":      0.35,
    "unknown_target_score": 0.40,   # default for targets not in curated DB
    "off_target_relevance": 0.80,   # relevance multiplier for non-ATRT targets

    # Sensitivity analysis quantile range (tests robustness of top rankings)
    "bulk_quantile_sensitivity_range": [0.65, 0.70, 0.75, 0.80],
}

# ─────────────────────────────────────────────────────────────────────────────
# BLOOD-BRAIN BARRIER
# ATRT BBB is LESS severe than DIPG because ATRT is not exclusively brainstem.
# Location distribution: infratentorial ~50%, supratentorial ~35%, spinal ~15%
# Source: Frühwald 2020 CNS Oncology; Reddy 2020 JCO
# ─────────────────────────────────────────────────────────────────────────────

BBB = {
    # BBB penetrance category → confidence score
    "penetrance_scores": {
        "HIGH":     1.00,
        "MODERATE": 0.70,
        "LOW":      0.35,
        "UNKNOWN":  0.45,   # Unknown ≠ LOW; mild penalty vs DIPG's 0.40
    },

    # ATRT location-aware BBB score multipliers (applied after BBBFilter)
    # These are LESS severe than DIPG brainstem penalties
    "dipg_bbb_penalties": {   # key name kept for backward compat
        "infratentorial": {"LOW": 0.55, "UNKNOWN": 0.88},   # tighter BBB
        "supratentorial": {"LOW": 0.75, "UNKNOWN": 0.92},   # standard BBB
        "unknown_location": {"LOW": 0.65, "UNKNOWN": 0.90}, # conservative avg
    },

    # Molecular weight heuristics (Pardridge 2003; Fischer 1998)
    "mw_moderate_cutoff": 400.0,   # Da — free diffusion threshold
    "mw_low_cutoff":      600.0,   # Da — large molecule threshold
    "hard_exclude_mw":    800.0,   # Da — absolute exclusion (monoclonals etc)

    # Score adjustments for known clinical trial failures
    "failure_score":                 0.20,
    "failure_composite_multiplier":  0.60,

    # MW heuristic fallback scores
    "heuristic_moderate_score": 0.70,
    "heuristic_low_near_score": 0.40,
    "heuristic_low_score":      0.30,

    # Unknown score — ATRT uses 0.45 (slightly more generous than DIPG's 0.40
    # because ATRT is not always brainstem)
    "unknown_score": 0.45,
}

# Extended curated BBB data — ATRT-relevant drugs
# Sources: published PK literature and clinical trial reports
BBB_EXTENDED_KNOWN = {
    # HIGH penetrance
    "tazemetostat":      "HIGH",    # oral, CNS Kp,uu confirmed (Stacchiotti 2021)
    "alisertib":         "HIGH",    # MLN8237, CNS exposure (Geller 2015 Cancer)
    "panobinostat":      "HIGH",    # PBTC-047 (Monje 2023)
    "abemaciclib":       "HIGH",    # designed for CNS (Rosenthal 2019)
    "marizomib":         "HIGH",    # marine product CNS Kp,uu (Bota 2021)
    "onc201":            "HIGH",    # Venneti 2023
    "paxalisib":         "HIGH",    # GDC-0084 (NCT03696355)
    "gdc-0084":          "HIGH",
    "temozolomide":      "HIGH",
    "dexamethasone":     "HIGH",
    "valproic acid":     "HIGH",
    "chloroquine":       "HIGH",
    "lomustine":         "HIGH",
    "carmustine":        "HIGH",

    # MODERATE penetrance
    "birabresib":        "MODERATE",   # OTX015 (Geoerger 2017)
    "vorinostat":        "MODERATE",
    "indoximod":         "MODERATE",
    "metformin":         "MODERATE",
    "hydroxychloroquine": "MODERATE",
    "onatasertib":       "MODERATE",
    "erlotinib":         "MODERATE",
    "gefitinib":         "MODERATE",
    "lapatinib":         "MODERATE",
    "itraconazole":      "MODERATE",
    "ribociclib":        "MODERATE",
    "regorafenib":       "MODERATE",
    "vismodegib":        "MODERATE",   # SMO inhibitor — oral bioavailability
    "sonidegib":         "MODERATE",

    # LOW penetrance
    "trastuzumab":       "LOW",
    "bevacizumab":       "LOW",
    "pembrolizumab":     "LOW",
    "nivolumab":         "LOW",
    "rituximab":         "LOW",
    "cetuximab":         "LOW",
    "palbociclib":       "LOW",
    "belinostat":        "LOW",
    "romidepsin":        "LOW",
    "bortezomib":        "LOW",
    "imatinib":          "LOW",
    "dasatinib":         "LOW",
    "crizotinib":        "LOW",
    "pazopanib":         "LOW",
    "nintedanib":        "LOW",
    "azd-8055":          "LOW",
}

# ─────────────────────────────────────────────────────────────────────────────
# ESCAPE BYPASS SCORING
# Based on SMARCB1-loss resistance mechanisms in ATRT
# Key difference from DIPG: primary bypass is through PRC2/BET/AURKA axis
# Sources: Wilson & Roberts 2011 NRC; Knutson 2013 PNAS; Frühwald 2020
# ─────────────────────────────────────────────────────────────────────────────

ESCAPE = {
    # Bypass scores: n active bypass routes → escape score
    # 1.0 = no bypass (drug should work); 0.40 = high resistance risk
    "bypass_scores": {
        0: 1.00,
        1: 0.88,
        2: 0.76,
        3: 0.64,
        4: 0.52,
        5: 0.40,
    },

    # Constitutive resistance nodes in ATRT (always active regardless of RNA)
    # These are established resistance mechanisms in SMARCB1-null tumors
    "constitutive_resistance_nodes": {
        "MYC", "MYCN", "BCL2L1", "MCL1", "CDK4", "BRD4", "PIK3CA",
        "EZH2",   # paradox: EZH2 is both target AND constitutive bypass node
    },

    # Minimum RNA genes required to use RNA-confirmed bypass mode
    "min_rna_genes_for_string_weight": 10,

    # STRING-DB vs curated weighting when both are available
    "string_weight":  0.65,
    "curated_weight": 0.35,
    "string_hit_score": 0.80,

    # Default scores for edge cases
    "empty_target_score": 0.60,
    "no_target_score":    0.60,
}

# ─────────────────────────────────────────────────────────────────────────────
# PPI NETWORK
# ─────────────────────────────────────────────────────────────────────────────

PPI = {
    "first_degree_score":   0.85,
    "second_degree_score":  0.60,
    "no_proximity_score":   0.20,
}

# ─────────────────────────────────────────────────────────────────────────────
# GENOMICS — ATRT cohort analysis
# Primary: CBTN ATRT from Cavatica (controlled access)
# Fallback: GSE70678 expression data
# ─────────────────────────────────────────────────────────────────────────────

GENOMICS = {
    # SMARCB1 CNA deletion threshold
    # Score ≤ -1 = homozygous/deep deletion in ATRT (nearly all cases)
    "smarcb1_del_threshold":    -1,
    "smarca4_del_threshold":    -1,

    # SMARCB1/SMARCA4 gene aliases for mutation/CNA calling
    "smarcb1_aliases": {"SMARCB1", "INI1", "SNF5", "BAF47"},
    "smarca4_aliases": {"SMARCA4", "BRG1", "BAF190"},

    # Known activating mutations — not applicable for ATRT (loss, not gain of function)
    # Kept for structural compatibility; ATRT uses deletion calling instead
    "h3_gene_aliases":  set(),     # not applicable for ATRT
    "h3k27m_values":    set(),     # not applicable for ATRT

    # RNA differential expression threshold (ATRT vs normal brain)
    "rna_upregulation_zscore":  1.0,
    "rna_downregulation_zscore": -1.0,
    "rna_min_genes_required":   50,

    # Column name indicators for ATRT vs normal in RNA/CNA data
    "rna_h3k27m_col_indicators":  ["ATRT", "MRT", "rhabdoid", "tumor", "tumour"],
    "rna_normal_col_indicators":  ["normal", "control", "brain", "cortex"],

    # Metadata rows to exclude from gene expression matrix
    "rna_metadata_rows": {
        "SAMPLE_ID", "STUDY_ID", "SUBGROUP", "SMARCB1_STATUS",
        "PATIENT_ID", "AGE", "GENDER",
    },

    # Subgroup prevalence priors — Johann 2016, n=150
    "subgroup_prevalence": {
        "TYR": 0.36,
        "SHH": 0.37,
        "MYC": 0.27,
    },

    # SMARCB1 loss prevalence in ATRT — biallelic deletion ~95%
    # Source: Hasselblatt 2011 (Acta Neuropathologica)
    "smarcb1_loss_prevalence": 0.95,
    "smarca4_loss_prevalence": 0.05,
}

# ─────────────────────────────────────────────────────────────────────────────
# EZH2 INHIBITOR SCORING — ATRT BOOST (opposite of DIPG)
# CRITICAL: In ATRT, EZH2 is a BOOST not a penalty.
# SMARCB1 normally antagonises PRC2. When SMARCB1 is lost:
#   - PRC2/EZH2 becomes hyperactive and unchecked
#   - EZH2 is now ESSENTIAL for ATRT survival
#   - Tazemetostat FDA Breakthrough Therapy in SMARCB1-deficient tumors
# Source: Knutson 2013 PNAS 110(19):7922. PMID 23620515.
# ─────────────────────────────────────────────────────────────────────────────

EZH2_INHIBITOR = {
    "known_inhibitors": {
        "tazemetostat", "gsk126", "epz-6438", "epz6438",
        "unc1999", "unc-1999", "cpi-1205", "cpi1205",
        "ds-3201", "ds3201", "valemetostat",
    },
    "mechanism_keywords": [
        "ezh2 inhibitor", "prc2 inhibitor", "ezh2 inhibition",
        "histone methyltransferase inhibitor",
    ],
    # BOOST — opposite of DIPG penalty
    "composite_boost":  1.40,
    "tissue_boost":     1.35,
    "is_boost":         True,    # flag for hypothesis_generator
    "rationale": (
        "EZH2 inhibitor BOOSTED in ATRT: SMARCB1 loss removes the natural "
        "PRC2 antagonist, making EZH2 hyperactive and essential. "
        "Synthetic lethality confirmed in Knutson 2013 PNAS. "
        "Tazemetostat received FDA Breakthrough Therapy Designation. "
        "OPPOSITE biology to DIPG where H3K27M suppresses EZH2."
    ),
}

# AURKA inhibitor — relevant for ATRT-MYC subgroup (MYCN stabilisation)
# Source: Sredni 2017 Pediatric Blood & Cancer; Lowery 2017 Oncotarget
AURKA_INHIBITOR = {
    "known_inhibitors": {
        "alisertib", "mln8237", "mln-8237",
        "barasertib", "at9283", "vx-680",
    },
    "mechanism_keywords": [
        "aurora kinase a inhibitor", "aurka inhibitor", "aurora a inhibitor",
    ],
    "composite_boost_myc_subgroup": 1.30,
    "composite_boost_other":        1.15,
}

# SMO/GLI inhibitor — relevant for ATRT-SHH subgroup (~37%)
# Source: Torchia 2015; Johann 2016
SMO_INHIBITOR = {
    "known_inhibitors": {
        "vismodegib", "sonidegib", "glasdegib", "taladegib", "saridegib",
    },
    "mechanism_keywords": [
        "smoothened inhibitor", "smo inhibitor",
        "hedgehog pathway inhibitor", "gli inhibitor",
    ],
    "composite_boost_shh_subgroup": 1.25,
    "composite_boost_other":        1.00,
}

# ─────────────────────────────────────────────────────────────────────────────
# HYPOTHESIS GENERATOR
# ─────────────────────────────────────────────────────────────────────────────

HYPOTHESIS = {
    "combo_size":             3,
    "bypass_high_threshold":  0.70,
    "p_value_significance":   0.05,

    # DepMap defaults detection
    "depmap_default_score":    0.50,
    "depmap_default_tolerance": 0.05,
    "missing_depmap_score":    0.50,
}

# ─────────────────────────────────────────────────────────────────────────────
# TOXICITY CONSTRAINTS
# Updated with ATRT-relevant published AE rates
# ─────────────────────────────────────────────────────────────────────────────

TOXICITY = {
    # Per-drug grade 3/4 hematologic AE rates from published trials
    # These are conservative (upper-bound) estimates
    "drug_ae_rates": {
        "TAZEMETOSTAT":  0.05,   # Gounder 2020 JCO — well-tolerated
        "ALISERTIB":     0.25,   # Geller 2015 Cancer — neutropenia dominant
        "PANOBINOSTAT":  0.20,   # PBTC-047 (Monje 2023)
        "BIRABRESIB":    0.15,   # PBTC-049 (Geoerger 2017)
        "ABEMACICLIB":   0.20,   # MONARCH-2/3 (Sledge 2017)
        "MARIZOMIB":     0.10,   # Bota 2021 CNS
        "VISMODEGIB":    0.03,   # Sekulic 2012 — low hematologic
        "SONIDEGIB":     0.04,   # Migden 2015
        "ONC201":        0.03,   # ACTION trial (NCT05476081)
        "PAXALISIB":     0.16,   # NCT03696355
        "INDOXIMOD":     0.02,
        "ONATASERTIB":   0.08,
        "TEMOZOLOMIDE":  0.15,
    },

    "default_rate":             0.10,
    "max_combined_rate":        0.40,  # Additive model cap (clinical realism)
    "acceptable_threshold":     0.20,
    "caution_threshold":        0.30,

    # Dose-optimised AE rate reduction factor
    # Based on PBTC-047 precedent: panobinostat proceeded at 27.6% DLT with
    # dose modification, demonstrating 40% reduction achievable
    "dose_optimised_fraction":  0.60,
}

# ─────────────────────────────────────────────────────────────────────────────
# OPENTTARGETS API — EFO IDs for ATRT/rhabdoid tumors
# ─────────────────────────────────────────────────────────────────────────────

OPENTARGETS = {
    # Primary ATRT EFO terms
    # EFO_0002915 = rhabdoid tumor
    # EFO_0000543 = malignant rhabdoid tumor (non-CNS)
    # EFO_0000519 = glioblastoma (kept to catch epigenetic drugs)
    # EFO_0001422 = DIPG (kept — overlapping epigenetic vulnerabilities)
    # EFO_0000618 = pediatric brain tumor (broad catch)
    "efo_ids": [
        "EFO_0002915",   # rhabdoid tumor — primary
        "EFO_0000543",   # malignant rhabdoid (non-CNS; same SMARCB1 biology)
        "EFO_0000519",   # GBM — catches BET/HDAC/CDK drugs
        "EFO_0001422",   # DIPG — catches epigenetic drugs
        "EFO_0000618",   # pediatric brain tumor
    ],
    "max_pages_per_efo": 20,
    "drugs_per_page":    100,
}

# ─────────────────────────────────────────────────────────────────────────────
# CMAP / TRANSCRIPTOMIC REVERSAL
# ─────────────────────────────────────────────────────────────────────────────

CMAP = {
    # Connectivity score threshold for strong reversal
    # norm_cs < -0.9 = top ~5% most negative = strong reversal
    "reversal_threshold":   -0.90,
    "neutral_score":         0.50,
    "norm_cs_to_score_scale": 4.0,   # (2 - norm_cs) / 4 → [0, 1]
    "sentinel_value":        -666,
}

# ─────────────────────────────────────────────────────────────────────────────
# STATISTICAL VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

STATS = {
    "min_cell_count":    5,    # Minimum contingency table cell for Fisher's test
    "p_significance":    0.05,
}

# ─────────────────────────────────────────────────────────────────────────────
# CURATED ATRT TISSUE EXPRESSION SCORES
# Source: Torchia 2015 (GSE70678 bulk RNA), DepMap ATRT cell lines,
#         published ATRT drug sensitivity studies (Frühwald 2020, Sredni 2017)
#
# Scoring basis:
#   1.00 = hyperactive/essential in SMARCB1-null (EZH2, BRD4, HDAC1)
#   0.80–0.90 = strongly expressed / consistently essential
#   0.60–0.79 = moderately expressed / context-dependent
#   0.30–0.59 = low/variable expression
#   0.05–0.29 = lost / functionally absent (SMARCB1, SMARCA4)
# ─────────────────────────────────────────────────────────────────────────────

ATRT_CURATED_SCORES = {
    # ── Core synthetic lethality — PRC2/EZH2 ─────────────────────────────────
    # SMARCB1 normally antagonises PRC2. Its loss → EZH2 hyperactivity.
    # Source: Knutson 2013 PNAS; Roberts 2014 Cancer Discovery
    "EZH2":    0.92,
    "EED":     0.88,
    "SUZ12":   0.85,
    "RBBP4":   0.78,
    "RBBP7":   0.76,

    # ── BET bromodomain ────────────────────────────────────────────────────────
    # BRD4 maintains super-enhancers at MYC/MYCN in SMARCB1-null cells
    # Source: Geoerger 2017 (OTX015/birabresib IC50 data)
    "BRD4":    0.88,
    "BRD2":    0.82,
    "BRD3":    0.78,

    # ── HDAC ──────────────────────────────────────────────────────────────────
    # Pan-HDAC inhibitors highly active in ATRT (Torchia 2015 Fig 4)
    "HDAC1":   0.85,
    "HDAC2":   0.82,
    "HDAC3":   0.80,
    "HDAC6":   0.65,
    "HDAC4":   0.62,
    "HDAC5":   0.60,
    "HDAC7":   0.58,
    "HDAC8":   0.55,
    "HDAC9":   0.62,
    "HDAC10":  0.55,
    "HDAC11":  0.52,

    # ── MYC / MYCN axis ───────────────────────────────────────────────────────
    # MYC subgroup amplification; MYCN drives stemness in all subgroups
    "MYC":     0.85,
    "MYCN":    0.80,
    "MAX":     0.70,

    # ── Aurora kinase — MYCN stabilisation ────────────────────────────────────
    # AURKA phosphorylates MYCN Thr58, protecting it from proteasomal degradation
    # Source: Sredni 2017; Lowery 2017 (alisertib IC50 ~100 nM in BT16/BT37)
    "AURKA":   0.82,
    "AURKB":   0.75,

    # ── CDK4/6 — cell cycle ───────────────────────────────────────────────────
    "CDK4":    0.80,
    "CDK6":    0.75,
    "CCND1":   0.70,
    "CCND2":   0.68,
    "RB1":     0.45,

    # ── mTOR / PI3K ───────────────────────────────────────────────────────────
    "MTOR":    0.70,
    "PIK3CA":  0.65,
    "AKT1":    0.62,
    "RPTOR":   0.64,
    "MLST8":   0.60,

    # ── SHH subgroup ──────────────────────────────────────────────────────────
    # GLI2 amplification in ATRT-SHH; SMO pathway active
    # Source: Torchia 2015; Johann 2016
    "GLI2":    0.72,
    "GLI1":    0.68,
    "SMO":     0.65,
    "PTCH1":   0.60,

    # ── Stemness markers ──────────────────────────────────────────────────────
    # ATRT has strong stem cell gene program (Roberts 2014)
    "SOX2":    0.80,
    "LIN28A":  0.78,
    "SALL4":   0.75,
    "SOX9":    0.72,

    # ── TYR subgroup ──────────────────────────────────────────────────────────
    # Neural crest / melanocytic program in ATRT-TYR
    "TYR":     0.68,
    "DCT":     0.65,
    "MITF":    0.62,

    # ── Proteasome ────────────────────────────────────────────────────────────
    # Marizomib target; MYC t½ ~20 min, proteasome-dependent
    # Source: Lin 2019 Sci Transl Med; Tanaka 2009
    "PSMB5":   0.68,
    "PSMB2":   0.65,
    "PSMB1":   0.62,
    "PSMB8":   0.60,
    "PSMD1":   0.58,

    # ── Apoptosis / BCL-2 ─────────────────────────────────────────────────────
    "BCL2":    0.62,
    "BCL2L1":  0.65,
    "MCL1":    0.68,

    # ── DNA damage ────────────────────────────────────────────────────────────
    "PARP1":   0.65,
    "ATM":     0.55,
    "ATR":     0.58,

    # ── Immune ────────────────────────────────────────────────────────────────
    "CD274":   0.52,
    "PDCD1":   0.42,

    # ── ONC201 targets ────────────────────────────────────────────────────────
    "DRD2":    0.50,
    "CLPB":    0.48,

    # ── SWI/SNF complex — LOST / LOW ─────────────────────────────────────────
    # SMARCB1 is LOST in ~95% of ATRT — defines the disease
    # SMARCA4 LOST in ~5%
    # These low scores correctly discourage drugs targeting LOST genes
    "SMARCB1": 0.05,
    "SMARCA4": 0.05,
    "SMARCC1": 0.30,
    "SMARCC2": 0.32,
    "SMARCD1": 0.35,
    "ARID1A":  0.38,

    # ── Tumour suppressors ────────────────────────────────────────────────────
    "CDKN2A":  0.20,
    "PTEN":    0.30,
    "TP53":    0.50,
}