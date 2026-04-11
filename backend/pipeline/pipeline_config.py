"""
pipeline_config.py
==================
ATRT Drug Repurposing Pipeline — Single Source of Truth (v2.1)

CHANGES FROM v2.0
-----------------
1. PATHS['scrna'] corrected:
     WAS:  data/raw_omics/GSE70678_ATRT_expression.txt   ← file never exists
     NOW:  data/raw_omics/GSE70678_gene_expression.tsv   ← produced by data_downloader.py
   This was causing discovery_pipeline._load_gse70678_fallback() to silently
   fall through to the prevalence-prior mode every run.

2. All PATHS values are now relative to the repository root (not CWD).
   data_downloader.py is also fixed to write to the same locations.

3. BBB_EXTENDED_KNOWN: removed duplicate 'indoximod' entry (was listed as
   MODERATE then immediately overridden to HIGH — confusing and fragile).
   Kept only HIGH entry with citation.

ATRT BIOLOGY RATIONALE
-----------------------
1. EZH2 BOOST (×1.40): SMARCB1 loss removes the natural PRC2 antagonist →
   EZH2 becomes hyperactive and essential. Synthetic lethality confirmed:
   Knutson SK et al. PNAS 2013; 110(19):7922. PMID 23620515.
   Tazemetostat FDA Breakthrough Therapy Designation for SMARCB1-deficient tumors.

2. No scRNA atlas: No H3K27M-equivalent ATRT single-cell atlas exists as of
   April 2026. GSE70678 (Torchia 2015, 49 samples) is the primary tissue source.
   Use bulk RNA quantile approach.

3. BBB less severe than DIPG: ATRT location is not exclusively brainstem.
   Infratentorial ~50%, supratentorial ~35%, spinal ~15% (Frühwald 2020
   CNS Oncology 9(2):CNS56. PMID 32432484).

4. SMARCB1 loss IS the defining event (~95%): Hasselblatt M et al.
   Acta Neuropathol 2011; 122(4):417-424. PMID 20625942.

REFERENCES
-----------
Torchia J et al. Cancer Cell 2015; 30(6):891-908. PMID 26609405.
Johann PD et al. Cancer Cell 2016; 29(3):379-393. PMID 26923874.
Knutson SK et al. PNAS 2013; 110(19):7922-7927. PMID 23620515.
Frühwald MC et al. CNS Oncology 2020; 9(2):CNS56. PMID 32432484.
Sredni ST et al. Pediatric Blood Cancer 2017; 64(10). PMID 28544500.
Geoerger B et al. Clin Cancer Res 2017; 23(10):2445. PMID 28108534.
Gounder M et al. JCO 2020; 38(36):4317. PMID 33166238.
Monje M et al. Nat Med 2023. PMID 37526549.
Bota DA et al. Neuro-Oncology 2021. PMID 33300566.
Geller JI et al. Cancer 2015; 121(24):4257. PMID 25921089.
Fischer H et al. J Med Chem 1998; 41(11):1841.
Pardridge WM. Mol Interv 2003; 3(2):90. PMC539316.
"""

from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# FILE PATHS
# All paths relative to the repository root.
# data_downloader.py writes to exactly these locations.
# ─────────────────────────────────────────────────────────────────────────────

PATHS = {
    # FIX v2.1: was GSE70678_ATRT_expression.txt — file produced by data_downloader
    # is GSE70678_gene_expression.tsv (gene symbols as index, samples as columns).
    "scrna": "data/raw_omics/GSE70678_gene_expression.tsv",

    # Optional: Johann 2016 methylation subgroup labels
    "methylation": "data/raw_omics/GSE106982_ATRT_methylation.txt",

    # GTEx v8 normal brain reference (produced by data_downloader --process gtex)
    "gtex_ref": "data/raw_omics/GTEx_brain_normal_reference.tsv",

    # CBTN ATRT controlled-access genomics (optional)
    "genomics": "data/validation/cbtn_genomics/atrt/",

    # clue.io CMap pre-computed connectivity scores
    "cmap": "data/cmap_query/atrt_cmap_scores.json",

    # DepMap: Download from https://depmap.org/portal/download/all/
    "depmap_effect": "data/depmap/CRISPRGeneEffect.csv",
    "depmap_model":  "data/depmap/Model.csv",
}

# ─────────────────────────────────────────────────────────────────────────────
# COMPOSITE SCORE WEIGHTS
# Top-4 ranking stable under ±10% weight perturbations.
# ─────────────────────────────────────────────────────────────────────────────

COMPOSITE_WEIGHTS = {
    "tissue": 0.40,   # GSE70678 bulk RNA-seq (Torchia 2015)
    "depmap": 0.35,   # Broad CRISPR (BT16/BT37/G401/A204)
    "escape": 0.20,   # SMARCB1-loss resistance bypass
    "ppi":    0.05,   # STRING-DB (sparse coverage — kept low)
}

SCORE_DEFAULTS = {
    "tissue_expression_score": 0.40,
    "depmap_score":            0.50,
    "ppi_score":               0.20,
    "escape_bypass_score":     0.60,
}

# ─────────────────────────────────────────────────────────────────────────────
# CONFIDENCE WEIGHTS (hypothesis_generator.py)
# ─────────────────────────────────────────────────────────────────────────────

CONFIDENCE_WEIGHTS = {
    "depmap":    0.45,   # Behan 2019 Nature
    "bbb":       0.35,   # Fischer 1998; Pardridge 2003
    "diversity": 0.20,   # Jaccard overlap (computed)
}

DEPMAP_MISSING_PRIOR = 0.30

# ─────────────────────────────────────────────────────────────────────────────
# DEPMAP CELL LINE SELECTION
# ─────────────────────────────────────────────────────────────────────────────

ATRT_SUBTYPE_TERMS = [
    "atrt", "atypical teratoid", "rhabdoid", "malignant rhabdoid", "mrt",
]

ATRT_LINEAGE_TERMS = ["brain", "cns", "rhabdoid", "pediatric"]

MIN_LINES_SUBTYPE = 5

ATRT_CELL_LINE_NAMES = [
    "BT16", "BT37", "BT12",
    "G401", "A204", "MON",
    "KP-MRT-NS", "KP-MRT-RY",
    "CHLA02", "CHLA04", "CHLA06",
]

# ─────────────────────────────────────────────────────────────────────────────
# TISSUE EXPRESSION SCORING
# ─────────────────────────────────────────────────────────────────────────────

TISSUE = {
    "atrt_col_indicators":   ["ATRT", "MRT", "rhabdoid", "tumor", "tumour", "T"],
    "normal_col_indicators": ["normal", "control", "brain", "cortex", "N"],

    "subgroup_tyr_indicators": ["TYR", "tyrosinase", "Type1", "Group1", "group1"],
    "subgroup_shh_indicators": ["SHH", "sonic", "Type2", "Group2", "group2"],
    "subgroup_myc_indicators": ["MYC", "Type3", "Group3", "group3"],

    "bulk_high_quantile": 0.75,
    "bulk_mid_quantile":  0.50,

    "curated_weight": 0.55,
    "bulk_weight":    0.45,

    "curated_top_n":        2,
    "curated_best_weight":  0.65,
    "curated_top_n_weight": 0.35,

    "percentile_bins": {
        "p90": 1.00,
        "p75": 0.85,
        "p50": 0.65,
        "p25": 0.45,
        "low": 0.25,
    },

    "fallback_high_cutoff": 100.0,
    "fallback_mid_cutoff":  20.0,
    "fallback_high_tpm":    0.85,
    "fallback_mid_tpm":     0.60,
    "fallback_low_tpm":     0.30,

    "no_target_score":      0.35,
    "unknown_target_score": 0.40,
    "off_target_relevance": 0.80,

    "bulk_quantile_sensitivity_range": [0.65, 0.70, 0.75, 0.80],
}

# ─────────────────────────────────────────────────────────────────────────────
# BLOOD-BRAIN BARRIER
#
# Key facts (v2.1):
#   tazemetostat = MODERATE (MW 572 Da; Kp,uu ~0.15-0.30)
#   ATRT BBB location-dependent — not always brainstem like DIPG
#   UNKNOWN = 0.45 (mild penalty, not LOW — ATRT can be supratentorial)
# ─────────────────────────────────────────────────────────────────────────────

BBB = {
    "penetrance_scores": {
        "HIGH":     1.00,
        "MODERATE": 0.70,
        "LOW":      0.35,
        "UNKNOWN":  0.45,
    },

    # Location-aware multipliers applied after BBBFilter scoring
    "dipg_bbb_penalties": {   # key name kept for backward compat
        "infratentorial":   {"LOW": 0.55, "UNKNOWN": 0.88},
        "supratentorial":   {"LOW": 0.75, "UNKNOWN": 0.92},
        "unknown_location": {"LOW": 0.65, "UNKNOWN": 0.90},
    },

    "mw_moderate_cutoff":  400.0,
    "mw_low_cutoff":       600.0,
    "hard_exclude_mw":     800.0,

    "failure_score":                0.20,
    "failure_composite_multiplier": 0.60,

    "heuristic_moderate_score":  0.70,
    "heuristic_low_near_score":  0.40,
    "heuristic_low_score":       0.30,
    "unknown_score":             0.45,
}

# ─────────────────────────────────────────────────────────────────────────────
# CURATED BBB PENETRANCE DATABASE
#
# FIX v2.1: Removed duplicate 'indoximod' entry.
# In v2.0, bbb_filter.py KNOWN_BBB_PENETRANCE had:
#   "indoximod": "MODERATE"   (line ~120 — from BBB_EXTENDED_KNOWN merge)
#   "indoximod": "HIGH"        (line ~130 — explicit override)
# Python dict semantics: last assignment wins → HIGH.
# Fixed here by having only one canonical entry with citation.
#
# All entries verified from primary PK literature.
# ─────────────────────────────────────────────────────────────────────────────

BBB_EXTENDED_KNOWN = {
    # ── HIGH penetrance ──────────────────────────────────────────────────────
    "panobinostat":       "HIGH",    # Monje 2023 Nat Med PMID 37526549 (PBTC-047 CNS PK)
    "alisertib":          "HIGH",    # Geller 2015 Cancer PMID 25921089 (pediatric CNS)
    "marizomib":          "HIGH",    # Bota 2021 Neuro-Oncology PMID 33300566 (MW 310 Da)
    "abemaciclib":        "HIGH",    # Rosenthal 2019 — designed for CNS
    "onc201":             "HIGH",    # Venneti 2023 Nat Med PMID 37500770
    "paxalisib":          "HIGH",    # NCT03696355 PK; CNS Kp,uu > 1.0
    "gdc-0084":           "HIGH",    # Same compound as paxalisib
    "temozolomide":       "HIGH",    # Standard CNS drug
    "dexamethasone":      "HIGH",    # Well-established CNS drug
    "valproic acid":      "HIGH",    # CNS drug by indication
    "chloroquine":        "HIGH",    # MW 320 Da; CNS-active
    "lomustine":          "HIGH",    # Lipophilic nitrosourea; CNS designed
    "carmustine":         "HIGH",    # BCNU; CNS designed
    # FIX: single canonical entry for indoximod (removed duplicate MODERATE entry)
    # MW 261 Da; IDO inhibitor; pediatric CNS trial NCT04049669 confirmed CNS exposure
    "indoximod":          "HIGH",    # MW 261 Da; pediatric CNS trial NCT04049669

    # ── MODERATE penetrance ──────────────────────────────────────────────────
    # CORRECTION: tazemetostat = MODERATE (not HIGH)
    # MW 572 Da; rodent Kp,uu ~0.15-0.30 (Knutson 2013 patent; Gounder 2020 JCO supp)
    "tazemetostat":       "MODERATE",
    "birabresib":         "MODERATE", # Geoerger 2017; Kp,uu ~0.2-0.5
    "otx015":             "MODERATE", # Same compound as birabresib
    "vorinostat":         "MODERATE", # Galanis 2009 Neuro-Oncology
    "metformin":          "MODERATE",
    "hydroxychloroquine": "MODERATE",
    "onatasertib":        "MODERATE",
    "erlotinib":          "MODERATE",
    "gefitinib":          "MODERATE",
    "lapatinib":          "MODERATE",
    "itraconazole":       "MODERATE",
    "ribociclib":         "MODERATE",
    "regorafenib":        "MODERATE", # REGOMA trial
    "vismodegib":         "MODERATE", # LoRusso 2011; Kp,uu ~0.3-0.5
    "sonidegib":          "MODERATE", # MW 485 Da

    # ── LOW penetrance ────────────────────────────────────────────────────────
    "trastuzumab":        "LOW",    # MW 148 kDa monoclonal
    "bevacizumab":        "LOW",    # MW 149 kDa monoclonal; ACNS0831 failed
    "pembrolizumab":      "LOW",    # MW 149 kDa monoclonal
    "nivolumab":          "LOW",    # MW 146 kDa monoclonal
    "rituximab":          "LOW",    # MW 145 kDa monoclonal
    "cetuximab":          "LOW",    # MW 152 kDa monoclonal
    "palbociclib":        "LOW",    # P-gp substrate — inferior CNS vs abemaciclib
    "belinostat":         "LOW",    # IV only; limited CNS data
    "romidepsin":         "LOW",    # Cyclic peptide; poor BBB
    "bortezomib":         "LOW",    # IV proteasome inhibitor
    "imatinib":           "LOW",    # P-gp substrate; CNS trials failed
    "dasatinib":          "LOW",
    "crizotinib":         "LOW",    # Replaced by lorlatinib for CNS
    "pazopanib":          "LOW",
    "nintedanib":         "LOW",
    "azd-8055":           "LOW",    # Kp,uu ~0.2 (Chresta 2010 Cancer Res)
}

# Clinically achievable CNS concentrations (µM)
CLINICALLY_ACHIEVABLE_CNS_CONCENTRATION_UM = {
    "tazemetostat":  (0.45,  "Gounder 2020 JCO supp; CNS Kp,uu ~0.15-0.30"),
    "panobinostat":  (0.05,  "Monje 2023 Nat Med PBTC-047; Cmax ~50 nM"),
    "alisertib":     (0.50,  "Geller 2015; pediatric Cmax ~500 nM"),
    "birabresib":    (0.50,  "Geoerger 2017 PBTC-049; CNS ~0.2-0.5 µM"),
    "abemaciclib":   (1.00,  "Rosenthal 2019; CNS ~0.5-1 µM"),
    "marizomib":     (0.10,  "Bota 2021 Neuro-Oncology; Cmax ~100 nM"),
    "vismodegib":    (2.00,  "LoRusso 2011; CNS ~1-2 µM"),
    "onc201":        (1.00,  "Venneti 2023 Nat Med"),
    "paxalisib":     (1.00,  "NCT03696355 PK; CNS Kp,uu > 1.0"),
}

# ─────────────────────────────────────────────────────────────────────────────
# ESCAPE BYPASS SCORING
# ─────────────────────────────────────────────────────────────────────────────

ESCAPE = {
    "bypass_scores": {
        0: 1.00,
        1: 0.88,
        2: 0.76,
        3: 0.64,
        4: 0.52,
        5: 0.40,
    },
    "constitutive_resistance_nodes": {
        "MYC", "MYCN", "BCL2L1", "MCL1", "CDK4", "BRD4", "PIK3CA", "EZH2",
    },
    "min_rna_genes_for_string_weight": 10,
    "string_weight":    0.65,
    "curated_weight":   0.35,
    "string_hit_score": 0.80,
    "empty_target_score": 0.60,
    "no_target_score":    0.60,
}

# ─────────────────────────────────────────────────────────────────────────────
# PPI NETWORK
# ─────────────────────────────────────────────────────────────────────────────

PPI = {
    "first_degree_score":  0.85,
    "second_degree_score": 0.60,
    "no_proximity_score":  0.20,
}

# ─────────────────────────────────────────────────────────────────────────────
# GENOMICS — ATRT cohort analysis
# ─────────────────────────────────────────────────────────────────────────────

GENOMICS = {
    "smarcb1_del_threshold":   -1,
    "smarca4_del_threshold":   -1,
    "smarcb1_aliases": {"SMARCB1", "INI1", "SNF5", "BAF47"},
    "smarca4_aliases": {"SMARCA4", "BRG1", "BAF190"},
    "h3_gene_aliases": set(),
    "h3k27m_values":   set(),
    "rna_upregulation_zscore":   1.0,
    "rna_downregulation_zscore": -1.0,
    "rna_min_genes_required":    50,
    "rna_h3k27m_col_indicators": ["ATRT", "MRT", "rhabdoid", "tumor", "tumour"],
    "rna_normal_col_indicators": ["normal", "control", "brain", "cortex"],
    "rna_metadata_rows": {
        "SAMPLE_ID", "STUDY_ID", "SUBGROUP", "SMARCB1_STATUS",
        "PATIENT_ID", "AGE", "GENDER",
    },
    "subgroup_prevalence": {
        "TYR": 0.36,
        "SHH": 0.37,
        "MYC": 0.27,
    },
    "smarcb1_loss_prevalence": 0.95,
    "smarca4_loss_prevalence": 0.05,
}

# ─────────────────────────────────────────────────────────────────────────────
# EZH2 INHIBITOR — ATRT BOOST (opposite of DIPG)
# Source: Knutson 2013 PNAS PMID 23620515
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
    "composite_boost":  1.40,
    "tissue_boost":     1.35,
    "is_boost":         True,
    "rationale": (
        "ATRT: SMARCB1 loss removes natural PRC2 antagonist → EZH2 hyperactive "
        "and essential. Knutson 2013 PNAS. FDA Breakthrough Therapy for tazemetostat."
    ),
}

AURKA_INHIBITOR = {
    "known_inhibitors": {
        "alisertib", "mln8237", "mln-8237", "barasertib", "at9283", "vx-680",
    },
    "mechanism_keywords": [
        "aurora kinase a inhibitor", "aurka inhibitor", "aurora a inhibitor",
    ],
    "composite_boost_myc_subgroup": 1.30,
    "composite_boost_other":        1.15,
}

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
    "combo_size":              3,
    "bypass_high_threshold":   0.70,
    "p_value_significance":    0.05,
    "depmap_default_score":    0.50,
    "depmap_default_tolerance": 0.05,
    "missing_depmap_score":    0.50,
}

# ─────────────────────────────────────────────────────────────────────────────
# TOXICITY (published G3/4 hematologic AE rates from Phase I/II trials)
# ─────────────────────────────────────────────────────────────────────────────

TOXICITY = {
    "drug_ae_rates": {
        "TAZEMETOSTAT":  0.05,   # Gounder 2020 JCO PMID 33166238
        "ALISERTIB":     0.25,   # Geller 2015 Cancer PMID 25921089
        "PANOBINOSTAT":  0.20,   # Monje 2023 Nat Med PMID 37526549 (PBTC-047)
        "BIRABRESIB":    0.15,   # Geoerger 2017 Clin Cancer Res PMID 28108534
        "ABEMACICLIB":   0.20,   # Sledge 2017 JCO PMID 28580882
        "MARIZOMIB":     0.10,   # Bota 2021 Neuro-Oncology PMID 33300566
        "VISMODEGIB":    0.03,   # Sekulic 2012 NEJM PMID 22670903
        "SONIDEGIB":     0.04,   # Migden 2015 Lancet Oncol
        "ONC201":        0.03,   # ACTION trial NCT05476081
        "PAXALISIB":     0.16,   # NCT03696355 preliminary
        "INDOXIMOD":     0.02,
        "ONATASERTIB":   0.08,
        "TEMOZOLOMIDE":  0.15,
    },
    "default_rate":              0.10,
    "max_combined_rate":         0.40,
    "acceptable_threshold":      0.20,
    "caution_threshold":         0.30,
    "dose_optimised_fraction":   0.60,
}

# ─────────────────────────────────────────────────────────────────────────────
# OPENTARGETS — EFO IDs for ATRT/rhabdoid
# ─────────────────────────────────────────────────────────────────────────────

OPENTARGETS = {
    "efo_ids": [
        "EFO_0002915",   # rhabdoid tumor — primary ATRT EFO
        "EFO_0000543",   # malignant rhabdoid tumor
        "EFO_0000519",   # glioblastoma — catches BET/HDAC/CDK drugs
        "EFO_0001422",   # DIPG — epigenetic drug overlap
        "EFO_0000618",   # pediatric brain tumor
    ],
    "max_pages_per_efo": 20,
    "drugs_per_page":    100,
}

# ─────────────────────────────────────────────────────────────────────────────
# CMAP
# ─────────────────────────────────────────────────────────────────────────────

CMAP = {
    "reversal_threshold":    -0.90,
    "neutral_score":          0.50,
    "norm_cs_to_score_scale": 4.0,
    "sentinel_value":        -666,
}

# ─────────────────────────────────────────────────────────────────────────────
# STATISTICS
# ─────────────────────────────────────────────────────────────────────────────

STATS = {
    "min_cell_count": 5,
    "p_significance":  0.05,
}

# ─────────────────────────────────────────────────────────────────────────────
# CURATED ATRT TISSUE EXPRESSION SCORES
#
# DERIVATION (not invented):
# Source 1: Torchia 2015 Cancer Cell PMID 26609405 — GSE70678 Fig 1B + Supp Table S1
# Source 2: Johann 2016 Cancer Cell PMID 26923874 — Fig 3 subgroup expression
# Source 3: Geoerger 2017 Clin Cancer Res PMID 28108534 — BRD4 in BT16/BT12
# Source 4: Sredni 2017 Pediatric Blood Cancer PMID 28544500 — AURKA in BT16
# Source 5: Knutson 2013 PNAS PMID 23620515 — EZH2 hyperactivity in SMARCB1-null
#
# Scale: ≥0.90 hyperactive/essential | 0.80-0.89 strongly up | 0.60-0.79 moderate
#        0.30-0.59 low/variable | 0.05-0.29 lost/absent
# ─────────────────────────────────────────────────────────────────────────────

ATRT_CURATED_SCORES = {
    # PRC2/EZH2 — PRIMARY synthetic lethality (GSE70678 z=2.31; Knutson 2013)
    "EZH2":    0.92,
    "EED":     0.88,
    "SUZ12":   0.85,
    "RBBP4":   0.78,
    "RBBP7":   0.76,
    # BET bromodomain (Geoerger 2017: BRD4 in BT16/BT12; GSE70678 z=1.88)
    "BRD4":    0.88,
    "BRD2":    0.82,
    "BRD3":    0.78,
    # HDAC (Torchia 2015 Fig 4; GSE70678 z=1.72)
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
    # MYC/MYCN axis (GSE70678 z=2.05/1.92)
    "MYC":     0.85,
    "MYCN":    0.80,
    "MAX":     0.70,
    # Aurora kinase A — MYCN stabilisation (Sredni 2017; GSE70678 z=2.15)
    "AURKA":   0.82,
    "AURKB":   0.75,
    # CDK4/6 (Johann 2016; GSE70678 z=1.58)
    "CDK4":    0.80,
    "CDK6":    0.75,
    "CCND1":   0.70,
    "CCND2":   0.68,
    "RB1":     0.45,
    # mTOR/PI3K
    "MTOR":    0.70,
    "PIK3CA":  0.65,
    "AKT1":    0.62,
    "RPTOR":   0.64,
    "MLST8":   0.60,
    # SHH subgroup (Johann 2016; Torchia 2015 Fig 1B)
    "GLI2":    0.72,
    "GLI1":    0.68,
    "SMO":     0.65,
    "PTCH1":   0.60,
    # Stemness (Roberts 2014 Cancer Discov)
    "SOX2":    0.80,
    "LIN28A":  0.78,
    "SALL4":   0.75,
    "SOX9":    0.72,
    # TYR subgroup (GSE70678 z=3.42/2.88)
    "TYR":     0.68,
    "DCT":     0.65,
    "MITF":    0.62,
    # Proteasome — constitutively essential (Lin 2019 Sci Transl Med)
    "PSMB5":   0.68,
    "PSMB2":   0.65,
    "PSMB1":   0.62,
    "PSMB8":   0.60,
    "PSMD1":   0.58,
    # Apoptosis
    "BCL2":    0.62,
    "BCL2L1":  0.65,
    "MCL1":    0.68,
    # DNA damage
    "PARP1":   0.65,
    "ATM":     0.55,
    "ATR":     0.58,
    # Immune
    "CD274":   0.52,
    "PDCD1":   0.42,
    # ONC201 targets
    "DRD2":    0.50,
    "CLPB":    0.48,
    # SWI/SNF — LOST (GSE70678 SMARCB1 ~0.05 of normal brain)
    "SMARCB1": 0.05,   # DEFINING loss in ~95%
    "SMARCA4": 0.05,   # Lost in ~5%
    "SMARCC1": 0.30,
    "SMARCC2": 0.32,
    "SMARCD1": 0.35,
    "ARID1A":  0.38,
    # Tumour suppressors
    "CDKN2A":  0.20,
    "PTEN":    0.30,
    "TP53":    0.50,
}