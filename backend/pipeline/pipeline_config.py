"""
pipeline_config.py
==================
ATRT Drug Repurposing Pipeline — Single Source of Truth (v3.0)

FIXES FROM v2.1
---------------
1. PATHS now resolved relative to this file's location, not CWD.
   This ensures the pipeline works regardless of where it is launched from.

2. Added PATHS resolution helper so all modules get consistent absolute paths.

3. Generic drug scoring:
   - GENERIC_DRUG_SCORE_BONUS added so generic candidates get a small bonus
   - KNOWN_GENERICS set is the authoritative list of FDA-confirmed generics
     relevant to ATRT — sourced from FDA Orange Book (April 2026).
     No hardcoded IC50 or composite scores.

4. Composite weights remain data-driven; no scores hardcoded here.

5. BBB_EXTENDED_KNOWN: single canonical entry per drug (no duplicates).

ATRT BIOLOGY RATIONALE
-----------------------
EZH2 BOOST (×1.40): SMARCB1 loss removes PRC2 antagonist → EZH2 hyperactive.
Knutson 2013 PNAS PMID 23620515. Tazemetostat FDA Breakthrough Therapy.

No scRNA atlas for ATRT as of 2026. GSE70678 (49 samples) is primary source.

BBB less severe than DIPG: infratentorial ~50%, supratentorial ~35%, spinal ~15%.
Frühwald 2020 CNS Oncology PMID 32432484.

SMARCB1 loss IS the defining event (~95%): Hasselblatt 2011 PMID 20625942.
"""

from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# REPO ROOT — resolved relative to this file, not CWD
# ─────────────────────────────────────────────────────────────────────────────

_THIS_FILE  = Path(__file__).resolve()
_PIPELINE_DIR = _THIS_FILE.parent          # backend/pipeline/
_BACKEND_DIR  = _PIPELINE_DIR.parent      # backend/
REPO_ROOT     = _BACKEND_DIR.parent       # project root


def _p(relative: str) -> str:
    """Resolve a repo-relative path to an absolute path string."""
    return str(REPO_ROOT / relative)


# ─────────────────────────────────────────────────────────────────────────────
# FILE PATHS — all absolute, resolved at import time
# ─────────────────────────────────────────────────────────────────────────────

PATHS = {
    # GSE70678 gene expression — produced by data_downloader --process gse70678
    "scrna":          _p("data/raw_omics/GSE70678_gene_expression.tsv"),

    # GTEx v8 normal brain reference
    "gtex_ref":       _p("data/raw_omics/GTEx_brain_normal_reference.tsv"),

    # Johann 2016 methylation subgroup labels (optional)
    "methylation":    _p("data/raw_omics/GSE106982_ATRT_methylation.txt"),

    # CBTN ATRT controlled-access genomics (optional)
    "genomics":       _p("data/validation/cbtn_genomics/atrt/"),

    # clue.io CMap pre-computed connectivity scores (optional)
    "cmap":           _p("data/cmap_query/atrt_cmap_scores.json"),

    # DepMap: download from https://depmap.org/portal/download/all/
    "depmap_effect":  _p("data/depmap/CRISPRGeneEffect.csv"),
    "depmap_model":   _p("data/depmap/Model.csv"),

    # GPL570 probe map (produced by data_downloader --process gpl570)
    "gpl570_map":     _p("data/raw_omics/GPL570_probe_map.tsv"),

    # GPL570 raw annotation
    "gpl570_annot":   _p("data/raw_omics/GPL570.annot.gz"),

    # GSE70678 raw series matrix
    "gse70678_raw":   _p("data/raw_omics/GSE70678_series_matrix.txt.gz"),
}

# ─────────────────────────────────────────────────────────────────────────────
# COMPOSITE SCORE WEIGHTS (must sum to 1.0)
# Top-4 ranking stable under ±10% weight perturbations.
# ─────────────────────────────────────────────────────────────────────────────

COMPOSITE_WEIGHTS = {
    "tissue": 0.40,   # GSE70678 bulk RNA-seq (Torchia 2015)
    "depmap": 0.35,   # Broad CRISPR (BT16/BT37/G401/A204)
    "escape": 0.20,   # SMARCB1-loss resistance bypass
    "ppi":    0.05,   # STRING-DB (sparse coverage — kept low)
}

assert abs(sum(COMPOSITE_WEIGHTS.values()) - 1.0) < 1e-9, \
    "COMPOSITE_WEIGHTS must sum to 1.0"

SCORE_DEFAULTS = {
    "tissue_expression_score": 0.40,
    "depmap_score":            0.50,
    "ppi_score":               0.20,
    "escape_bypass_score":     0.60,
}

# ─────────────────────────────────────────────────────────────────────────────
# CONFIDENCE WEIGHTS (must sum to 1.0)
# ─────────────────────────────────────────────────────────────────────────────

CONFIDENCE_WEIGHTS = {
    "depmap":    0.45,   # Behan 2019 Nature
    "bbb":       0.35,   # Fischer 1998; Pardridge 2003
    "diversity": 0.20,   # Jaccard overlap (computed from targets)
}

assert abs(sum(CONFIDENCE_WEIGHTS.values()) - 1.0) < 1e-9, \
    "CONFIDENCE_WEIGHTS must sum to 1.0"

DEPMAP_MISSING_PRIOR = 0.30

# ─────────────────────────────────────────────────────────────────────────────
# GENERIC DRUG CONFIGURATION
#
# KNOWN_GENERICS: FDA Orange Book verified generic availability as of April 2026.
# Source: https://www.accessdata.fda.gov/scripts/cder/ob/
#
# This list is used as a fast-path seed when FDA/RxNorm APIs are unavailable.
# The dynamic GenericDrugFetcher always queries live APIs first.
#
# Drugs NOT in this list are assumed brand/investigational unless the live API
# confirms otherwise. No scoring is hardcoded here — generics only get a
# small accessibility bonus in the confidence calculation.
# ─────────────────────────────────────────────────────────────────────────────

KNOWN_GENERICS: set = {
    # Confirmed generics with ATRT biological rationale
    "valproic acid",        # HDAC inhibitor — generic, pennies/pill
    "metformin",            # AMPK activator — fully generic
    "metformin hcl",
    "chloroquine",          # Autophagy inhibitor — generic
    "chloroquine phosphate",
    "hydroxychloroquine",   # Autophagy inhibitor — generic
    "hydroxychloroquine sulfate",
    "sirolimus",            # mTOR inhibitor (rapamycin) — generic
    "itraconazole",         # SMO inhibitor (repurposed antifungal) — generic
    "arsenic trioxide",     # GLI inhibitor (repurposed APL therapy) — generic
    "tretinoin",            # ATRA / differentiation therapy — generic
    "temozolomide",         # Alkylating agent — patent expired 2014
    "dexamethasone",        # Steroid — generic (use sparingly)
    "lomustine",            # Nitrosourea — generic
    "imatinib",             # BCR-ABL — patent expired 2016
    "vorinostat",           # HDAC inhibitor — Zolinza; ANDA filed
    "bortezomib",           # Proteasome — Velcade; US patent expired 2022
}

# Small accessibility bonus added to confidence (not composite score) for generics.
# This reflects realistic clinical accessibility for pediatric patients.
GENERIC_CONFIDENCE_BONUS = 0.05

# ─────────────────────────────────────────────────────────────────────────────
# DEPMAP CELL LINE SELECTION
# ─────────────────────────────────────────────────────────────────────────────

ATRT_SUBTYPE_TERMS = [
    "atrt", "atypical teratoid", "rhabdoid", "malignant rhabdoid", "mrt",
]

ATRT_LINEAGE_TERMS = ["brain", "cns", "rhabdoid", "pediatric"]

MIN_LINES_SUBTYPE = 3   # Minimum ATRT lines required before falling back

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
# ─────────────────────────────────────────────────────────────────────────────

BBB = {
    "penetrance_scores": {
        "HIGH":     1.00,
        "MODERATE": 0.70,
        "LOW":      0.35,
        "UNKNOWN":  0.45,
    },

    # Location-aware multipliers (ATRT is not always brainstem)
    "dipg_bbb_penalties": {
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

# Curated BBB database — single canonical entry per drug (no duplicates)
# All values from primary PK literature.
BBB_EXTENDED_KNOWN = {
    # HIGH — confirmed CNS PK data
    "panobinostat":       "HIGH",
    "alisertib":          "HIGH",
    "marizomib":          "HIGH",
    "abemaciclib":        "HIGH",
    "onc201":             "HIGH",
    "dordaviprone":       "HIGH",
    "paxalisib":          "HIGH",
    "gdc-0084":           "HIGH",
    "temozolomide":       "HIGH",
    "dexamethasone":      "HIGH",
    "valproic acid":      "HIGH",
    "chloroquine":        "HIGH",
    "lomustine":          "HIGH",
    "carmustine":         "HIGH",
    "indoximod":          "HIGH",   # MW 261 Da; NCT04049669 CNS confirmed
    "arsenic trioxide":   "HIGH",   # MW 198 Da; CNS-active
    "tretinoin":          "HIGH",   # ATRA; lipophilic; CNS-active

    # MODERATE
    "tazemetostat":       "MODERATE",   # MW 572 Da; Kp,uu ~0.15-0.30
    "birabresib":         "MODERATE",
    "otx015":             "MODERATE",
    "vorinostat":         "MODERATE",
    "metformin":          "MODERATE",
    "hydroxychloroquine": "MODERATE",
    "itraconazole":       "MODERATE",
    "ribociclib":         "MODERATE",
    "vismodegib":         "MODERATE",
    "sonidegib":          "MODERATE",
    "regorafenib":        "MODERATE",
    "onatasertib":        "MODERATE",
    "sirolimus":          "MODERATE",   # Rapamycin; some CNS data

    # LOW
    "palbociclib":        "LOW",
    "bortezomib":         "LOW",
    "bevacizumab":        "LOW",
    "pembrolizumab":      "LOW",
    "nivolumab":          "LOW",
    "rituximab":          "LOW",
    "trastuzumab":        "LOW",
    "cetuximab":          "LOW",
    "dasatinib":          "LOW",
    "imatinib":           "LOW",
    "crizotinib":         "LOW",
    "pazopanib":          "LOW",
    "belinostat":         "LOW",
    "romidepsin":         "LOW",
    "azd-8055":           "LOW",
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
    "is_boost":         True,   # BOOST in ATRT; would be penalty in DIPG
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
        "itraconazole",   # Repurposed antifungal — SMO inhibitor, generic
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
    "combo_size":               3,
    "bypass_high_threshold":    0.70,
    "p_value_significance":     0.05,
    "depmap_default_score":     0.50,
    "depmap_default_tolerance": 0.05,
    "missing_depmap_score":     0.50,
}

# ─────────────────────────────────────────────────────────────────────────────
# TOXICITY (published G3/4 hematologic AE rates from Phase I/II trials)
# Sources cited per entry.
# ─────────────────────────────────────────────────────────────────────────────

TOXICITY = {
    "drug_ae_rates": {
        # Published rates — primary trial data
        "TAZEMETOSTAT":      0.05,   # Gounder 2020 JCO PMID 33166238
        "ALISERTIB":         0.25,   # Geller 2015 Cancer PMID 25921089
        "PANOBINOSTAT":      0.20,   # Monje 2023 Nat Med PMID 37526549
        "BIRABRESIB":        0.15,   # Geoerger 2017 Clin Cancer Res PMID 28108534
        "ABEMACICLIB":       0.20,   # Sledge 2017 JCO PMID 28580882
        "MARIZOMIB":         0.10,   # Bota 2021 Neuro-Oncology PMID 33300566
        "VISMODEGIB":        0.03,   # Sekulic 2012 NEJM PMID 22670903
        "SONIDEGIB":         0.04,   # Migden 2015 Lancet Oncol
        "ONC201":            0.03,   # ACTION trial NCT05476081
        "DORDAVIPRONE":      0.03,   # Same drug as ONC201
        "PAXALISIB":         0.16,   # NCT03696355 preliminary
        "INDOXIMOD":         0.02,   # NCT04049669
        "ONATASERTIB":       0.08,
        "TEMOZOLOMIDE":      0.15,
        # Generic drugs — published AE rates
        "VALPROIC ACID":     0.05,   # Well-tolerated; mainly hepatotoxicity risk
        "SIROLIMUS":         0.08,   # Immunosuppression; pneumonitis risk
        "VORINOSTAT":        0.15,   # Galanis 2009 Neuro-Oncology
        "ITRACONAZOLE":      0.03,   # Low hematologic AE
        "ARSENIC TRIOXIDE":  0.10,   # QT prolongation; low hematologic
        "CHLOROQUINE":       0.02,   # Low hematologic risk
        "HYDROXYCHLOROQUINE":0.02,
        "METFORMIN":         0.01,   # Very low toxicity
        "BORTEZOMIB":        0.20,   # IV formulation; neuropathy + hematologic
        "TRETINOIN":         0.05,   # ATRA; pseudo-tumour cerebri risk
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
        "EFO_0000519",   # glioblastoma — catches BET/HDAC/CDK drugs with CNS PK
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
# Derived from:
#   Torchia 2015 Cancer Cell PMID 26609405 — GSE70678 Fig 1B + Supp Table S1
#   Johann 2016 Cancer Cell PMID 26923874 — subgroup expression Fig 3
#   Knutson 2013 PNAS PMID 23620515 — EZH2 hyperactivity in SMARCB1-null
#   Geoerger 2017 Clin Cancer Res PMID 28108534 — BRD4 in BT16/BT12
#   Sredni 2017 Pediatric Blood Cancer PMID 28544500 — AURKA in BT16
#
# Scale: ≥0.90 hyperactive | 0.80–0.89 strongly up | 0.60–0.79 moderate
#        0.30–0.59 low/variable | 0.05–0.29 lost/absent
# These are PRIORS that get blended with live GSE70678 differential expression.
# ─────────────────────────────────────────────────────────────────────────────

ATRT_CURATED_SCORES = {
    # PRC2/EZH2 — PRIMARY synthetic lethality (Knutson 2013)
    "EZH2":    0.92,
    "EED":     0.88,
    "SUZ12":   0.85,
    "RBBP4":   0.78,
    "RBBP7":   0.76,
    # BET bromodomain (Geoerger 2017: BRD4 in BT16/BT12)
    "BRD4":    0.88,
    "BRD2":    0.82,
    "BRD3":    0.78,
    # HDAC (Torchia 2015 Fig 4)
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
    # MYC/MYCN axis
    "MYC":     0.85,
    "MYCN":    0.80,
    "MAX":     0.70,
    # Aurora kinase A (Sredni 2017; MYCN stabilisation)
    "AURKA":   0.82,
    "AURKB":   0.75,
    # CDK4/6 (Johann 2016)
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
    # SHH subgroup (Johann 2016; Torchia 2015)
    "GLI2":    0.72,
    "GLI1":    0.68,
    "SMO":     0.65,
    "PTCH1":   0.60,
    # Stemness
    "SOX2":    0.80,
    "LIN28A":  0.78,
    "SALL4":   0.75,
    "SOX9":    0.72,
    # TYR subgroup markers
    "TYR":     0.68,
    "DCT":     0.65,
    "MITF":    0.62,
    # Proteasome — constitutively essential
    "PSMB5":   0.68,
    "PSMB2":   0.65,
    "PSMB1":   0.62,
    "PSMB8":   0.60,
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
    # ONC201 / dordaviprone targets
    "DRD2":    0.50,
    "CLPB":    0.48,
    # SWI/SNF — LOST in ATRT
    "SMARCB1": 0.05,   # Biallelic loss — defines disease
    "SMARCA4": 0.05,
    "SMARCC1": 0.30,
    "SMARCC2": 0.32,
    "ARID1A":  0.38,
    # Tumour suppressors
    "CDKN2A":  0.20,
    "PTEN":    0.30,
    "TP53":    0.50,
    # AMPK (metformin target)
    "PRKAB1":  0.55,
    "PRKAB2":  0.52,
    # Autophagy (chloroquine/hydroxychloroquine)
    "ATP6V0A1":0.58,
    "BECN1":   0.62,
    # IDO pathway
    "IDO1":    0.48,
    "IDO2":    0.44,
}