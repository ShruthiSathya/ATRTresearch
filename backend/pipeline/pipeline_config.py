from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# DATA PATHS
# ─────────────────────────────────────────────────────────────────────────────

PATHS = {
    "scrna":    "data/raw_omics/GSE102130_K27Mproject.RSEM.vh20170621.txt",
    "genomics": "data/validation/cbtn_genomics/",
    "results":  "results/",
}

# ─────────────────────────────────────────────────────────────────────────────
# COMPOSITE SCORE WEIGHTS  (must sum to 1.0)
# ─────────────────────────────────────────────────────────────────────────────

COMPOSITE_WEIGHTS = {
    # tissue 0.40 — disease-specific GSC context (Filbin 2018, GSE102130)
    # depmap 0.35 — externally validated functional dependency
    #               (Behan et al. 2019 Nature; Tsherniak et al. 2017 Cell)
    #               INCREASED from 0.30: PPI is too sparse (490/557 drugs score 0.20)
    #               so redistributing its contribution to the strongest signal.
    # escape 0.20 — computed resistance estimate, informative but indirect
    # ppi    0.05 — reduced from 0.10: only 67/557 drugs get non-floor PPI scores
    #               due to curated neighbor coverage gaps. At 0.10 this was
    #               contributing near-zero discriminative power for 88% of drugs.
    #               Sensitivity analysis: top-2 stable under this change.
    "tissue":  0.40,
    "depmap":  0.35,
    "escape":  0.20,
    "ppi":     0.05,
}

assert abs(sum(COMPOSITE_WEIGHTS.values()) - 1.0) < 1e-9, \
    "COMPOSITE_WEIGHTS must sum to 1.0"

# Default values when a score key is absent from a candidate dict
SCORE_DEFAULTS = {
    "tissue_expression_score": 0.10,
    "depmap_score":            0.10,
    "ppi_score":               0.10,
    "escape_bypass_score":     0.40,
}

# ─────────────────────────────────────────────────────────────────────────────
# CONFIDENCE SCORE WEIGHTS  (must sum to 1.0)
# ─────────────────────────────────────────────────────────────────────────────

CONFIDENCE_WEIGHTS = {
    "depmap":    0.45,
    "bbb":       0.35,
    "diversity": 0.20,
}

assert abs(sum(CONFIDENCE_WEIGHTS.values()) - 1.0) < 1e-9, \
    "CONFIDENCE_WEIGHTS must sum to 1.0"

DEPMAP_MISSING_PRIOR = 0.30

# ─────────────────────────────────────────────────────────────────────────────
# TISSUE EXPRESSION SCORING
# ─────────────────────────────────────────────────────────────────────────────

TISSUE = {
    "curated_weight": 0.50,
    "sc_weight":      0.50,

    # GSC stem-cell quantile — p85 is the primary value.
    # Sensitivity analysis runs p75/p80/p85/p90 and checks rank stability.
    # See tissue_expression.py: TissueExpressionScorer.run_quantile_sensitivity()
    "gsc_stem_quantile": 0.85,

    # Range for sensitivity analysis
    "gsc_quantile_sensitivity_range": [0.75, 0.80, 0.85, 0.90],

    "stem_markers": ["SOX2", "NES", "PROM1", "CD44"],

    "percentile_bins": {
        "p90": 1.00,
        "p75": 0.82,
        "p50": 0.62,
        "p25": 0.42,
        "low": 0.18,
    },

    "fallback_high_tpm":    1.0,
    "fallback_mid_tpm":     0.6,
    "fallback_low_tpm":     0.2,
    "fallback_high_cutoff": 2.0,
    "fallback_mid_cutoff":  0.5,

    "off_target_relevance": 0.35,
    "unknown_target_score": 0.40,

    "curated_best_weight":  0.60,
    "curated_top_n":        3,
    "curated_top_n_weight": 0.40,

    "no_target_score": 0.40,
}

assert abs(TISSUE["curated_weight"] + TISSUE["sc_weight"] - 1.0) < 1e-9, \
    "TISSUE curated_weight + sc_weight must sum to 1.0"

# ─────────────────────────────────────────────────────────────────────────────
# EZH2 INHIBITOR PENALTY
# H3K27M is a dominant-negative inhibitor of PRC2/EZH2. Residual EZH2
# activity is LOW in H3K27M DIPG (Bender et al. 2014, Cancer Cell).
# EZH2 inhibitors (tazemetostat, GSK126) are therefore non-rational in
# H3K27M-mutant DIPG — they inhibit an already-inhibited enzyme.
# This penalty is applied in dipg_specialization.py.
# ─────────────────────────────────────────────────────────────────────────────

EZH2_INHIBITOR = {
    # Drugs known to be EZH2 inhibitors — penalised in H3K27M DIPG
    "known_inhibitors": {
        "tazemetostat", "gsk126", "epz-6438", "epz6438",
        "unc1999", "unc-1999", "cpi-1205", "cpi1205",
        "ds-3201", "ds3201", "valemetostat",
    },
    # Keywords in mechanism string that identify EZH2 inhibitors
    "mechanism_keywords": [
        "ezh2 inhibitor", "prc2 inhibitor", "ezh2 inhibition",
        "histone methyltransferase inhibitor",
    ],
    # Score penalty: multiply tissue and composite score by this factor
    "composite_penalty": 0.50,
    # Rationale logged when penalty is applied
    "rationale": (
        "EZH2 inhibitor penalised in H3K27M DIPG: H3K27M dominant-negatively "
        "inhibits PRC2/EZH2 (Bender et al. 2014). Residual EZH2 activity is "
        "already suppressed — EZH2 inhibitors have no additional effect. "
        "Tazemetostat FDA-approved for EZH2-mutant FL but contraindicated "
        "mechanistically in H3K27M DIPG."
    ),
}

# ─────────────────────────────────────────────────────────────────────────────
# BBB KNOWN PENETRANCE — AZD-8055 and other UNKNOWN drugs
# Fills gaps in bbb_filter.py curated list using published PK data.
# ─────────────────────────────────────────────────────────────────────────────

BBB_EXTENDED_KNOWN = {
    # AZD-8055: mTOR kinase inhibitor. CNS Kp,uu ~0.2 in preclinical rodent PK
    # (Chresta et al. 2010, Mol Cancer Ther). Low but not negligible penetrance.
    "azd-8055":    "LOW",
    # CC-115: dual mTOR/DNA-PK inhibitor. Phase I CNS data limited; MW ~370 Da
    # suggests moderate heuristic, but no confirmed CNS exposure data.
    "cc-115":      "MODERATE",
    # Crizotinib: ALK/MET inhibitor. CNS penetrance is LOW (designed before
    # lorlatinib). Replaced by lorlatinib for CNS disease. MW=450 Da but P-gp substrate.
    # Source: Shaw et al. 2013 NEJM; Costa et al. 2015 J Thorac Oncol.
    "crizotinib":  "LOW",
    # Crenolanib: FLT3/PDGFRA inhibitor. MW=519 Da; limited CNS PK data.
    "crenolanib":  "LOW",
    # Dovitinib: multi-kinase. MW=392 Da heuristic → MODERATE, but no CNS data.
    "dovitinib":   "MODERATE",
    # Vatalanib: VEGFR inhibitor. MW=347 Da; some CNS exposure data in GBM trials.
    "vatalanib":   "MODERATE",
    # Nintedanib: FGFR/VEGFR/PDGFR. MW=540 Da; P-gp substrate; LOW CNS penetrance.
    "nintedanib":  "LOW",
    # Regorafenib: multi-kinase. MW=483 Da; some CNS data from REGOMA trial.
    "regorafenib": "MODERATE",
    # Pazopanib: VEGFR/PDGFR/KIT. MW=473 Da; poor CNS penetrance (P-gp substrate).
    "pazopanib":   "LOW",
    # Sunitinib: multi-kinase. MW=399 Da; some CNS exposure but inconsistent.
    "sunitinib":   "LOW",
    # Tovetumab (Olaratumab): anti-PDGFRA monoclonal antibody. MW ~148 kDa → LOW.
    "tovetumab":   "LOW",
    "olaratumab":  "LOW",
    # Depatuxizumab mafodotin: anti-EGFR ADC. MW ~152 kDa → LOW.
    "depatuxizumab mafodotin": "LOW",
    "depatuxizumab": "LOW",
}

# ─────────────────────────────────────────────────────────────────────────────
# ESCAPE BYPASS SCORING
# ─────────────────────────────────────────────────────────────────────────────

ESCAPE = {
    "bypass_scores": {
        0: 1.00,
        1: 0.85,
        2: 0.72,
        3: 0.58,
        4: 0.40,
    },

    "no_target_score":    0.75,
    "empty_target_score": 0.70,

    # Data-driven escape: weight live RNA-upregulated bypass hits more than curated
    # FIX: was 0.60/0.40 split. Now 0.70/0.30 — live RNA data from GSE115397
    # is more specific than the hand-written bypass map.
    "string_weight":  0.70,
    "curated_weight": 0.30,

    "string_hit_score": 0.40,

    # Constitutive resistance nodes — hardwired biology regardless of RNA expression
    # RAD51/BRCA1: H3K27M DIPG is HR-proficient — PARP-i have no synthetic lethality
    "constitutive_resistance_nodes": {
        "PIK3CA", "MTOR", "MYC", "MYCN", "BCL2L1", "MCL1", "CDK4", "BRD4",
        "RAD51", "BRCA1",
    },

    # Minimum upregulated gene hits from RNA data before we weight string_weight
    # If RNA data produced 0 upregulated genes (file not loaded), fall back to
    # curated-only scoring with equal weights
    "min_rna_genes_for_string_weight": 10,
}

# ─────────────────────────────────────────────────────────────────────────────
# TOXICITY CONSTRAINT
# ─────────────────────────────────────────────────────────────────────────────

TOXICITY = {
    "max_combined_rate": 0.40,
    "default_rate":      0.05,

    # Thresholds calibrated for pediatric DIPG (PBTC-047/049 precedent)
    "acceptable_threshold": 0.20,
    "caution_threshold":    0.33,

    # Toxicity confidence interval model:
    # The additive AE rate model is CONSERVATIVE (assumes independence).
    # Real combination toxicity with dose reduction is typically lower.
    # Report confidence interval: [multiplier_conservative, multiplier_optimistic]
    # where optimistic assumes 60% of additive rate (dose-optimised).
    "dose_optimised_fraction": 0.60,
}

# ─────────────────────────────────────────────────────────────────────────────
# BBB FILTER
# ─────────────────────────────────────────────────────────────────────────────

BBB = {
    "hard_exclude_mw": 800.0,
    "mw_moderate_cutoff": 400.0,
    "mw_low_cutoff":      600.0,

    "penetrance_scores": {
        "HIGH":    1.0,
        "MODERATE": 0.6,
        "LOW":     0.2,
        "UNKNOWN": 0.4,   # FIX: was 0.5 (same as MODERATE). Unknown ≠ MODERATE.
                          # 0.4 penalises lack of data without assuming worst case.
    },

    "failure_score": 0.05,
    "failure_composite_multiplier": 0.10,

    "heuristic_moderate_score": 0.6,
    "heuristic_low_near_score": 0.3,
    "heuristic_low_score":      0.1,
    "unknown_score":            0.4,  # matches penetrance_scores["UNKNOWN"]
}

# ─────────────────────────────────────────────────────────────────────────────
# PPI NETWORK SCORES
# ─────────────────────────────────────────────────────────────────────────────

PPI = {
    "first_degree_score":  0.85,
    "second_degree_score": 0.60,
    "no_proximity_score":  0.20,
}

# ─────────────────────────────────────────────────────────────────────────────
# HYPOTHESIS GENERATOR THRESHOLDS
# ─────────────────────────────────────────────────────────────────────────────

HYPOTHESIS = {
    "p_value_significance":     0.05,
    "bypass_high_threshold":    0.50,
    "depmap_default_score":     0.50,
    "depmap_default_tolerance": 0.01,
    "depmap_missing_score":     0.10,
    "combo_size":               3,
    "missing_depmap_score":     0.10,
}

# ─────────────────────────────────────────────────────────────────────────────
# GENOMIC VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

GENOMICS = {
    "rna_upregulation_zscore": 1.0,
    "cna_deletion_threshold":  -1,
    "rna_min_genes_required":  50,

    "h3k27m_values":   {"K28M", "K27M"},
    "h3_gene_aliases": {"H3-3A", "H3F3A", "HIST1H3B", "H33A", "H3.3A"},

    "rna_h3k27m_col_indicators": ["Pons", "K27M", "DIPG", "H3K27M", "pontine"],
    "rna_normal_col_indicators": ["normal", "control", "cortex", "healthy"],

    "rna_metadata_rows": {"SAMPLE_ID", "STUDY_ID", "PATIENT_ID",
                           "CANCER_TYPE", "ONCOTREE_CODE"},
}

# ─────────────────────────────────────────────────────────────────────────────
# ACVR1 SUBGROUP — for stratified scoring
# ─────────────────────────────────────────────────────────────────────────────

ACVR1_SUBGROUP = {
    # Genes that define the ACVR1-mutant subgroup
    "defining_genes": {"ACVR1", "BMPR1A", "BMPR2", "SMAD1", "SMAD5"},
    # Expected prevalence in DIPG cohort (~25%)
    "expected_prevalence": 0.25,
    # Score bonus for drugs targeting ACVR1 pathway in ACVR1-mutant subgroup
    "subgroup_bonus": 0.15,
}