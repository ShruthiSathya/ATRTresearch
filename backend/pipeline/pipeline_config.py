
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
    # Weights reflect signal reliability hierarchy for CNS drug prioritisation.
    # No single paper establishes this exact split; reasoning documented here:
    # tissue 0.40 — disease-specific GSC context (Filbin 2018, GSE102130)
    # depmap 0.30 — externally validated functional dependency
    #               (Behan et al. 2019 Nature; Tsherniak et al. 2017 Cell)
    # escape 0.20 — computed resistance estimate, informative but indirect
    # ppi    0.10 — network proximity is weakest causal signal
    # Sensitivity analysis: top-2 stable under all ±10% perturbations.
    # #3 position contested between Abemaciclib/Marizomib — both reported.
    "tissue":  0.40,
    "depmap":  0.30,
    "escape":  0.20,
    "ppi":     0.10,
}

assert abs(sum(COMPOSITE_WEIGHTS.values()) - 1.0) < 1e-9, \
    "COMPOSITE_WEIGHTS must sum to 1.0"

# Default values when a score key is absent from a candidate dict
# These represent "no data" — conservative neutral priors
SCORE_DEFAULTS = {
    "tissue_expression_score": 0.10,
    "depmap_score":            0.10,
    "ppi_score":               0.10,
    "escape_bypass_score":     0.40,   # Neutral — unknown escape risk
}

# ─────────────────────────────────────────────────────────────────────────────
# CONFIDENCE SCORE WEIGHTS  (must sum to 1.0)
# ─────────────────────────────────────────────────────────────────────────────

CONFIDENCE_WEIGHTS = {
    # depmap 0.45 — strongest external signal; Broad CRISPR screens are the
    #               most reproducible functional data in oncology
    #               (Behan et al. 2019; Meyers et al. 2017 Nat Genet)
    # bbb    0.35 — near-binary requirement for CNS drugs; non-crossing drugs
    #               are essentially non-starters (Pardridge 2003 PMC539316;
    #               Fischer et al. 1998 — 400 Da MW rule)
    # diversity 0.20 — mechanistic diversity reduces combination toxicity risk;
    #               lower weight as Jaccard overlap is a proxy, not direct measure
    "depmap":    0.45,
    "bbb":       0.35,
    "diversity": 0.20,
}

assert abs(sum(CONFIDENCE_WEIGHTS.values()) - 1.0) < 1e-9, \
    "CONFIDENCE_WEIGHTS must sum to 1.0"

# Default depmap component when data not loaded (conservative prior)
DEPMAP_MISSING_PRIOR = 0.30

# ─────────────────────────────────────────────────────────────────────────────
# TISSUE EXPRESSION SCORING
# ─────────────────────────────────────────────────────────────────────────────

TISSUE = {
    # Blending: curated DIPG knowledge vs real scRNA-seq signal
    # Equal weighting valid because reference is H3K27M DIPG (GSE102130)
    "curated_weight": 0.50,
    "sc_weight":      0.50,

    # Quantile used to select stem-like (GSC) cells from the scRNA matrix
    "gsc_stem_quantile": 0.85,

    # Stem cell marker genes used to identify GSC cells
    "stem_markers": ["SOX2", "NES", "PROM1", "CD44"],

    # Percentile bin boundaries → score mapping
    # Keys are percentile thresholds (p90, p75, p50, p25); value = assigned score
    "percentile_bins": {
        "p90": 1.00,
        "p75": 0.82,
        "p50": 0.62,
        "p25": 0.42,
        "low": 0.18,
    },

    # Fallback scores when no percentile data computed
    "fallback_high_tpm":   1.0,   # expr >= 2.0 TPM
    "fallback_mid_tpm":    0.6,   # expr >= 0.5 TPM
    "fallback_low_tpm":    0.2,
    "fallback_high_cutoff": 2.0,
    "fallback_mid_cutoff":  0.5,

    # DIPG-relevance discount for genes not in the curated list
    # Prevents housekeeping genes dominating due to ubiquitous expression
    "off_target_relevance": 0.35,

    # Neutral score for unknown targets
    "unknown_target_score": 0.40,

    # Curated score blend: best target vs mean of top-N
    "curated_best_weight":    0.60,
    "curated_top_n":          3,
    "curated_top_n_weight":   0.40,

    # Score for drug with no targets listed
    "no_target_score": 0.40,
}

assert abs(TISSUE["curated_weight"] + TISSUE["sc_weight"] - 1.0) < 1e-9, \
    "TISSUE curated_weight + sc_weight must sum to 1.0"

# ─────────────────────────────────────────────────────────────────────────────
# ESCAPE BYPASS SCORING
# ─────────────────────────────────────────────────────────────────────────────

ESCAPE = {
    # Scores by number of active bypass nodes detected
    "bypass_scores": {
        0: 1.00,   # No active bypass — drug should work
        1: 0.85,
        2: 0.72,
        3: 0.58,
        4: 0.40,   # 4+ bypass routes — high resistance risk (used as floor)
    },

    # Score when drug has no targets in the bypass map
    "no_target_score":    0.75,
    # Score when drug has no targets at all
    "empty_target_score": 0.70,

    # Blend when live STRING-DB escape hits are found
    "string_weight":  0.60,
    "curated_weight": 0.40,

    # Score used for string_score when live hits found
    # (reflects that live hit = confirmed active bypass = bad)
    "string_hit_score": 0.40,

    # Genes that are constitutively active resistance nodes in DIPG
    # regardless of RNA expression (hardwired biology)
    # RAD51/BRCA1 added: H3K27M DIPG is HR-proficient — PARP inhibitors
    # have no synthetic lethality rationale in this tumour type.
    "constitutive_resistance_nodes": {
        "PIK3CA", "MTOR", "MYC", "MYCN", "BCL2L1", "MCL1", "CDK4", "BRD4",
        "RAD51", "BRCA1",
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# TOXICITY CONSTRAINT
# ─────────────────────────────────────────────────────────────────────────────

TOXICITY = {
    # Maximum combined hematologic AE rate before multiplier floors
    "max_combined_rate": 0.40,

    # Default rate for drugs without published data
    "default_rate": 0.05,

    # Thresholds calibrated for pediatric DIPG — no curative standard of care.
    # PBTC-047 (panobinostat): 8/29 DLTs = 27.6% (Monje et al. 2023, PMID 37526549)
    # PBTC-049 (birabresib):   grade 3/4 haem AEs in ~15-20% (Hennika et al.)
    # These trials proceeded at rates that would be HIGH_RISK by adult standards,
    # establishing that the pediatric CNS oncology field accepts higher AE rates
    # when there is no alternative. Thresholds set relative to published PBTC data.
    "acceptable_threshold": 0.20,   # ≤20%: single-agent equivalent risk
    "caution_threshold":    0.33,   # 20-33%: manageable with dose reduction (PBTC precedent)
    # >33% → HIGH_RISK: requires prospective toxicity management plan
}

# ─────────────────────────────────────────────────────────────────────────────
# BBB FILTER
# ─────────────────────────────────────────────────────────────────────────────

BBB = {
    # Hard molecular weight cutoff — above this, drug cannot cross BBB
    "hard_exclude_mw": 800.0,

    # MW heuristic thresholds — Pardridge 2003 (PMC539316); Fischer et al. 1998
    # < 400 Da: free diffusion through BBB (Lipinski rule of five)
    # 400-600 Da: reduced penetrance
    # > 600 Da: minimal CNS penetration without active transport
    "mw_moderate_cutoff": 400.0,
    "mw_low_cutoff":      600.0,

    # Numeric scores per penetrance category
    "penetrance_scores": {
        "HIGH":    1.0,
        "MODERATE": 0.6,
        "LOW":     0.2,
        "UNKNOWN": 0.5,
    },

    # Score for known GBM clinical trial failures
    "failure_score": 0.05,

    # Multiplier applied to composite score for known failures
    "failure_composite_multiplier": 0.10,

    # Heuristic scores
    "heuristic_moderate_score": 0.6,
    "heuristic_low_near_score": 0.3,
    "heuristic_low_score":      0.1,
    "unknown_score":            0.5,
}

# ─────────────────────────────────────────────────────────────────────────────
# PPI NETWORK SCORES
# ─────────────────────────────────────────────────────────────────────────────

PPI = {
    "first_degree_score":  0.85,   # Target is 1st-degree STRING-DB neighbor of disease gene
    "second_degree_score": 0.60,   # Target is 2nd-degree neighbor (curated only)
    "no_proximity_score":  0.20,   # No network proximity found
}

# ─────────────────────────────────────────────────────────────────────────────
# HYPOTHESIS GENERATOR THRESHOLDS
# ─────────────────────────────────────────────────────────────────────────────

HYPOTHESIS = {
    # p-value threshold for calling genomic co-occurrence significant
    "p_value_significance": 0.05,

    # Escape bypass score minimum for a drug to be considered "bypass-resistant"
    "bypass_high_threshold": 0.50,

    # DepMap default score — scores near this value indicate data not loaded
    "depmap_default_score":    0.50,
    "depmap_default_tolerance": 0.01,

    # Fallback depmap score when data confirmed missing
    "depmap_missing_score": 0.10,

    # Top N diverse drugs to select for combination hypothesis
    "combo_size": 3,

    # Default component scores when a key is absent from candidate dict
    "missing_depmap_score": 0.10,
}

# ─────────────────────────────────────────────────────────────────────────────
# GENOMIC VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

GENOMICS = {
    # log2FC > 1.0 = ≥2-fold higher in tumour vs normal.
    # Standard DEG threshold — Love et al. 2014 (DESeq2, PMC4302049);
    # Used across GSE50021 (Grasso/Monje, 35 DIPG + 10 normal, microarray)
    # and GSE115397 (Nagaraja et al., 5 H3K27M pons + 3 normal, RNA-seq).
    "rna_upregulation_zscore": 1.0,

    # CNA score threshold for calling CDKN2A deleted
    "cna_deletion_threshold": -1,

    # Minimum gene rows to treat RNA file as a real expression matrix
    "rna_min_genes_required": 50,

    # H3K27M mutation identifiers
    "h3k27m_values": {"K28M", "K27M"},
    "h3_gene_aliases": {"H3-3A", "H3F3A", "HIST1H3B", "H33A", "H3.3A"},

    # Substrings identifying H3K27M tumour columns in a reference cohort
    # RNA file (e.g. GSE115397) when sample IDs don't match the genomic cohort
    "rna_h3k27m_col_indicators": ["Pons", "K27M", "DIPG", "H3K27M", "pontine"],

    # Substrings identifying normal/control columns — excluded from upregulation
    "rna_normal_col_indicators": ["normal", "control", "cortex", "healthy"],

    # Metadata row names to skip when counting gene rows
    "rna_metadata_rows": {"SAMPLE_ID", "STUDY_ID", "PATIENT_ID",
                           "CANCER_TYPE", "ONCOTREE_CODE"},
}