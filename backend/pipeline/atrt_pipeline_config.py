"""
atrt_pipeline_config.py
========================
ATRT-specific configuration to be merged into pipeline_config.py.

INSTRUCTIONS FOR INTEGRATION
------------------------------
Add the following blocks to your existing pipeline_config.py.
Each block is clearly labelled with where it goes.

DATA SOURCES TO DOWNLOAD BEFORE RUNNING
-----------------------------------------

1. GSE70678 — Torchia et al. 2015 Cell. ATRT RNA-seq cohort.
   URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70678
   File needed: GSE70678_series_matrix.txt.gz (or processed expression matrix)
   Place at: data/raw_omics/GSE70678_ATRT_expression.txt
   Contains: 49 ATRT samples + subgroup labels (TYR/SHH/MYC) + normals.
   Platform: GPL6244 (Affymetrix HuGene 1.0 ST)

2. GSE106982 — Johann et al. 2016 Cancer Cell. Methylation + RNA subgroups.
   URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106982
   File needed: GSE106982_series_matrix.txt.gz
   Place at: data/raw_omics/GSE106982_ATRT_methylation.txt
   Contains: 150 ATRT samples with subgroup calls. Use for prevalence priors.

3. DepMap — already downloaded, but re-filter for ATRT/rhabdoid:
   OncotreeSubtype values to include: "MRT", "ATRT", "PRCC" (proximal renal)
   Cell lines: BT16, BT37, G401, A204, KP-MRT-NS, MON, CHLA02, CHLA04, CHLA06, BT12
   No new download needed — add "MRT" and "ATRT" to GBM_SUBTYPE_TERMS in depmap_essentiality.py

4. CBTN ATRT cohort (optional, controlled access):
   URL: https://portal.kidsfirstdrc.org
   Study: PBTA ATRT samples (same Cavatica portal as your DIPG data)
   ~30-40 ATRT samples with RNA-seq and methylation
   Place genomic files at: data/validation/cbtn_genomics/atrt/

5. CMap query for ATRT signature:
   Use GSE70678 to derive upregulated genes in ATRT vs normal brain.
   Run prepare_cmap_query.py with ATRT_INDICATORS instead of H3K27M_INDICATORS.
   Save to: data/cmap_query/atrt_cmap_scores.json
"""

# ─────────────────────────────────────────────────────────────────────────────
# ADD TO pipeline_config.py: PATHS section
# ─────────────────────────────────────────────────────────────────────────────

ATRT_PATHS = {
    "scrna_atrt":     "data/raw_omics/GSE70678_ATRT_expression.txt",
    "methylation":    "data/raw_omics/GSE106982_ATRT_methylation.txt",
    "genomics_atrt":  "data/validation/cbtn_genomics/atrt/",
    "cmap_atrt":      "data/cmap_query/atrt_cmap_scores.json",
}

# ─────────────────────────────────────────────────────────────────────────────
# ADD TO pipeline_config.py: DEPMAP section
# ─────────────────────────────────────────────────────────────────────────────

# Add these to GBM_SUBTYPE_TERMS in depmap_essentiality.py
ATRT_SUBTYPE_TERMS = [
    "atrt",
    "atypical teratoid",
    "rhabdoid",
    "malignant rhabdoid",
    "mrt",
]

# Known ATRT/rhabdoid DepMap cell line IDs
# Verify these against your Model.csv — IDs change between DepMap releases
ATRT_CELL_LINE_NAMES = [
    "BT16", "BT37", "BT12",
    "G401", "A204", "MON",
    "KP-MRT-NS", "KP-MRT-RY",
    "CHLA02", "CHLA04", "CHLA06",
]

# ─────────────────────────────────────────────────────────────────────────────
# ADD TO pipeline_config.py: COMPOSITE_WEIGHTS section
# ─────────────────────────────────────────────────────────────────────────────

# ATRT uses the same weight architecture as DIPG but with different tissue source
ATRT_COMPOSITE_WEIGHTS = {
    # tissue 0.40 — GSE70678 ATRT bulk RNA expression (no scRNA atlas exists yet)
    # depmap 0.35 — Broad CRISPR, rhabdoid/ATRT lines (BT16, BT37, G401, A204)
    # escape 0.20 — SMARCB1-loss resistance bypass
    # ppi    0.05 — STRING-DB (same as DIPG)
    "tissue": 0.40,
    "depmap": 0.35,
    "escape": 0.20,
    "ppi":    0.05,
}

# ─────────────────────────────────────────────────────────────────────────────
# ADD TO pipeline_config.py: BBB section
# ─────────────────────────────────────────────────────────────────────────────

# ATRT BBB penalties are LESS severe than DIPG.
# ATRT can be supratentorial (~35%) where BBB is standard, or
# infratentorial/posterior fossa (~50%) where it's tighter.
# Without location data, use moderate penalty.
ATRT_BBB_PENALTIES = {
    "infratentorial": {"LOW": 0.55, "UNKNOWN": 0.88},  # Similar to DIPG
    "supratentorial": {"LOW": 0.75, "UNKNOWN": 0.92},  # Standard BBB
    "unknown_location": {"LOW": 0.65, "UNKNOWN": 0.90},  # Conservative average
}

# ─────────────────────────────────────────────────────────────────────────────
# ADD TO pipeline_config.py: EZH2 section
# ─────────────────────────────────────────────────────────────────────────────

# CRITICAL: Override DIPG's EZH2_INHIBITOR penalty for ATRT
# In ATRT, EZH2 is a BOOST not a penalty.
# This must be handled in the disease-dispatch logic.
EZH2_ATRT_BOOST = {
    "known_inhibitors": {
        "tazemetostat", "gsk126", "epz-6438", "epz6438",
        "unc1999", "unc-1999", "cpi-1205", "cpi1205",
        "ds-3201", "ds3201", "valemetostat",
    },
    "composite_boost": 1.40,
    "tissue_boost": 1.35,
    "rationale": (
        "SMARCB1 loss → PRC2/EZH2 hyperactivity → EZH2 is ESSENTIAL in ATRT. "
        "Synthetic lethality confirmed Knutson 2013 PNAS. "
        "FDA Breakthrough Therapy designation for tazemetostat in SMARCB1-deficient tumors."
    ),
}

# ─────────────────────────────────────────────────────────────────────────────
# ADD TO pipeline_config.py: GENOMICS section
# ─────────────────────────────────────────────────────────────────────────────

ATRT_GENOMICS = {
    # SMARCB1 loss calling from mutation/CNA data
    "smarcb1_del_threshold":    -1,     # CNA score ≤ -1 = homozygous deletion
    "smarca4_del_threshold":    -1,
    "rna_upregulation_zscore":  1.0,    # z-score > 1.0 = upregulated in ATRT vs normal
    "rna_min_genes_required":   50,

    # Column indicators for ATRT vs normal in RNA data (GSE70678)
    "rna_atrt_col_indicators":  ["ATRT", "MRT", "rhabdoid", "tumor", "T"],
    "rna_normal_col_indicators":["normal", "control", "brain", "N"],

    # Subgroup column indicators in GSE70678
    "subgroup_indicators": {
        "TYR": ["TYR", "tyrosinase", "Type1", "Group1"],
        "SHH": ["SHH", "sonic", "Type2", "Group2"],
        "MYC": ["MYC", "Type3", "Group3"],
    },

    # SMARCB1 / SMARCA4 gene aliases for mutation calling
    "smarcb1_aliases": {"SMARCB1", "INI1", "SNF5", "BAF47"},
    "smarca4_aliases": {"SMARCA4", "BRG1", "BAF190"},

    # Metadata rows to skip in RNA matrix
    "rna_metadata_rows": {"SAMPLE_ID", "STUDY_ID", "SUBGROUP", "SMARCB1_STATUS"},
}

# ─────────────────────────────────────────────────────────────────────────────
# ADD TO pipeline_config.py: TISSUE section
# ─────────────────────────────────────────────────────────────────────────────

ATRT_TISSUE = {
    # No scRNA atlas for ATRT (unlike DIPG Filbin 2018).
    # Use bulk RNA (GSE70678) with curated cell-line scores as fallback.
    # Quantile threshold for high expression in ATRT bulk RNA cohort.
    "bulk_high_quantile": 0.75,   # p75 of GSE70678 cohort = high expression
    "bulk_mid_quantile":  0.50,
    "curated_weight":     0.60,   # Higher weight on curated (less scRNA data)
    "bulk_weight":        0.40,

    # Subgroup-specific expression multipliers
    # If subgroup is known, upweight subgroup-relevant targets
    "subgroup_multipliers": {
        "TYR": {"TYR": 1.5, "DCT": 1.5, "MITF": 1.5},
        "SHH": {"GLI2": 1.5, "SMO": 1.4, "PTCH1": 1.3},
        "MYC": {"MYC": 1.6, "MYCN": 1.5, "AURKA": 1.4},
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# CURATED ATRT TISSUE EXPRESSION SCORES
# Source: Torchia 2015 (bulk RNA), DepMap ATRT cell lines, published studies
# ─────────────────────────────────────────────────────────────────────────────

ATRT_GSC_CURATED_SCORES = {
    # Core synthetic lethality targets — SMARCB1-loss dependency
    "EZH2":    0.92,   # BOOSTED: hyperactive in SMARCB1-null. Knutson 2013.
    "EED":     0.88,
    "SUZ12":   0.85,
    # BET bromodomain
    "BRD4":    0.88,
    "BRD2":    0.82,
    "BRD3":    0.78,
    # HDAC
    "HDAC1":   0.85,
    "HDAC2":   0.82,
    "HDAC3":   0.80,
    "HDAC6":   0.65,
    # MYC axis
    "MYC":     0.85,
    "MYCN":    0.80,
    "MAX":     0.70,
    # Aurora kinase — MYCN stabiliser
    "AURKA":   0.82,
    "AURKB":   0.75,
    # CDK4/6
    "CDK4":    0.80,
    "CDK6":    0.75,
    "CCND1":   0.70,
    # mTOR/PI3K
    "MTOR":    0.70,
    "PIK3CA":  0.65,
    "AKT1":    0.62,
    # SHH subgroup
    "GLI2":    0.72,
    "GLI1":    0.68,
    "SMO":     0.65,
    # Stemness
    "SOX2":    0.80,
    "LIN28A":  0.78,
    "SALL4":   0.75,
    "SOX9":    0.72,
    # Proteasome
    "PSMB5":   0.68,
    "PSMB2":   0.65,
    "PSMB1":   0.62,
    # Apoptosis
    "BCL2":    0.62,
    "BCL2L1":  0.65,
    "MCL1":    0.68,
    # Tumour suppressors (lost/low)
    "SMARCB1": 0.05,   # LOST — defines the disease
    "SMARCA4": 0.05,   # LOST in ~5% of ATRT
    "CDKN2A":  0.20,
    "PTEN":    0.30,
    "RB1":     0.40,
    # DNA damage
    "PARP1":   0.65,
    "ATM":     0.55,
    # Immune
    "CD274":   0.52,
    # ONC201 target
    "DRD2":    0.50,
    "CLPB":    0.48,
    # SMARCB1 re-expression targets (research context)
    "SMARCC1": 0.30,
    "SMARCC2": 0.32,
}