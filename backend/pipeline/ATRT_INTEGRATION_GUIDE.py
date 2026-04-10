"""
EXISTING FILE MODIFICATIONS FOR ATRT SUPPORT
=============================================
These are the targeted changes needed in your existing pipeline files.
Each section is labelled with the file and the exact location of the change.

──────────────────────────────────────────────────────────────────────────────
FILE 1: backend/pipeline/depmap_essentiality.py
──────────────────────────────────────────────────────────────────────────────

CHANGE: Add ATRT/rhabdoid cell line terms to GBM_SUBTYPE_TERMS list.
LOCATION: Top of file, in the GBM_SUBTYPE_TERMS list (around line 17).

BEFORE:
    GBM_SUBTYPE_TERMS = [
        "glioblastoma",
        "diffuse intrinsic pontine",
        ...
    ]

AFTER:
    GBM_SUBTYPE_TERMS = [
        "glioblastoma",
        "diffuse intrinsic pontine",
        "diffuse midline glioma",
        "high grade glioma",
        "pediatric gbm",
        "anaplastic glioma",
        "grade iv glioma",
        "gliosarcoma",
        "giant cell glioblastoma",
        # ATRT/rhabdoid additions
        "atrt",
        "atypical teratoid",
        "rhabdoid",
        "malignant rhabdoid",
        "mrt",
    ]

WHY: DepMap has ATRT/rhabdoid lines (BT16, BT37, G401, A204, CHLA02/04/06)
     that score with OncotreeSubtype "MRT" or "ATRT". These are the most
     relevant cell lines for SMARCB1-null biology.

──────────────────────────────────────────────────────────────────────────────
FILE 2: backend/pipeline/tissue_expression.py
──────────────────────────────────────────────────────────────────────────────

CHANGE 1: Import ATRT curated scores.
LOCATION: Top of file, after DIPG_GSC_CURATED_SCORES definition.

ADD after the DIPG_GSC_CURATED_SCORES dict:

    # ATRT curated scores — imported from atrt_pipeline_config
    try:
        from .atrt_pipeline_config import ATRT_GSC_CURATED_SCORES
    except ImportError:
        ATRT_GSC_CURATED_SCORES = {}

CHANGE 2: Make TissueExpressionScorer disease-aware.
LOCATION: In _curated_score() method.

BEFORE:
    def _curated_score(self, drug: Dict) -> float:
        targets = drug.get("targets", [])
        if not targets:
            return TISSUE["no_target_score"]
        scores = [DIPG_GSC_CURATED_SCORES.get(t.upper(), TISSUE["unknown_target_score"])
                  for t in targets]

AFTER:
    def _curated_score(self, drug: Dict) -> float:
        targets = drug.get("targets", [])
        if not targets:
            return TISSUE["no_target_score"]

        # Select curated score database based on disease
        if "atrt" in (self.disease_name or "").lower() or "rhabdoid" in (self.disease_name or "").lower():
            curated_db = ATRT_GSC_CURATED_SCORES if ATRT_GSC_CURATED_SCORES else DIPG_GSC_CURATED_SCORES
        else:
            curated_db = DIPG_GSC_CURATED_SCORES

        scores = [curated_db.get(t.upper(), TISSUE["unknown_target_score"])
                  for t in targets]

CHANGE 3: Add ATRT bulk RNA loader.
LOCATION: In _load_sc_data() method, before the scRNA loading block.

ADD this block at the start of _load_sc_data():

    async def _load_sc_data(self, quantile: float = None):
        q = quantile or TISSUE["gsc_stem_quantile"]
        if self.is_ready and quantile is None:
            return

        # ATRT: use bulk RNA (GSE70678) instead of scRNA (no ATRT atlas exists)
        if "atrt" in (self.disease_name or "").lower():
            await self._load_atrt_bulk_rna()
            return

        # existing DIPG scRNA loading continues below...

    async def _load_atrt_bulk_rna(self):
        '''Load GSE70678 bulk RNA-seq for ATRT tissue scoring.'''
        from .atrt_pipeline_config import ATRT_PATHS, ATRT_TISSUE, ATRT_GENOMICS
        bulk_path = Path(ATRT_PATHS["scrna_atrt"])

        if not bulk_path.exists():
            logger.warning(
                "GSE70678 not found at %s. Using curated ATRT scores only. "
                "Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70678",
                bulk_path,
            )
            return

        try:
            import pandas as pd
            logger.info("Loading GSE70678 ATRT bulk RNA-seq...")
            df = pd.read_csv(bulk_path, sep="\\t", index_col=0)
            df.index = df.index.str.upper().str.strip()

            atrt_indicators = ATRT_GENOMICS["rna_atrt_col_indicators"]
            normal_indicators = ATRT_GENOMICS["rna_normal_col_indicators"]

            atrt_cols = [c for c in df.columns
                         if any(ind.lower() in c.lower() for ind in atrt_indicators)
                         and not any(n.lower() in c.lower() for n in normal_indicators)]

            if not atrt_cols:
                logger.warning("No ATRT columns found in GSE70678. Check column names.")
                return

            # Compute mean expression across ATRT samples
            mean_expr = df[atrt_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
            self.gsc_target_scores = mean_expr.to_dict()

            # Compute percentiles for _expression_to_score()
            all_vals = sorted(mean_expr.dropna().values)
            n = len(all_vals)
            self._gsc_all_values_sorted = all_vals
            self._gsc_p25 = all_vals[int(n * 0.25)]
            self._gsc_p50 = all_vals[int(n * 0.50)]
            self._gsc_p75 = all_vals[int(n * 0.75)]
            self._gsc_p90 = all_vals[int(n * 0.90)]

            self.is_ready = True
            logger.info(
                "✅ GSE70678 loaded: %d ATRT samples, %d genes. "
                "p25=%.2f, p50=%.2f, p75=%.2f, p90=%.2f",
                len(atrt_cols), len(all_vals),
                self._gsc_p25, self._gsc_p50, self._gsc_p75, self._gsc_p90,
            )
        except Exception as e:
            logger.error("GSE70678 loading failed: %s", e)

──────────────────────────────────────────────────────────────────────────────
FILE 3: backend/pipeline/bbb_filter.py
──────────────────────────────────────────────────────────────────────────────

NO CHANGES NEEDED to the BBBFilter class itself.
The ATRT BBB penalty is applied in run_atrt_pipeline.py using
ATRT_BBB_PENALTIES from atrt_pipeline_config.py, which is location-aware.

The key difference: ATRT does NOT get the severe 0.50x brainstem penalty
that DIPG gets. Instead it uses ATRT_BBB_PENALTIES["unknown_location"]
which is {"LOW": 0.65, "UNKNOWN": 0.90} — more moderate.

──────────────────────────────────────────────────────────────────────────────
FILE 4: backend/pipeline/pipeline_config.py
──────────────────────────────────────────────────────────────────────────────

CHANGE: The EZH2_INHIBITOR block defines a PENALTY for DIPG.
This is correct for DIPG and must NOT be changed.
The ATRT EZH2 BOOST is handled separately in atrt_specialization.py.

However, add a disease-dispatch note comment:

ADD at the end of pipeline_config.py:

    # ─────────────────────────────────────────────────────────────────────────
    # DISEASE DISPATCH NOTE
    # EZH2 treatment differs by disease:
    #   DIPG:  EZH2 inhibitors PENALISED (×0.50) — H3K27M suppresses PRC2
    #   ATRT:  EZH2 inhibitors BOOSTED  (×1.40) — SMARCB1 loss needs EZH2
    # This is handled in dipg_specialization.py vs atrt_specialization.py.
    # Never apply the DIPG EZH2 penalty when disease is ATRT.
    # ─────────────────────────────────────────────────────────────────────────

──────────────────────────────────────────────────────────────────────────────
FILE 5: backend/pipeline/discovery_pipeline.py  (ProductionPipeline.run())
──────────────────────────────────────────────────────────────────────────────

CHANGE: The run() method currently hard-applies the DIPG BBB penalty.
For ATRT, skip the dipg_bbb_penalties block and use ATRT_BBB_PENALTIES instead.

LOCATION: In the ProductionPipeline.run() method, the DIPG BBB penalty block.

BEFORE:
    # DIPG brainstem BBB penalty
    dipg_bbb_penalties = {"LOW": 0.50, "UNKNOWN": 0.85}
    n_penalised = 0
    for c in candidates:
        bbb = c.get("bbb_penetrance", "UNKNOWN")
        if bbb in dipg_bbb_penalties:
            c["score"] = round(c["score"] * dipg_bbb_penalties[bbb], 4)

AFTER:
    # Disease-aware BBB penalty
    is_atrt = any(k in disease_name.lower() for k in ("atrt", "rhabdoid"))
    if is_atrt:
        from .atrt_pipeline_config import ATRT_BBB_PENALTIES
        bbb_penalties = ATRT_BBB_PENALTIES["unknown_location"]
        penalty_label = "ATRT (location-unknown)"
    else:
        # DIPG: severe brainstem penalty
        bbb_penalties = {"LOW": 0.50, "UNKNOWN": 0.85}
        penalty_label = "DIPG brainstem"

    n_penalised = 0
    for c in candidates:
        bbb = c.get("bbb_penetrance", "UNKNOWN")
        if bbb in bbb_penalties:
            c["score_before_bbb"] = c["score"]
            c["score"] = round(c["score"] * bbb_penalties[bbb], 4)
            n_penalised += 1

    if n_penalised:
        logger.info(
            "%s BBB penalty: %d candidates penalised", penalty_label, n_penalised
        )

──────────────────────────────────────────────────────────────────────────────
FILE 6: backend/pipeline/hypothesis_generator.py
──────────────────────────────────────────────────────────────────────────────

CHANGE: The hypothesis generator currently surfaces DIPG-specific notes.
Make the EZH2 note disease-aware.

LOCATION: In the generate() method, the ezh2_warnings block.

BEFORE:
    ezh2_warnings = [
        c.get("name", "?") for c in top_3
        if c.get("dipg_components", {}).get("is_ezh2_inhibitor")
    ]

AFTER:
    # DIPG: EZH2 inhibitors in combo are a WARNING (non-rational)
    # ATRT: EZH2 inhibitors in combo are EXPECTED (synthetic lethality)
    is_atrt_disease = any(
        k in (genomic_stats or {}).get("disease", "").lower()
        for k in ("atrt", "rhabdoid")
    )

    if is_atrt_disease:
        ezh2_warnings = []  # Not a warning in ATRT
        ezh2_boosts = [
            c.get("name", "?") for c in top_3
            if c.get("atrt_components", {}).get("ezh2_boosted")
        ]
    else:
        ezh2_warnings = [
            c.get("name", "?") for c in top_3
            if c.get("dipg_components", {}).get("is_ezh2_inhibitor")
        ]
        ezh2_boosts = []

──────────────────────────────────────────────────────────────────────────────
NEW FILES TO ADD TO backend/pipeline/
──────────────────────────────────────────────────────────────────────────────

Copy these files from /home/claude/ to backend/pipeline/:
  - atrt_specialization.py
  - atrt_pipeline_config.py
  - run_atrt_pipeline.py
  - published_ic50_atrt_validation.py

──────────────────────────────────────────────────────────────────────────────
DATA DOWNLOAD COMMANDS
──────────────────────────────────────────────────────────────────────────────

# GSE70678 — primary ATRT RNA-seq
mkdir -p data/raw_omics
# Download via GEOquery (R) or GEO FTP:
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70678/matrix/GSE70678_series_matrix.txt.gz
# Rename/process to: data/raw_omics/GSE70678_ATRT_expression.txt

# GSE106982 — subgroup methylation (optional but useful for subgroup priors)
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106982/matrix/

# Python GEO download (no R required):
# pip install GEOparse
# python -c "
# import GEOparse
# gse = GEOparse.get_GEO('GSE70678', destdir='data/raw_omics/')
# # Extract expression matrix and save as TSV
# "

# Run ATRT pipeline:
# python -m backend.pipeline.run_atrt_pipeline --disease atrt --top_n 20

──────────────────────────────────────────────────────────────────────────────
VALIDATION CHECKLIST
──────────────────────────────────────────────────────────────────────────────

After implementing, verify these outputs:
1. Tazemetostat should appear in top-5 (EZH2 inhibitor, SMARCB1 synthetic lethality)
2. Alisertib should appear in top-10 (AURKA inhibitor)
3. Panobinostat should appear in top-5 (pan-HDAC, ATRT-validated)
4. Birabresib should appear in top-5 (BET bromodomain)
5. Tazemetostat score should be HIGHER in ATRT than in DIPG (opposite of DIPG penalty)
6. EZH2 inhibitors should show "ezh2_boosted: True" in atrt_components
7. BT16/G401 IC50 validation should annotate top candidates
"""