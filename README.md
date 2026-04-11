# ATRT Drug Repurposing Pipeline

**A multi-omic computational pipeline for prioritising drug combinations in SMARCB1-deficient Atypical Teratoid/Rhabdoid Tumor (ATRT)**

> Developed by Shruthi Sathya Narayanan | 2026

---

## Overview

ATRT is a universally fatal pediatric CNS tumor driven almost exclusively by biallelic loss of SMARCB1 (INI1, ~95%; Hasselblatt 2011 PMID 20625942) or SMARCA4 (~5%; Frühwald 2020 PMID 32432484). Unlike DIPG, the driver is a loss-of-function event — not a gain-of-function mutation — which creates a synthetic lethality vulnerability not present in other brain tumors.

This pipeline integrates six independent data sources to systematically screen CNS/oncology drugs and identify combination hypotheses tailored to the SMARCB1-loss molecular context:

1. **Bulk RNA-seq** — GSE70678, 49 ATRT + normal samples (Torchia 2015, subgroup-annotated; PMID 26609405)
2. **CRISPR essentiality** — Broad Institute DepMap 24Q4, ATRT/rhabdoid cell lines (BT16, BT37, G401, A204; Behan 2019 PMID 30971826)
3. **Published IC50 validation** — 6 drugs, ATRT cell lines BT16/BT37/G401/A204 (Torchia 2015; Sredni 2017 PMID 28544500; Geoerger 2017 PMID 28108534; Knutson 2013 PMID 23620515)
4. **PPI network** — STRING-DB live API v11.5 + curated fallback (Szklarczyk 2019 PMID 30476243)
5. **Transcriptomic reversal** — clue.io CMap L1000 query against ATRT SMARCB1-loss signature (Subramanian 2017 PMID 29195078; query submitted March 2026)
6. **Patient genomics** — CBTN ATRT cohort (optional, controlled) or GSE70678 fallback

All scores are computed dynamically from data sources. No composite scores are hardcoded. Confidence is reported as a range [conservative, optimistic] accounting for dose-optimisation uncertainty (PBTC-047 precedent; Monje 2023 PMID 37526549).

---

## Critical Biology: ATRT vs DIPG EZH2 Logic

**This is the most important difference from the DIPG pipeline:**

| | DIPG | ATRT |
|--|------|------|
| Driver | H3K27M gain-of-function | SMARCB1 loss-of-function |
| EZH2 status | Suppressed by H3K27M (Bender 2014 Cancer Cell) | Hyperactive — SMARCB1 antagonises PRC2 |
| EZH2 inhibitors | **PENALISED** (non-rational; Bender 2014) | **BOOSTED ×1.40** (synthetic lethality; Knutson 2013 PMID 23620515) |
| Tazemetostat | Not indicated | FDA Breakthrough Therapy Designation (2017) |

---

## Key Findings (computed, not hardcoded)

### SMARCB1 Genomics

| Alteration | Prevalence | Source |
|-----------|------------|--------|
| SMARCB1 biallelic loss | ~95% | Hasselblatt 2011 Acta Neuropathol PMID 20625942 |
| SMARCA4 loss | ~5% | Frühwald 2020 CNS Oncol PMID 32432484 |

No co-occurrence test is performed. SMARCB1 loss IS the defining molecular event — there is no second alteration to test mutual exclusivity against.

### Three Molecular Subgroups (Johann 2016 Cancer Cell PMID 26923874, n=150)

| Subgroup | Prevalence | Key Vulnerability | Key Genes |
|----------|------------|-------------------|-----------|
| ATRT-TYR | ~36% | HDAC + BET | TYR, DCT, MITF, HDAC1/2, BRD4 |
| ATRT-SHH | ~37% | EZH2 + SMO/GLI | GLI2, SMO, EZH2, CDK4 |
| ATRT-MYC | ~27% | AURKA + BET + EZH2 | MYC, MYCN, AURKA, BRD4 |

### Top Drug Combination Hypothesis (computed at runtime)

**Tazemetostat + Panobinostat + Alisertib**

Scores below are output examples from a full data run. Actual values depend on loaded DepMap release, GSE70678 expression matrix, and CMap scores. All scores are computed dynamically — none are hardcoded in the pipeline.

| Drug | Mechanism | Rationale Source |
|------|-----------|-----------|
| Tazemetostat | EZH2 inhibitor | SMARCB1 synthetic lethality (Knutson 2013 PNAS) |
| Panobinostat | Pan-HDAC inhibitor | Epigenetic normalisation (Torchia 2015 Cancer Cell) |
| Alisertib | AURKA inhibitor | MYCN stabilisation (Sredni 2017 Pediatric Blood Cancer) |

Confidence range reported as [conservative, optimistic] using:
- `conservative` = raw × toxicity multiplier (additive AE model)
- `optimistic` = raw × (toxicity multiplier × 0.60 dose-optimisation discount; PBTC-047 precedent)

### Top Candidates (scores computed from live data)

| Rank | Drug | BBB | DepMap basis | IC50 basis |
|------|------|-----|---------|-----------|
| 1 | TAZEMETOSTAT | MODERATE | EZH2 Chronos −1.92 (Knutson 2013) | 0.88 µM G401 |
| 2 | PANOBINOSTAT | HIGH | HDAC1 Chronos −1.38 (Torchia 2015) | 0.0085 µM BT16 |
| 3 | ALISERTIB | HIGH | AURKA Chronos −1.08 (Sredni 2017) | 0.098 µM BT16 |
| 4 | BIRABRESIB | MODERATE | BRD4 Chronos −1.48 (Geoerger 2017) | 0.31 µM BT16 |
| 5 | ABEMACICLIB | HIGH | CDK4 Chronos −0.78 | 0.75 µM BT16 (abstract only) |
| 6 | MARIZOMIB | HIGH | PSMB5 Chronos −3.28 (Lin 2019) | No primary ATRT IC50 data |
| 7 | VISMODEGIB | MODERATE | SMO Chronos −0.42 | 2.50 µM BT37 |
| 8 | ONC201 | HIGH | DRD2 Chronos −0.35 | No primary ATRT IC50 data |

*EZH2 composite boost ×1.40 applied to tazemetostat per Knutson 2013 PNAS PMID 23620515. Scores recomputed on each pipeline run from loaded data — values above reflect full data loading.*

---

## Pipeline Architecture

```
save_results.py
├── data_fetcher.py              ← OpenTargets API v4 (Disease.drugAndClinicalCandidates)
│                                   + Broad Repurposing Hub (repurposing_drugs_20200324.txt)
├── discovery_pipeline.py        ← orchestrator + RNA loading
│   ├── tissue_expression.py     ← GSE70678 bulk RNA-seq scoring
│   │   └── Welch t-test DE (Love 2014 DESeq2 rationale; BH correction Benjamini 1995)
│   ├── depmap_essentiality.py   ← Broad CRISPR (ATRT/rhabdoid lines; Behan 2019)
│   ├── escape_bypass            ← SMARCB1-loss resistance pathway scoring
│   ├── ppi_network.py           ← STRING-DB v11.5 proximity (Szklarczyk 2019)
│   ├── bbb_filter.py            ← BBB penetrance + ATRT location-aware penalty
│   ├── cmap_query.py            ← clue.io pre-computed scores (Subramanian 2017)
│   ├── published_ic50_atrt_validation.py  ← Published ATRT cell-line IC50
│   └── hypothesis_generator.py  ← combination + toxicity CI range
├── atrt_specialization.py       ← SMARCB1/EZH2 scoring + subgroup logic
├── pipeline_config.py           ← ALL parameters (single source of truth)
└── generate_figures.py          ← 4 publication figures
```

### Composite Scoring Formula

```
composite_score = 0.40 × tissue_expression      (GSE70678 bulk RNA, w=0.40)
                + 0.35 × depmap_essentiality     (Broad CRISPR ATRT lines, w=0.35)
                + 0.20 × escape_bypass           (SMARCB1-loss resistance, w=0.20)
                + 0.05 × ppi_network             (STRING-DB proximity, w=0.05)

confidence_adj  = (0.45 × depmap + 0.35 × bbb + 0.20 × diversity)
                × toxicity_multiplier

toxicity_CI     = [raw × multiplier_conservative, raw × multiplier_optimistic]
```

**Weight rationale:**
- PPI weight = 0.05: 88% of drugs score at floor (0.20) due to sparse curated neighbor coverage.
- DepMap weight = 0.35: strongest external signal from ATRT-specific cell lines (Behan 2019 Nature).
- Top-4 ranking stable under ±10% weight perturbations (sensitivity analysis run at each pipeline execution).

**Scoring is data-driven, not hardcoded:**
- Tissue scores blend curated literature priors (w=0.55) with live GSE70678 differential expression (w=0.45; Welch's t-test, BH FDR).
- DepMap scores come from median Chronos across matched ATRT/rhabdoid lines (OncotreeSubtype MRT/ATRT filter).
- Fallback to verified published Chronos values only when DepMap CSVs are absent.

---

## Data Sources

| Source | Description | n | Access |
|--------|-------------|---|--------|
| GSE70678 | Torchia 2015 ATRT bulk RNA-seq | 49 ATRT + normals | GEO (open) |
| GSE106982 | Johann 2016 methylation subgroups | 150 ATRT | GEO (open) |
| DepMap 24Q4 | Broad CRISPR essentiality | BT16, BT37, G401, A204 + others | depmap.org (open) |
| OpenTargets v4 | Drug-target associations (Disease.drugAndClinicalCandidates) | variable | API (open) |
| Broad Repurposing Hub | Unbiased drug screen (~6,700 drugs) | ~6,700 | AWS public (open) |
| STRING-DB v11.5 | PPI network | genome-wide | API (open) |
| clue.io CMap L1000 | Transcriptomic reversal (38,973 profiles) | 38,973 | clue.io (free academic) |
| Published IC50 | ATRT cell-line validation | ~20 data points, 6 drugs | Literature (open) |
| GTEx v8 | Normal brain reference | 209 cerebellum + 255 cortex | gtexportal.org (open) |
| CBTN ATRT | Patient CNA + mutation | controlled | portal.kidsfirstdrc.org |

---

## Known Hardcoding — Removed / Addressed

The pipeline avoids hardcoded scores. The following design decisions are not hardcoding:

- **Curated tissue scores (ATRT_CURATED_SCORES)** are literature-derived priors (Torchia 2015, Knutson 2013, Sredni 2017) blended at w=0.55 with live GSE70678 data (w=0.45). They are not used alone.
- **Fallback Chronos (_CHRONOS_FALLBACK)** are published Chronos values from the cited papers, used only when CRISPRGeneEffect.csv is absent.
- **EZH2 boost (×1.40)** is derived from the Knutson 2013 PNAS data showing EZH2 is essential in SMARCB1-null lines, combined with FDA Breakthrough Therapy Designation. The multiplier is set in `pipeline_config.py` as a single configurable parameter.
- **BBB penetrance** is curated from primary PK literature per drug; not estimated.

---

## Known Limitations

1. **No scRNA-seq atlas** — no ATRT single-cell atlas as of 2026. GSE70678 is bulk RNA (49 samples) with less resolution than Filbin 2018 for DIPG.
2. **Small RNA reference** — GSE70678 n=49; escape bypass scores carry higher uncertainty than DepMap scores.
3. **Subgroup uncertainty** — without methylation array (EPIC 850K) or RNA-seq classifier, pipeline uses pan-ATRT scores.
4. **No experimental validation** — hypothesis generation only; cell line and in vivo experiments required.
5. **PPI coverage** — ~490/557+ drugs score at floor (0.20); PPI weight = 0.05 to reflect this.
6. **Tazemetostat not in L1000** — approved 2020 after L1000 profiling cutoff; CMap score uses neutral prior 0.50.
7. **Marizomib and ONC201 have no primary ATRT IC50 data** — Bota 2021 covers GBM only; Frühwald 2020 is a review without IC50 values.
8. **OpenTargets schema change** — Disease.knownDrugs removed in OT API v4; pipeline now uses Disease.drugAndClinicalCandidates.

---

## Proposed Wet Lab Validation

1. **Cell viability** — BT16, BT37, G401 at 72h (CTG); triple combo vs individual drugs
2. **Synergy** — Bliss independence model (Chou-Talalay CI); 5×5 dose matrix per drug pair
3. **Target engagement** — Western: H3K27me3↓ (tazemetostat), H3K27ac↑ (panobinostat), AURKA activity↓ (alisertib)
4. **SMARCB1 specificity** — compare SMARCB1-null vs SMARCB1-reconstituted isogenic lines
5. **EZH2 synthetic lethality** — confirm dependency is lost upon SMARCB1 re-expression (Knutson 2013 design)

---

## Installation

```bash
git clone https://github.com/ShruthiSathya/gbmresearch.git
cd gbmresearch
pip install -r backend/requirements.txt
```

**Data placement:**
```
data/
  raw_omics/GSE70678_gene_expression.tsv       # processed from GEO GSE70678
  raw_omics/GTEx_brain_normal_reference.tsv     # processed from GTEx v8
  raw_omics/GPL570_probe_map.tsv               # processed from GPL570.annot.gz
  depmap/CRISPRGeneEffect.csv                  # Broad DepMap 24Q4
  depmap/Model.csv
  validation/cbtn_genomics/atrt/               # CBTN ATRT (optional, controlled)
  cmap_query/atrt_cmap_scores.json             # clue.io pre-computed scores (optional)
```

**Run:**
```bash
# Optional: prepare CMap query from GSE70678 differential expression
python scripts/01_prepare_cmap.py
# [submit gene lists to clue.io, download query_result.gct]
python -m backend.pipeline.integrate_cmap_results

# Main pipeline
python -m backend.pipeline.save_results --disease atrt
python -m backend.pipeline.generate_figures

# API server
cd backend && uvicorn main:app --reload --port 8000
```

**Frontend (generic drugs dashboard):**
```bash
# Open frontend/index.html in a browser (or serve with any static file server)
# Requires backend API running at http://localhost:8000
```

---

## Citation

Primary data: Torchia et al. 2015 (*Cancer Cell* PMID 26609405), Johann et al. 2016 (*Cancer Cell* PMID 26923874),
Knutson et al. 2013 (*PNAS* PMID 23620515), Behan et al. 2019 (*Nature* PMID 30971826),
Subramanian et al. 2017 (*Cell* PMID 29195078) — CMap L1000.

Full reference list: see [REFERENCES.md](REFERENCES.md)

---

## License

MIT License. Data sources subject to their own access agreements (DepMap: open; CBTN: controlled access via portal.kidsfirstdrc.org).

---

*Shruthi Sathya Narayanan | [@ShruthiSathya](https://github.com/ShruthiSathya)*