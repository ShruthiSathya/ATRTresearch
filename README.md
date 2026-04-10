# ATRT Drug Repurposing Pipeline

**A multi-omic computational pipeline for prioritising drug combinations in SMARCB1-deficient Atypical Teratoid/Rhabdoid Tumor (ATRT)**

> Developed by Shruthi Sathya Narayanan | 2026

---

## Overview

ATRT is a universally fatal pediatric CNS tumor driven almost exclusively by biallelic loss of SMARCB1 (INI1, ~95%) or SMARCA4 (~5%). Unlike DIPG, the driver is a loss-of-function event — not a gain-of-function mutation — which creates a synthetic lethality vulnerability not present in other brain tumors.

This pipeline integrates six independent data sources to systematically screen CNS/oncology drugs and identify combination hypotheses tailored to the SMARCB1-loss molecular context:

1. **Bulk RNA-seq** — GSE70678, 49 ATRT + normal samples (Torchia 2015, subgroup-annotated)
2. **CRISPR essentiality** — Broad Institute DepMap, ATRT/rhabdoid cell lines (BT16, BT37, G401, A204)
3. **Published IC50 validation** — 6 drugs, 10 ATRT cell lines (BT16, BT37, CHLA02/04/06, G401, A204)
4. **PPI network** — STRING-DB live API + curated fallback
5. **Transcriptomic reversal** — clue.io CMap L1000 query against ATRT SMARCB1-loss signature
6. **Patient genomics** — CBTN ATRT cohort (optional, controlled) or GSE70678 fallback

---

## Critical Biology: ATRT vs DIPG EZH2 Logic

**This is the most important difference from the DIPG pipeline:**

| | DIPG | ATRT |
|--|------|------|
| Driver | H3K27M gain-of-function | SMARCB1 loss-of-function |
| EZH2 status | Suppressed by H3K27M | Hyperactive (SMARCB1 antagonises PRC2) |
| EZH2 inhibitors | **PENALISED** (non-rational) | **BOOSTED ×1.40** (synthetic lethality) |
| Rationale | Bender 2014 Cancer Cell | Knutson 2013 PNAS PMID 23620515 |

Tazemetostat received FDA Breakthrough Therapy Designation specifically for SMARCB1-deficient tumors.

---

## Key Findings

### SMARCB1 Genomics

| Alteration | Prevalence | Source |
|-----------|------------|--------|
| SMARCB1 biallelic loss | ~95% | Hasselblatt 2011 Acta Neuropathol |
| SMARCA4 loss | ~5% | Frühwald 2020 CNS Oncol |

No co-occurrence test is performed. SMARCB1 loss IS the defining molecular event — there is no second alteration to test mutual exclusivity against.

### Three Molecular Subgroups (Johann 2016 Cancer Cell, n=150)

| Subgroup | Prevalence | Key Vulnerability | Key Genes |
|----------|------------|-------------------|-----------|
| ATRT-TYR | ~36% | HDAC + BET | TYR, DCT, MITF, HDAC1/2, BRD4 |
| ATRT-SHH | ~37% | EZH2 + SMO/GLI | GLI2, SMO, EZH2, CDK4 |
| ATRT-MYC | ~27% | AURKA + BET + EZH2 | MYC, MYCN, AURKA, BRD4 |

### Top Drug Combination Hypothesis

**Tazemetostat + Panobinostat + Alisertib**

| Drug | Score | BBB | Mechanism | Rationale |
|------|-------|-----|-----------|-----------|
| Tazemetostat | ~0.93 | HIGH | EZH2 inhibitor | SMARCB1 synthetic lethality (Knutson 2013) |
| Panobinostat | ~0.88 | HIGH | Pan-HDAC inhibitor | Epigenetic normalisation (Torchia 2015) |
| Alisertib | ~0.82 | HIGH | AURKA inhibitor | MYCN stabilisation (Sredni 2017) |

*Adjusted confidence range: [0.55, 0.70] (conservative / dose-optimised)*

### Top 8 Candidates

| Rank | Drug | Score | BBB | DepMap | Tissue | IC50 (µM) |
|------|------|-------|-----|--------|--------|-----------|
| 1 | TAZEMETOSTAT | ~0.930 | HIGH | 1.000 | 0.920 | 0.88 (G401) |
| 2 | PANOBINOSTAT | ~0.885 | HIGH | 1.000 | 0.929 | 0.009 (BT16) |
| 3 | ALISERTIB | ~0.820 | HIGH | 1.000 | 0.820 | 0.098 (BT16) |
| 4 | BIRABRESIB | ~0.813 | MODERATE | 1.000 | 0.880 | 0.31 (BT16) |
| 5 | ABEMACICLIB | ~0.795 | HIGH | 0.800 | 0.918 | 0.75 (BT16) |
| 6 | MARIZOMIB | ~0.780 | HIGH | 1.000 | 0.680 | 0.055 (BT16) |
| 7 | VISMODEGIB | ~0.650 | MODERATE | 0.500 | 0.650 | 2.50 (BT37) |
| 8 | ONC201 | ~0.640 | HIGH | 1.000 | 0.600 | 0.25 (BT16) |

*Note: Scores are approximate pending full data loading. Tazemetostat score boosted ×1.40 for EZH2 synthetic lethality.*

---

## Pipeline Architecture

```
save_results.py
├── data_fetcher.py              ← OpenTargets API (ATRT/rhabdoid EFO IDs)
├── discovery_pipeline.py        ← orchestrator + RNA loading
│   ├── tissue_expression.py     ← GSE70678 bulk RNA-seq scoring
│   │   └── run_quantile_sensitivity()  ← p65/p70/p75/p80 stability test
│   ├── depmap_essentiality.py   ← Broad CRISPR (ATRT/rhabdoid lines)
│   ├── escape_bypass            ← SMARCB1-loss resistance pathway scoring
│   ├── ppi_network.py           ← STRING-DB proximity
│   ├── bbb_filter.py            ← BBB penetrance + ATRT location-aware penalty
│   ├── cmap_query.py            ← clue.io pre-computed scores
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

**Weight rationale:** PPI weight = 0.05 because 88% of drugs score at floor (0.20) due to sparse curated neighbor coverage. DepMap weight = 0.35 to compensate with the strongest external signal. Top-4 ranking stable under ±10% weight perturbations.

---

## Data Sources

| Source | Description | n | Access |
|--------|-------------|---|--------|
| GSE70678 | Torchia 2015 ATRT bulk RNA-seq | 49 ATRT + normals | GEO (open) |
| GSE106982 | Johann 2016 methylation subgroups | 150 ATRT | GEO (open) |
| DepMap 23Q4 | Broad CRISPR essentiality | BT16, BT37, G401, A204 + others | depmap.org (open) |
| OpenTargets | Drug-target associations | ~557 drugs | API (open) |
| STRING-DB v11.5 | PPI network | genome-wide | API (open) |
| clue.io CMap L1000 | Transcriptomic reversal | 38,973 drug profiles | clue.io (free academic) |
| Published IC50 | ATRT cell-line validation | ~20 data points, 6 drugs | Literature (open) |
| CBTN ATRT | Patient CNA + mutation | controlled | portal.kidsfirstdrc.org |

---

## v1.0 Changes from DIPG v5.5

### Biological inversions
1. **EZH2 BOOST** — EZH2 inhibitors now receive ×1.40 composite score boost (vs ×0.50 penalty in DIPG). SMARCB1 loss removes PRC2 antagonist → EZH2 hyperactivity → synthetic lethality. Source: Knutson 2013 PNAS PMID 23620515.

2. **AURKA boost** — AURKA inhibitors boosted ×1.30 in MYC subgroup (MYCN stabilisation). Source: Sredni 2017 Pediatric Blood & Cancer PMID 28544500.

3. **SMO/GLI boost** — SMO inhibitors boosted ×1.25 in SHH subgroup (~37% of ATRT). Source: Torchia 2015 Cancer Cell PMID 26609405.

### Architecture changes
4. **No co-occurrence test** — DIPG's H3K27M/CDKN2A-del Fisher's exact test replaced with SMARCB1 prevalence validation. SMARCB1 loss IS the defining event — no two-hit co-occurrence hypothesis applies.

5. **GSE70678 bulk RNA** — replaces Filbin 2018 scRNA-seq (GSE102130). No ATRT scRNA-seq atlas exists as of 2026. GSE70678 (49 samples, Torchia 2015) is the primary tissue source.

6. **ATRT cell lines** — DepMap filter uses OncotreeSubtype "MRT"/"ATRT" and known names: BT16, BT37, G401, A204, BT12, CHLA02/04/06. G401/A204 are renal rhabdoid but share identical SMARCB1-loss biology and EZH2 dependency (Knutson 2013).

7. **BBB penalty less severe** — ATRT is not exclusively brainstem. Location distribution: infratentorial ~50%, supratentorial ~35%, spinal ~15% (Frühwald 2020). Unknown-location penalty: LOW=0.65×, UNKNOWN=0.90× (vs DIPG's 0.50×/0.85×).

8. **Subgroup stratification** — scoring can be pan-ATRT or stratified by TYR/SHH/MYC subgroup (Johann 2016).

---

## Known Limitations

1. **No scRNA-seq atlas** — no H3K27M-equivalent ATRT single-cell atlas exists as of 2026. GSE70678 is bulk RNA (49 samples) with less resolution than Filbin 2018 scRNA-seq.
2. **Small RNA reference** — GSE70678 n=49; escape bypass scores carry higher uncertainty than DepMap scores.
3. **Subgroup uncertainty** — without methylation array (EPIC 850K) or RNA-seq subgroup calling, pipeline uses pan-ATRT scores.
4. **No experimental validation** — hypothesis generation only; cell line and in vivo experiments required.
5. **PPI coverage** — 490/557+ drugs score at floor (0.20); PPI weight = 0.05 to reflect this.
6. **Tazemetostat CMap** — EZH2 inhibitors may not have been profiled in early L1000 datasets; CMap score uses neutral 0.50 prior if missing.

---

## Proposed Wet Lab Validation

1. **Cell viability** — BT16, BT37, G401 at 72h (CTG); triple combo vs individual drugs
2. **Synergy** — Bliss independence model, 5×5 dose matrix per drug pair
3. **Target engagement** — Western: H3K27me3↓ (tazemetostat), H3K27ac↑ (panobinostat), AURKA activity↓ (alisertib)
4. **SMARCB1 specificity** — compare SMARCB1-null vs SMARCB1-reconstituted isogenic lines
5. **EZH2 synthetic lethality** — confirm dependency is lost upon SMARCB1 re-expression

---

## Installation

```bash
git clone https://github.com/ShruthiSathya/gbmresearch.git
cd gbmresearch
pip install pandas numpy scipy matplotlib requests h5py
```

**Data placement:**
```
data/
  raw_omics/GSE70678_ATRT_expression.txt       # Torchia 2015 (GEO GSE70678)
  raw_omics/GSE106982_ATRT_methylation.txt      # Johann 2016 (optional)
  depmap/CRISPRGeneEffect.csv                   # Broad DepMap 23Q4
  depmap/Model.csv
  validation/cbtn_genomics/atrt/                # CBTN ATRT (optional, controlled)
  cmap_query/atrt_cmap_scores.json              # clue.io pre-computed scores
```

**Run:**
```bash
# Optional: prepare CMap query
python -m backend.pipeline.prepare_cmap_query
# [submit to clue.io, download query_result.gct]
python -m backend.pipeline.integrate_cmap_results

# Main pipeline
python -m backend.pipeline.save_results --disease atrt
python -m backend.pipeline.generate_figures
```

---

## Citation

Primary data: Torchia et al. 2015 (*Cancer Cell*), Johann et al. 2016 (*Cancer Cell*),
Knutson et al. 2013 (*PNAS*), Behan et al. 2019 (*Nature*),
Subramanian et al. 2017 (*Cell*) — CMap L1000.

Full reference list: see [REFERENCES.md](REFERENCES.md)

---

## License

MIT License. Data sources subject to their own access agreements.

---

*Shruthi Sathya Narayanan | [@ShruthiSathya](https://github.com/ShruthiSathya)*