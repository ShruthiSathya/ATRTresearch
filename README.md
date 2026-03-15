# GBM/DIPG Drug Repurposing Pipeline

**A multi-omic computational pipeline for prioritising drug combinations in H3K27M-mutant Diffuse Intrinsic Pontine Glioma (DIPG)**

> Developed by Shruthi Sathya Narayanan | March 2026 

---

## Overview

DIPG is a universally fatal pediatric brainstem cancer with no curative treatment. The H3K27M histone mutation, present in ~80% of cases, fundamentally rewires the epigenetic landscape and creates therapeutic vulnerabilities not present in adult GBM.

This pipeline integrates six independent data sources to systematically screen 557 CNS/oncology drugs and identify combination hypotheses tailored to the H3K27M molecular context:

1. **Patient genomics** — 184 PNOC/PBTA DIPG samples (CNA + mutation data)
2. **Single-cell RNA-seq** — Filbin 2018 H3K27M scRNA-seq, 609 stem-like GSC cells
3. **CRISPR essentiality** — Broad Institute DepMap, 52 GBM cell lines
4. **RNA reference cohorts** — GSE50021 (35 DIPG) + GSE115397 (5 H3K27M pons)
5. **PPI network** — STRING-DB live API + curated fallback
6. **Transcriptomic reversal** — clue.io CMap L1000 query against H3K27M signature (38,973 drug profiles)

---

## Key Findings

### Genomic Co-occurrence (n=184 PNOC/PBTA patients)

| Alteration | n | % |
|-----------|---|---|
| H3K27M mutation | 95 | 52% |
| CDKN2A deletion | 50 | 27% |
| Double-hit (both) | 14 | 7.6% |

Fisher's exact test (one-sided): **p = 7.55×10⁻⁵** — H3K27M and CDKN2A deletion co-occur significantly less than expected by chance, indicating mutual exclusivity and alternative oncogenic mechanisms.

### Top Drug Combination Hypothesis

**Birabresib + Panobinostat + Marizomib**

| Drug | Score | BBB | Mechanism | CMap norm_cs |
|------|-------|-----|-----------|--------------|
| Birabresib | 0.883 | MODERATE | BET bromodomain inhibitor (BRD4) | -1.08 (reverser) |
| Panobinostat | 0.855 | HIGH | Pan-HDAC inhibitor | -1.17 (reverser) |
| Marizomib | 0.823 | HIGH | Proteasome inhibitor (PSMB5) | N/A (not in L1000) |

*Adjusted confidence: [0.55, 0.70] (raw 0.92 × toxicity multiplier 0.60 conservative / 0.76 dose-optimised — HIGH_RISK: 40% combined hematologic AE)*
*Toxicity flag: HIGH_RISK (40% combined hematologic AE rate, additive model)*
*Confidence is reported as a range: conservative (additive toxicity model) to optimistic (dose-optimised, 60% AE reduction based on PBTC-047 precedent)*

### Top 8 Candidates (557 drugs screened)

| Rank | Drug | Score | BBB | DepMap | Tissue | CMap | IC50 (µM) |
|------|------|-------|-----|--------|--------|------|-----------|
| 1 | BIRABRESIB | 0.883 | MODERATE | 1.000 | 0.937 | -1.08 | 0.18 (SU-DIPG-IV) |
| 2 | PANOBINOSTAT | 0.855 | HIGH | 1.000 | 0.929 | -1.17 | 0.005 (DIPG4) |
| 3 | MARIZOMIB | 0.823 | HIGH | 1.000 | 0.770 | N/A | 0.032 (SU-DIPG-IV) |
| 4 | ABEMACICLIB | 0.813 | HIGH | 0.800 | 0.918 | -1.51 | 0.52 (SU-DIPG-IV) |
| 5 | INDOXIMOD | 0.795 | MODERATE | 1.000 | 0.627 | -1.54 | — |
| 6 | ONATASERTIB | 0.760 | MODERATE | 1.000 | 0.670 | N/A | — |
| 7 | REGORAFENIB | 0.675 | MODERATE | 0.500 | 0.886 | N/A | — |
| 8 | DOVITINIB | 0.675 | MODERATE | 0.500 | 0.876 | N/A | — |

*Note: Drugs with BBB=LOW (Crizotinib, AZD-8055) receive a 0.50× DIPG brainstem penetrance penalty and are excluded from top candidates.*

---

## Pipeline Architecture

```
save_results.py
├── data_fetcher.py              ← OpenTargets API (557 drugs)
├── discovery_pipeline.py        ← orchestrator + RNA loading
│   ├── tissue_expression.py     ← scRNA-seq + curated DIPG scoring
│   │   └── run_quantile_sensitivity()  ← GSC p75/p80/p85/p90 stability test
│   ├── depmap_essentiality.py   ← Broad CRISPR (52 GBM lines)
│   ├── escape_bypass            ← RNA-confirmed resistance pathway scoring
│   ├── ppi_network.py           ← STRING-DB proximity (extended proteasome neighbors)
│   ├── bbb_filter.py            ← BBB penetrance + DIPG brainstem penalty
│   ├── cmap_query.py            ← clue.io pre-computed scores (38,973 drugs)
│   ├── published_ic50_validation.py  ← Published DIPG cell-line IC50 anchoring
│   └── hypothesis_generator.py  ← combination + toxicity CI range
├── dipg_specialization.py       ← H3K27M/ACVR1 scoring + EZH2 inhibitor penalty
├── pipeline_config.py           ← ALL parameters (single source of truth)
└── generate_figures.py          ← 4 publication figures
```

### Composite Scoring Formula

```
composite_score = 0.40 × tissue_expression      (Filbin 2018 scRNA-seq, w=0.40)
                + 0.35 × depmap_essentiality     (Broad CRISPR, Behan 2019, w=0.35)
                + 0.20 × escape_bypass           (RNA-informed resistance, w=0.20)
                + 0.05 × ppi_network             (STRING-DB proximity, w=0.05)

confidence_adj  = (0.45 × depmap + 0.35 × bbb + 0.20 × diversity)
                × toxicity_multiplier

toxicity_CI     = [raw × multiplier_conservative, raw × multiplier_optimistic]
                = [raw × (1 - combined_AE_rate), raw × (1 - combined_AE_rate × 0.60)]
```

**Weight changes from v5.4:** DepMap increased from 0.30 → 0.35; PPI reduced from 0.10 → 0.05. Rationale: PPI scored 490/557 drugs at the floor value (0.20) due to sparse curated neighbor coverage, contributing near-zero discriminative signal for 88% of candidates. DepMap weight increased to compensate with the strongest external signal.

Top-2 ranking stable under all ±10% weight perturbations (sensitivity analysis). Top-4 stable under DIPG BBB penalty.

---

## v5.5 Changes (March 2026)

### Accuracy fixes
1. **EZH2 inhibitor penalty** — EZH2 inhibitors now receive 0.50× composite score penalty in H3K27M DIPG. H3K27M dominant-negatively inhibits PRC2/EZH2 (Bender et al. 2014); EZH2 inhibitors are mechanistically non-rational. Previously the specialization module was boosting these drugs while the tissue scorer was penalising them — a direct contradiction now resolved.

2. **Marizomib PPI score fixed** — PSMB5/PSMB2/PSMB1 now have curated neighbors (TP53, MYC, H3F3A, MCL1) based on published proteasome biology (MYC t½ ~20 min; Lin et al. 2019). Marizomib PPI score: 0.20 → 0.85.

3. **UNKNOWN BBB score** — Changed from 0.50 to 0.40. Unknown ≠ moderate; lack of PK data should carry a mild penalty.

4. **BBB extended known list** — AZD-8055 (LOW), crizotinib (LOW), pazopanib (LOW), regorafenib (MODERATE), nintedanib (LOW), and others filled from published PK literature. Previously these were UNKNOWN and received an inflated 0.50 score.

5. **DIPG brainstem BBB penalty** — LOW BBB drugs receive 0.50× score penalty; UNKNOWN drugs receive 0.85× penalty. Brainstem BBB is tighter than cortex. Crizotinib and AZD-8055 correctly excluded from top candidates.

6. **PPI weight** — Reduced from 0.10 to 0.05; DepMap increased from 0.30 to 0.35.

7. **Escape bypass** — Now RNA-confirmed when ≥10 upregulated genes available from patient RNA data. Previously used curated fallback exclusively (less precise).

8. **Toxicity confidence range** — Adjusted confidence now reported as [conservative, optimistic] range rather than a single point estimate. Additive toxicity model is conservative; dose-optimised estimate accounts for PBTC-047 precedent.

9. **GSC quantile sensitivity** — `run_quantile_sensitivity()` tests p75/p80/p85/p90 thresholds and reports top-2 ranking stability.

### New data streams
10. **CMap transcriptomic reversal** — clue.io L1000 query against 1,513 H3K27M upregulated genes. 38,973 drug profiles integrated as pre-computed norm_cs scores. All testable top candidates confirmed as H3K27M signature reversers (norm_cs < -0.9, FDR -log10 = 15.65).

11. **Published IC50 validation** — `published_ic50_validation.py` anchors tissue scores against 18 published DIPG cell-line data points across 6 drugs (Birabresib, Panobinostat, Marizomib, Abemaciclib, ONC201, Paxalisib). Zero discordances between computational scores and experimental IC50 data.

12. **ACVR1 direct mutation calling** — Pipeline now attempts direct ACVR1 mutation calling (R206H, G328V, G328E, G356D, R258G) from mutations.txt. Falls back to 25% prevalence prior when not available.

---

## CMap Validation Results

All testable top pipeline candidates confirmed as H3K27M signature reversers:

| Drug | CMap name | norm_cs | Reverser |
|------|-----------|---------|----------|
| Birabresib | OTX-015 | -1.08 | ✓ |
| Panobinostat | PANOBINOSTAT | -1.17 | ✓ |
| Abemaciclib | CDK4 (shRNA proxy) | -1.51 | ✓ |
| ONC201 | TIC-10 | -1.03 | ✓ |
| AZD-8055 | AZD-8055 | -1.59 | ✓ |
| Indoximod | INDOXIMOD | -1.54 | ✓ |
| Paxalisib | GDC-0941 (PI3K proxy) | -1.32 | ✓ |
| Marizomib | — | N/A | Not in L1000 |

*norm_cs < -0.9 = strong reversal of H3K27M disease signature. FDR -log10 = 15.65 for all significant hits.*

---

## Data Sources

| Source | Description | n | Access |
|--------|-------------|---|--------|
| PNOC/PBTA | Patient CNA + mutation | 184 samples | Cavatica (controlled) |
| GSE102130 | Filbin 2018 H3K27M scRNA-seq | 4,058 cells | GEO (open) |
| GSE50021 | Grasso 2015 DIPG microarray | 35 DIPG + 10 normal | GEO (open) |
| GSE115397 | Nagaraja 2018 RNA-seq | 5 DIPG + 3 normal | GEO (open) |
| DepMap 23Q4 | Broad CRISPR essentiality | 52 GBM lines | depmap.org (open) |
| OpenTargets | Drug-target associations | 557 drugs | API (open) |
| STRING-DB v11.5 | PPI network | genome-wide | API (open) |
| clue.io CMap L1000 | Transcriptomic reversal | 38,973 drug profiles | clue.io (free academic) |
| Published IC50 | DIPG cell-line validation | 18 data points, 6 drugs | Literature (open) |

---

## Known Limitations

1. **Cohort mismatch** — genomics (PNOC/PBTA) and RNA reference (GSE50021 + GSE115397) are separate cohorts; cross-cohort use justified by shared H3K27M DIPG context
2. **No experimental validation** — hypothesis generation only; cell line and in vivo experiments required
3. **ACVR1 subgroup** — direct mutation calling unavailable from current data export; estimated n=24 using 25% prevalence prior (Taylor 2014). Full MAF from CBTTC would enable direct calling.
4. **Composite weights** — rationale-based hierarchy, not empirically derived; sensitivity analysis confirms top-4 stability
5. **CMap RNA reference** — n=5 H3K27M samples (GSE115397) is small; escape bypass scores carry higher uncertainty than DepMap/tissue scores
6. **Marizomib CMap** — not profiled in L1000 (marine natural product, profiled after dataset generation); CMap score uses neutral 0.50 prior
7. **PPI coverage** — 490/557 drugs score at floor (0.20) due to sparse curated neighbor coverage; PPI weight reduced to 0.05 to reflect this

---

## Proposed Wet Lab Validation

Minimum experiment set to test the top hypothesis:

1. **Cell viability** — SU-DIPG-IV, SU-DIPG-XIII at 72h (CTG assay); triple combo vs individual drugs
2. **Synergy** — Bliss independence model, 5×5 dose matrix per drug pair
3. **Target engagement** — Western: BRD4↓ (Birabresib), H3K27ac↑ (Panobinostat), 20S activity↓ (Marizomib)
4. **H3K27M specificity** — compare H3K27M+ vs H3K27M− isogenic lines
5. **CMap validation** — RNA-seq after drug treatment; confirm reversal of H3K27M signature

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
  raw_omics/GSE102130_K27Mproject.RSEM.vh20170621.txt   # Filbin scRNA-seq
  depmap/CRISPRGeneEffect.csv                            # Broad DepMap
  depmap/Model.csv
  validation/cbtn_genomics/mutations.txt                 # PNOC/PBTA (controlled)
  validation/cbtn_genomics/cna.txt
  validation/cbtn_genomics/rna_zscores.txt               # Processed RNA
  cmap_query/cmap_scores_pipeline.json                   # clue.io pre-computed scores
```

**Run:**
```bash
# Optional: prepare and run clue.io query (one-time)
python -m backend.pipeline.prepare_cmap_query
# [submit query at clue.io, download arfs/TAG/query_result.gct]
python -m backend.pipeline.integrate_cmap_results

# Main pipeline
python -m backend.pipeline.save_results
python -m backend.pipeline.generate_figures
```

---

## Citation

Primary data citations: Filbin et al. 2018 (*Science*), Grasso et al. 2015 (*Nature Medicine*),
Behan et al. 2019 (*Nature*), Monje et al. 2023 (*Nature Medicine*),
Subramanian et al. 2017 (*Cell*) — CMap L1000.

Full reference list: see [REFERENCES.md](REFERENCES.md)

---

## License

MIT License. Data sources subject to their own access agreements.
clue.io CMap data subject to Broad Institute terms of use.

---

*Shruthi Sathya Narayanan | [@ShruthiSathya](https://github.com/ShruthiSathya)*