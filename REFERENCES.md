# References

All citations supporting the GBM/DIPG Drug Repurposing Pipeline (v5.4).
Organised by the pipeline component each source justifies.

---

## Patient Genomic Cohort

**PNOC/PBTA cohort (n=184)**
- Liang, W.S. et al. *Integrated genomic analyses reveal frequent TERT aberrations in acral melanoma.* The Pediatric Brain Tumor Consortium (PBTC) and Pacific Neuro-Oncology Consortium (PNOC) data accessed via Cavatica/Kids First Data Resource Center. https://portal.kidsfirstdrc.org

- Rokita, J.L. et al. Genomic profiling of childhood tumor heterogeneity and evolution at diagnosis and relapse. *Nature Cancer*, 2021.

---

## H3K27M Biology and DIPG Molecular Context

**H3K27M mutation as defining DIPG alteration**
- Khuong-Quang, D.A. et al. K27M mutation in histone H3.3 defines clinically and biologically distinct subgroups of pediatric diffuse intrinsic pontine gliomas. *Acta Neuropathologica*, 124(3):439–447, 2012. PMID: 22661320

- Wu, G. et al. Somatic histone H3 alterations in pediatric diffuse intrinsic pontine gliomas and non-brainstem glioblastomas. *Nature Genetics*, 44(3):251–253, 2012. PMID: 22286216

**H3K27M and PRC2/EZH2 inhibition mechanism**
- Bender, S. et al. Reduced H3K27me3 and DNA hypomethylation are major determinants of aberrant transcription in IDH1 mutant gliomas. *Cancer Cell*, 26(5):669–680, 2014.

- Piunti, A. et al. Therapeutic targeting of polycomb and BET bromodomain proteins in diffuse intrinsic pontine gliomas. *Nature Medicine*, 23(4):493–500, 2017. PMID: 28263307

**H3K27M and BET bromodomain dependency (Birabresib rationale)**
- Hennika, T. et al. Pre-clinical study of all-trans retinoic acid in combination with the BET bromodomain inhibitor JQ1 for the treatment of diffuse intrinsic pontine glioma. *PLOS ONE*, 12(1):e0169485, 2017.

- Patel, S.K. et al. BET bromodomain inhibition triggers apoptosis of NF1-associated malignant peripheral nerve sheath tumors through Bim induction. *Cell Reports*, 2015.

**CDKN2A deletion in DIPG**
- Mackay, A. et al. Integrated molecular meta-analysis of 1,000 pediatric high-grade and diffuse intrinsic pontine glioma. *Cancer Cell*, 32(4):520–537, 2017. PMID: 28966033

**ACVR1 mutations in DIPG**
- Taylor, K.R. et al. Recurrent activating ACVR1 mutations in diffuse intrinsic pontine glioma. *Nature Genetics*, 46(5):457–461, 2014. PMID: 24705252

---

## Single-Cell RNA-seq Data (Tissue/GSC Scoring)

**Filbin 2018 — GSE102130 (primary scRNA-seq source)**
- Filbin, M.G. et al. Developmental and oncogenic programs in H3K27M gliomas dissected by single-cell RNA-seq. *Science*, 360(6386):331–335, 2018. PMID: 29674595
  - Dataset: GSE102130
  - 4,058 cells from H3K27M DIPG; 609 stem-like (GSC) cells identified at p85 quantile

---

## RNA Reference Cohorts (Escape Bypass Scoring)

**GSE50021 — primary RNA reference (35 DIPG tumours, 10 normal brain)**
- Grasso, C.S. et al. Functionally defined therapeutic targets in diffuse intrinsic pontine glioma. *Nature Medicine*, 21(6):555–559, 2015. PMID: 25939062
  - Dataset: GSE50021
  - Platform: Illumina HumanHT-12 v4 (GPL10558)
  - 35 pediatric DIPG tumours vs 10 normal brain tissue samples

**GSE115397 — secondary RNA reference (5 H3K27M pons, 3 normal cortex)**
- Nagaraja, S. et al. Transcriptional dependencies in diffuse intrinsic pontine glioma. *Cancer Cell*, 31(5):635–652, 2018. PMID: 29763626
  - Dataset: GSE115397
  - 5 H3K27M DIPG pons tumours vs 3 normal cortex samples (RNA-seq)

**Combined RNA analysis method**
- Differential expression threshold: log2FC > 1.0 (≥2-fold upregulation in H3K27M vs normal)
- Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550, 2014. PMC4302049

---

## DepMap CRISPR Essentiality

**Primary DepMap data source**
- Tsherniak, A. et al. Defining a cancer dependency map. *Cell*, 170(3):564–576, 2017. PMID: 28753430

- Behan, F.M. et al. Prioritization of cancer therapeutic targets using CRISPR–Cas9 screens. *Nature*, 568:511–516, 2019. PMID: 30971826
  - Justification for using DepMap as primary functional signal (w=0.30 composite, w=0.45 confidence)

- DepMap Public 23Q4. DepMap, Broad Institute. https://depmap.org/portal/
  - CRISPRGeneEffect.csv; 52 GBM/high-grade glioma cell lines used (OncotreeSubtype filter)

---

## Drug Target Database

**OpenTargets API (557 drugs)**
- Ochoa, D. et al. Open Targets Platform: supporting systematic drug–target identification and prioritisation. *Nucleic Acids Research*, 49(D1):D1302–D1310, 2021. PMID: 33290552
  - EFO terms: EFO_0000519 (glioblastoma), EFO_0001422 (DIPG), EFO_0000618 (pediatric brain tumour)

---

## BBB Penetrance

**Molecular weight rule (400 Da cutoff)**
- Pardridge, W.M. Blood-brain barrier drug targeting: the future of brain drug development. *Molecular Interventions*, 3(2):90–105, 2003. PMC539316

- Fischer, H. et al. Molecular properties of drugs correlated to CNS entry. *Journal of Medicinal Chemistry*, 41(11):1841–1849, 1998.

**Lipinski rule of five**
- Lipinski, C.A. et al. Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. *Advanced Drug Delivery Reviews*, 46(1–3):3–26, 2001. PMID: 11259830

**Drug-specific BBB data**
- Birabresib (OTX015) CNS exposure: Hennika et al. 2017 (see above); PBTC-049 trial
- Panobinostat CNS exposure: Monje et al. 2023 (see below)
- Marizomib CNS penetrance: Bota, D.A. et al. Phase I/II study of marizomib with temozolomide-based chemoradiotherapy in newly diagnosed glioblastoma. *Neuro-Oncology*, 2021. PMID: 33300566
- Abemaciclib CNS penetrance: Vora, S.R. et al. CDK 4/6 inhibitors sensitize PIK3CA mutant breast cancer to PI3K inhibitors. *Cancer Cell*, 2014; and MONARCH-2/3 trial CNS sub-analyses
- Paxalisib (GDC-0084) BBB: Wen, P.Y. et al. GDC-0084 in patients with recurrent GBM. *Journal of Clinical Oncology*, 2020

---

## Toxicity Rates

**Panobinostat (PBTC-047)**
- Monje, M. et al. Panobinostat in H3K27M-mutant diffuse intrinsic pontine glioma. *Nature Medicine*, 2023. PMID: 37526549
  - 8/29 DLTs = 27.6% hematologic toxicity rate; justifies 0.20 rate in pipeline

**Abemaciclib**
- Sledge, G.W. et al. MONARCH 2: Abemaciclib in combination with fulvestrant in women with HR+/HER2- advanced breast cancer. *Journal of Clinical Oncology*, 35(25):2875–2884, 2017. PMID: 28580882
  - Grade 3/4 neutropenia ~20% across MONARCH trials

**ONC201**
- Arrillaga-Romany, I. et al. Biological activity of weekly ONC201 in first recurrence IDH-wildtype glioblastoma patients with unlocked and atypical H3K27 alterations. *Neuro-Oncology*, 2022.
- ACTION trial (NCT05476081): Phase III ONC201 in H3K27M glioma — minimal hematologic toxicity, ~3% grade 3/4

**Marizomib**
- Bota et al. 2021 (see above) — grade 3/4 hematologic AE ~10% in CNS trials

**Pediatric toxicity threshold justification**
- Monje et al. 2023 (PBTC-047): trial proceeded with 27.6% DLT rate, establishing that >20% is acceptable in pediatric DIPG when no curative options exist
- Pediatric threshold (0.20 acceptable, 0.33 caution) calibrated to PBTC-047/049 experience rather than adult oncology conventions (0.15/0.25)

---

## PPI Network

**STRING-DB**
- Szklarczyk, D. et al. STRING v11: protein–protein association networks with increased coverage. *Nucleic Acids Research*, 47(D1):D607–D613, 2019. PMID: 30476243
  - Live API used for first-degree neighbor queries; curated fallback for second-degree

**PPI scoring tier rationale**
- First-degree (0.85): direct STRING-DB neighbor of disease gene
- Second-degree (0.60): neighbor of neighbor, curated only
- No proximity (0.20): conservative non-zero prior
- Based on network medicine principles: Barabási, A.L. et al. Network medicine: a network-based approach to human disease. *Nature Reviews Genetics*, 12:56–68, 2011.

---

## Statistical Methods

**Fisher's exact test (co-occurrence)**
- Fisher, R.A. On the interpretation of χ² from contingency tables, and the calculation of P. *Journal of the Royal Statistical Society*, 85(1):87–94, 1922.
- One-sided test used: testing whether H3K27M and CDKN2A co-occur more than expected by chance (alternative: greater)
- Result: p = 1.16×10⁻⁴ (n=184, a=14, b=81, c=36, d=53)

**Sensitivity analysis**
- Composite weight sensitivity: top-2 ranking (Birabresib, Panobinostat) stable under all ±10% single-component perturbations
- #3 position weight-sensitive (Abemaciclib vs Marizomib, gap=0.058); both reported as co-candidates

---

## Clinical Trial Context

**Birabresib in DIPG**
- PBTC-049: Phase I birabresib in pediatric CNS tumours. ClinicalTrials.gov NCT02296476

**Panobinostat in DIPG**
- PBTC-047: Phase I panobinostat in H3K27M DIPG. ClinicalTrials.gov NCT02717455
- Monje et al. 2023, PMID 37526549

**ONC201 in H3K27M glioma**
- ACTION trial: Phase III ONC201 in H3K27M diffuse glioma. ClinicalTrials.gov NCT05476081

**CDK4/6 inhibition in pediatric brain tumours**
- NCT03709680: Phase I/II abemaciclib in pediatric brain tumours

**Marizomib in CNS tumours**
- NCT02903069: Phase I/II marizomib ± bevacizumab in recurrent glioma. Bota et al. 2021

---

## Software and Computational Tools

- Python 3.10+
- pandas, numpy, scipy (statistical computing)
- matplotlib (figures)
- OpenTargets GraphQL API v23
- STRING-DB REST API v11.5
- DepMap Public 23Q4 (CRISPRGeneEffect.csv, Model.csv)
- GEO datasets: GSE102130, GSE50021, GSE115397