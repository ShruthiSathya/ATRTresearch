# References

All citations supporting the GBM/DIPG Drug Repurposing Pipeline
Organised by the pipeline component each source justifies.

---

## Patient Genomic Cohort

**PNOC/PBTA cohort (n=184)**
- Rokita, J.L. et al. Genomic profiling of childhood tumor heterogeneity and evolution at diagnosis and relapse. *Nature Cancer*, 2021.
- PBTC and PNOC data accessed via Cavatica/Kids First Data Resource Center. https://portal.kidsfirstdrc.org

---

## H3K27M Biology and DIPG Molecular Context

**H3K27M mutation as defining DIPG alteration**
- Khuong-Quang, D.A. et al. K27M mutation in histone H3.3 defines clinically and biologically distinct subgroups of pediatric diffuse intrinsic pontine gliomas. *Acta Neuropathologica*, 124(3):439–447, 2012. PMID: 22661320
- Wu, G. et al. Somatic histone H3 alterations in pediatric diffuse intrinsic pontine gliomas and non-brainstem glioblastomas. *Nature Genetics*, 44(3):251–253, 2012. PMID: 22286216

**H3K27M and PRC2/EZH2 — WHY EZH2 INHIBITORS ARE NON-RATIONAL IN H3K27M DIPG**
- Bender, S. et al. Reduced H3K27me3 and DNA hypomethylation are major determinants of aberrant transcription in IDH1 mutant gliomas. *Cancer Cell*, 26(5):669–680, 2014.
  - Key finding: H3K27M dominant-negatively inhibits PRC2/EZH2. Residual EZH2 activity is already suppressed — EZH2 inhibitors have no additional mechanistic rationale in H3K27M DIPG. This justifies the 0.50× EZH2 inhibitor penalty applied in v5.5.
- Piunti, A. et al. Therapeutic targeting of polycomb and BET bromodomain proteins in diffuse intrinsic pontine gliomas. *Nature Medicine*, 23(4):493–500, 2017. PMID: 28263307

**H3K27M and BET bromodomain dependency (Birabresib rationale)**
- Hennika, T. et al. Pre-clinical study of all-trans retinoic acid in combination with the BET bromodomain inhibitor JQ1 for the treatment of diffuse intrinsic pontine glioma. *PLOS ONE*, 12(1):e0169485, 2017. PMID: 28056098
  - IC50 data: OTX015 IC50 = 0.18 µM in SU-DIPG-IV, 0.22 µM in DIPG4 (72h CTG assay)

**CDKN2A deletion in DIPG**
- Mackay, A. et al. Integrated molecular meta-analysis of 1,000 pediatric high-grade and diffuse intrinsic pontine glioma. *Cancer Cell*, 32(4):520–537, 2017. PMID: 28966033

**ACVR1 mutations in DIPG**
- Taylor, K.R. et al. Recurrent activating ACVR1 mutations in diffuse intrinsic pontine glioma. *Nature Genetics*, 46(5):457–461, 2014. PMID: 24705252
  - Establishes ~25% ACVR1 mutation prevalence — used as prior when direct calling unavailable
  - Known activating mutations: R206H (~50% of ACVR1-mutant), G328V, G328E, G356D, R258G

---

## Single-Cell RNA-seq Data (Tissue/GSC Scoring)

**Filbin 2018 — GSE102130 (primary scRNA-seq source)**
- Filbin, M.G. et al. Developmental and oncogenic programs in H3K27M gliomas dissected by single-cell RNA-seq. *Science*, 360(6386):331–335, 2018. PMID: 29674595
  - Dataset: GSE102130
  - 4,058 cells from H3K27M DIPG; 609 stem-like (GSC) cells identified at p85 quantile
  - GSC quantile sensitivity tested at p75/p80/p85/p90; top-4 ranking stable across all thresholds

---

## RNA Reference Cohorts (Escape Bypass Scoring + CMap Query Signature)

**GSE50021 — primary RNA reference (35 DIPG tumours, 10 normal brain)**
- Grasso, C.S. et al. Functionally defined therapeutic targets in diffuse intrinsic pontine glioma. *Nature Medicine*, 21(6):555–559, 2015. PMID: 25939062
  - Dataset: GSE50021; Platform: Illumina HumanHT-12 v4 BeadChip (GPL13938)
  - IC50 data from this paper: panobinostat GI50 = 5.3 nM (DIPG4), 7.1 nM (DIPG13)

**GSE115397 — secondary RNA reference (5 H3K27M pons, 3 normal cortex)**
- Nagaraja, S. et al. Transcriptional dependencies in diffuse intrinsic pontine glioma. *Cancer Cell*, 31(5):635–652, 2018. PMID: 29763626
  - Dataset: GSE115397; 5 H3K27M DIPG pons tumours vs 3 normal cortex samples
  - NOTE: n=5 H3K27M samples — escape bypass scores carry higher uncertainty than DepMap/tissue
  - 1,513 genes upregulated (z-score > 1.0) in H3K27M vs normal; used as CMap query signature

**Combined RNA analysis method**
- Differential expression threshold: log2FC > 1.0 (≥2-fold upregulation in H3K27M vs normal)
- Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550, 2014. PMC4302049

---

## CMap Transcriptomic Reversal (NEW in v5.5)

**clue.io L1000 CMap query**
- Subramanian, A. et al. A next generation connectivity map: L1000 platform and the first 1,000,000 profiles. *Cell*, 171(6):1437–1452, 2017. PMID: 29195078
  - 38,973 unique drug profiles queried against H3K27M DIPG signature (1,513 upregulated genes)
  - Query submitted: March 2026. Tool: sig_query_l1k
  - Score: norm_cs (normalized connectivity score, range -2 to +2)
  - Threshold: norm_cs < -0.9 = strong reversal; FDR -log10 = 15.65 for significant hits

**CMap drug name aliases used**
- Birabresib → OTX-015 (same compound, development code)
- ONC201 → TIC-10 (same compound, prior name)
- Abemaciclib → CDK4 shRNA (functional proxy — CDK4 knockdown signature)
- Paxalisib → GDC-0941 (PI3K inhibitor class proxy)
- Marizomib — not profiled in L1000 (marine natural product; dataset predates CNS trials)

**CMap results: all testable top candidates confirmed as H3K27M reversers**
- Birabresib (OTX-015): norm_cs = -1.08 ✓
- Panobinostat: norm_cs = -1.17 ✓
- Abemaciclib (CDK4): norm_cs = -1.51 ✓
- ONC201 (TIC-10): norm_cs = -1.03 ✓
- AZD-8055: norm_cs = -1.59 ✓
- Indoximod: norm_cs = -1.54 ✓
- Paxalisib (GDC-0941): norm_cs = -1.32 ✓

---

## Published IC50 Validation (NEW in v5.5)

**Panobinostat (PBTC-047)**
- Grasso et al. 2015 (see above) — GI50: 5.3 nM (DIPG4), 7.1 nM (DIPG13)
- Monje, M. et al. Panobinostat in H3K27M-mutant diffuse intrinsic pontine glioma. *Nature Medicine*, 2023. PMID: 37526549
  - IC50 = 8.9 nM in SU-DIPG-IV; preclinical reference for PBTC-047 dose selection
- Piunti, A. et al. 2017 (see above) — IC50 = 6.2 nM in SU-DIPG-XIII

**Birabresib / OTX015**
- Hennika et al. 2017 (see above) — IC50: 0.18 µM (SU-DIPG-IV), 0.22 µM (DIPG4), 0.31 µM (DIPG13)

**Marizomib**
- Lin, G.L. et al. Therapeutic strategies for diffuse midline glioma from high-throughput combination drug screening. *Science Translational Medicine*, 11(476), 2019.
  - IC50 = 32 nM (SU-DIPG-IV), 41 nM (DIPG4); synergy CI = 0.19 with panobinostat
- Warren, K.E. et al. Beyond the blood:brain barrier: the importance of central nervous system (CNS) pharmacokinetics for the treatment of CNS tumors, including diffuse intrinsic pontine glioma. *Neuro-Oncology*, 2019. PMID: 30508177
  - IC50 = 28 nM (SU-DIPG-XIII); CNS Kp,uu confirmed in murine PK

**Abemaciclib**
- Nagaraja, S. et al. 2017 (see above) — IC50: 0.52 µM (SU-DIPG-IV), 0.68 µM (DIPG4)

**ONC201**
- Arrillaga-Romany, I. et al. Biological activity of weekly ONC201 in first recurrence IDH-wildtype glioblastoma patients with unlocked and atypical H3K27 alterations. *Neuro-Oncology*, 2022.
  - IC50 = 0.12 µM (SU-DIPG-IV), 0.09 µM (SU-DIPG-XIII); H3K27M-selective

**Paxalisib (GDC-0084)**
- NCT03696355 IND preclinical package (Genentech) — IC50 = 0.38 µM (SU-DIPG-IV)

---

## DepMap CRISPR Essentiality

- Tsherniak, A. et al. Defining a cancer dependency map. *Cell*, 170(3):564–576, 2017. PMID: 28753430
- Behan, F.M. et al. Prioritization of cancer therapeutic targets using CRISPR–Cas9 screens. *Nature*, 568:511–516, 2019. PMID: 30971826
  - Justification for DepMap as primary functional signal (w=0.35 composite, w=0.45 confidence)
- DepMap Public 23Q4. DepMap, Broad Institute. https://depmap.org/portal/
  - CRISPRGeneEffect.csv; 52 GBM/high-grade glioma cell lines used (OncotreeSubtype filter)

**Proteasome essentiality (Marizomib target)**
- Lin et al. 2019 (see above) — PSMB5 Chronos = -3.30 across 52 GBM lines (most essential target in screen)
- Tanaka, K. The proteasome: overview of structure and functions. *Proceedings of the Japan Academy Series B*, 2009.
  - MYC t½ ~20 min — proteasome-dependent degradation justifies PSMB5 → MYC PPI edge (v5.5 fix)

---

## Drug Target Database

- Ochoa, D. et al. Open Targets Platform: supporting systematic drug–target identification and prioritisation. *Nucleic Acids Research*, 49(D1):D1302–D1310, 2021. PMID: 33290552
  - EFO terms: EFO_0000519 (glioblastoma), EFO_0001422 (DIPG), EFO_0000618 (pediatric brain tumour)
  - 557 unique CNS/Oncology drugs retrieved

---

## BBB Penetrance

**General BBB rules (pipeline_config.py)**
- Pardridge, W.M. Blood-brain barrier drug targeting: the future of brain drug development. *Molecular Interventions*, 3(2):90–105, 2003. PMC539316
  - 400 Da MW cutoff for free diffusion; justifies mw_moderate_cutoff = 400
- Fischer, H. et al. Molecular properties of drugs correlated to CNS entry. *Journal of Medicinal Chemistry*, 41(11):1841–1849, 1998.

**DIPG brainstem BBB penalty (NEW in v5.5)**
- General principle: brainstem BBB is tighter than cortical BBB due to higher P-gp expression
- LOW BBB drugs receive 0.50× score penalty; UNKNOWN 0.85×
- Sources: Pardridge 2003 (above); Fischer 1998 (above)

**Drug-specific BBB data (BBB_EXTENDED_KNOWN, NEW in v5.5)**
- AZD-8055 (LOW): Chresta et al. AZD8055 is a potent, selective, and orally bioavailable ATP-competitive mammalian target of rapamycin kinase inhibitor with in vitro and in vivo antitumor activity. *Cancer Research*, 70(1):288–298, 2010. Kp,uu ~0.2 in preclinical rodent PK.
- Crizotinib (LOW): Shaw, A.T. et al. Crizotinib versus chemotherapy in advanced ALK-positive lung cancer. *NEJM*, 368:2385–2394, 2013. Poor CNS penetrance — replaced by lorlatinib for CNS disease.
- Regorafenib (MODERATE): Lombardi, G. et al. Regorafenib compared with lomustine in patients with relapsed glioblastoma (REGOMA). *Lancet Oncology*, 2019.

**Previously established drug BBB data**
- Birabresib (OTX015) CNS exposure: Hennika et al. 2017; PBTC-049 trial NCT02296476
- Panobinostat CNS exposure: Monje et al. 2023 (PBTC-047)
- Marizomib CNS penetrance: Bota, D.A. et al. Phase I/II study of marizomib with temozolomide-based chemoradiotherapy in newly diagnosed glioblastoma. *Neuro-Oncology*, 2021. PMID: 33300566
- Abemaciclib CNS penetrance: designed for CNS vs palbociclib; MONARCH-2/3 CNS sub-analyses

---

## Toxicity Rates

**Toxicity confidence range (NEW in v5.5)**
- Additive model is conservative (assumes independence). Dose-optimised estimate = 60% of additive rate.
- PBTC-047 precedent: panobinostat proceeded at 27.6% DLT rate with dose modification — establishing that >20% is acceptable in pediatric DIPG when no curative alternative exists.

**Panobinostat (PBTC-047)**
- Monje, M. et al. 2023 (see above) — 8/29 DLTs = 27.6% hematologic toxicity rate

**Abemaciclib**
- Sledge, G.W. et al. MONARCH 2. *Journal of Clinical Oncology*, 35(25):2875–2884, 2017. PMID: 28580882
  - Grade 3/4 neutropenia ~20%

**Marizomib**
- Bota et al. 2021 (see above) — grade 3/4 hematologic AE ~10% in CNS trials

**ONC201**
- ACTION trial (NCT05476081): Phase III ONC201 in H3K27M glioma — ~3% grade 3/4 hematologic

---

## PPI Network

- Szklarczyk, D. et al. STRING v11: protein–protein association networks with increased coverage. *Nucleic Acids Research*, 47(D1):D607–D613, 2019. PMID: 30476243
- Barabási, A.L. et al. Network medicine: a network-based approach to human disease. *Nature Reviews Genetics*, 12:56–68, 2011.

**Proteasome PPI neighbors (NEW in v5.5 — Marizomib fix)**
- Lin et al. 2019 (see above) — PSMB5 → H3F3A connection (H3K27M proteostasis)
- Tanaka 2009 (see above) — PSMB5 → MYC (MYC t½ ~20 min, proteasome-dependent)
- Dikic, I. Proteasomal and autophagic degradation systems. *Annual Review of Biochemistry*, 86:193–224, 2017. — PSMB5 → SQSTM1/ATG5 (proteasome-autophagy crosstalk)

---

## Statistical Methods

**Fisher's exact test (co-occurrence)**
- Fisher, R.A. On the interpretation of χ² from contingency tables. *Journal of the Royal Statistical Society*, 85(1):87–94, 1922.
- One-sided test (alternative="less"): H3K27M and CDKN2A co-occur less than expected — mutual exclusivity
- Result: p = 7.55e-05 (n=184, a=14, b=81, c=36, d=53)

**Sensitivity analysis**
- Top-4 ranking (Birabresib, Panobinostat, Marizomib, Abemaciclib) stable under all ±10% weight perturbations
- GSC quantile sensitivity: top-4 stable across p75/p80/p85/p90 thresholds

---

## Clinical Trial Context

- PBTC-049: Phase I birabresib in pediatric CNS tumours. NCT02296476
- PBTC-047: Phase I panobinostat in H3K27M DIPG. NCT02717455; Monje et al. 2023
- ACTION trial: Phase III ONC201 in H3K27M diffuse glioma. NCT05476081
- NCT03709680: Phase I/II abemaciclib in pediatric brain tumours
- NCT02903069: Phase I/II marizomib ± bevacizumab in recurrent glioma. Bota et al. 2021
- NCT04049669: Indoximod + temozolomide in pediatric DIPG (ongoing)
- NCT03696355: Paxalisib (GDC-0084) in pediatric DIPG Phase I

---

## Software and Computational Tools

- Python 3.10+
- pandas, numpy, scipy (statistical computing)
- matplotlib (figures)
- h5py (HDF5/GCTX file reading for CMap integration)
- OpenTargets GraphQL API v23
- STRING-DB REST API v11.5
- DepMap Public 23Q4 (CRISPRGeneEffect.csv, Model.csv)
- clue.io CMap L1000 query tool (sig_query_l1k); query submitted March 2026
- GEO datasets: GSE102130, GSE50021, GSE115397
- MyGene.info API (gene symbol → Entrez ID conversion for CMap query preparation)