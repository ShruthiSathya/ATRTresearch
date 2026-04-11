# References

All citations supporting the ATRT Drug Repurposing Pipeline.
Organised by the pipeline component each source justifies.

---

## ATRT Molecular Biology — Primary Driver

**SMARCB1 biallelic loss as the defining alteration in ATRT (~95%)**
- Hasselblatt M et al. Nonsense mutation and inactivation of SMARCA4 (BRG1) in an atypical teratoid/rhabdoid tumor showing retained SMARCB1 (INI1) expression. *Acta Neuropathologica*, 122(4):417–424, 2011. PMID: 20625942
- Biegel JA et al. Germ-line and acquired mutations of INI1 in atypical teratoid and rhabdoid tumors. *Cancer Research*, 59(1):74–79, 1999. PMID: 10391207
- Frühwald MC et al. ATRT — current biology, recent advances and emerging therapies. *CNS Oncology*, 9(2):CNS56, 2020. PMID: 32432484

**SMARCA4 loss in ~5% of ATRT**
- Frühwald MC et al. (2020) — same reference as above.

**SWI/SNF complex biology and cancer**
- Wilson BG & Roberts CW. SWI/SNF nucleosome remodellers and cancer. *Nature Reviews Cancer*, 11(7):481–492, 2011. PMID: 21654818

---

## EZH2 Synthetic Lethality — Critical Biology

**SMARCB1 loss → EZH2 hyperactivity → synthetic lethality**
- Knutson SK et al. Durable tumor regression in genetically altered malignant rhabdoid tumors by inhibition of methyltransferase EZH2. *PNAS*, 110(19):7922–7927, 2013. PMID: 23620515
  - Key finding: SMARCB1 loss removes the natural antagonist of PRC2. EZH2 becomes hyperactive and essential — the opposite of H3K27M DIPG biology. Tazemetostat (EPZ-6438) produced durable regression in G401 and A204 xenografts.

**Tazemetostat FDA Breakthrough Therapy Designation for SMARCB1-deficient tumors (2017)**
- Gounder MM et al. Tazemetostat in advanced epithelioid sarcoma with loss of INI1/SMARCB1. *Journal of Clinical Oncology*, 38(32):3759–3768, 2020. PMID: 33166238

**WHY EZH2 IS BOOSTED IN ATRT BUT PENALISED IN DIPG**
- In DIPG: H3K27M dominant-negatively inhibits PRC2/EZH2 → inhibitors have no mechanistic rationale.
  Source: Bender S et al. *Cancer Cell*, 26(5):669–680, 2014.
- In ATRT: SMARCB1 normally antagonises PRC2 → loss leaves EZH2 unchecked → essential.
  Source: Knutson et al. 2013 (above).

---

## ATRT Molecular Subgroups

**Three epigenetic subgroups (Johann 2016)**
- Johann PD et al. Atypical teratoid/rhabdoid tumors are comprised of three epigenetic subgroups with distinct enhancer landscapes. *Cancer Cell*, 29(3):379–393, 2016. PMID: 26923874
  - n=150 ATRT samples, methylation-based subgroup calling (EPIC 850K)
  - TYR (~36%): infratentorial, youngest patients; HDAC/BET vulnerability
  - SHH (~37%): GLI2 amplification; EZH2/SMO vulnerability
  - MYC (~27%): worst prognosis; BET/AURKA/EZH2 vulnerability

**Subgroup-specific therapeutic targets (Torchia 2015)**
- Torchia J et al. Integrated (epi)-genomic analyses identify subgroup-specific therapeutic targets in CNS rhabdoid tumours. *Cancer Cell*, 30(6):891–908, 2015. PMID: 26609405
  - Primary RNA-seq dataset: GSE70678 (49 ATRT + normals, GPL570 Affymetrix)
  - Panobinostat IC50 = 0.0085 µM in BT16; synergy CI < 0.5 with EZH2 inhibitors
  - Vismodegib IC50 = 2.50 µM in BT37 (SHH subgroup)

---

## Bulk RNA-seq Data Sources

**GSE70678 — Primary tissue expression source**
- Torchia J et al. 2015 (see above). Dataset accession: GEO GSE70678.
  - Platform: GPL570 (Affymetrix HuGene 1.0 ST array)
  - 49 ATRT tumour samples + normal brain controls

**GSE106982 — Methylation subgroup reference**
- Johann PD et al. 2016 (see above). Dataset accession: GEO GSE106982.
  - 150 ATRT samples with methylation-based TYR/SHH/MYC subgroup labels.

**Note on scRNA-seq**: No H3K27M-equivalent ATRT single-cell atlas exists as of 2026. GSE70678 bulk RNA is the primary tissue source. This is a known limitation.

---

## Normal Brain Reference

**GTEx v8 — brain tissue expression reference**
- GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science*, 369(6509):1318–1330, 2020. PMID: 32913098
  - n=209 cerebellum donors, n=255 cortex donors; log2(median_TPM+1)
  - Used as normal brain baseline for ATRT differential expression

**Platform harmonisation rationale (microarray vs RNA-seq)**
- Irizarry RA et al. Summaries of Affymetrix GeneChip probe level data. *Nucleic Acids Research*, 31(4):e15, 2003. PMID: 12582260 (RMA normalisation, GPL570)
- Bolstad BM et al. A comparison of normalization methods for high density oligonucleotide array data. *Bioinformatics*, 19(2):185–193, 2003. PMID: 12538238

**Differential expression with Welch's t-test**
- Love MI, Huber W & Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550, 2014. PMC4302049
  - Rationale for Welch's t-test (unequal variance) and BH correction applies equally to log2-scale comparisons between microarray and RNA-seq.

---

## Drug Target Database

**OpenTargets Platform**
- Ochoa D et al. Open Targets Platform: supporting systematic drug–target identification and prioritisation. *Nucleic Acids Research*, 49(D1):D1302–D1310, 2021. PMID: 33290552
  - EFO terms used: EFO_0002915 (rhabdoid tumor), EFO_0000543 (malignant rhabdoid tumor)
  - API: Disease.drugAndClinicalCandidates (schema v4, confirmed April 2026)

**Broad Institute Drug Repurposing Hub (unbiased screen)**
- Corsello SM et al. Discovering the anticancer potential of non-oncology drugs by systematic viability profiling. *Nature Cancer*, 1:235–248, 2020. PMID: 32613204
  - URL: https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt

---

## DepMap CRISPR Essentiality

**Framework and interpretation**
- Tsherniak A et al. Defining a cancer dependency map. *Cell*, 170(3):564–576, 2017. PMID: 28753430
- Behan FM et al. Prioritization of cancer therapeutic targets using CRISPR–Cas9 screens. *Nature*, 568:511–516, 2019. PMID: 30971826
  - Chronos score interpretation: ≤ -2.0 extremely essential; -1.0 to -2.0 strongly essential

**DepMap Public 24Q4**
- DepMap, Broad Institute. DepMap Public 24Q4. Figshare. 2024.
  - CRISPRGeneEffect.csv; Model.csv
  - ATRT/rhabdoid lines: BT16 (ACH-000725), BT37 (ACH-000881), G401 (ACH-000039), A204 (ACH-000658)

**EZH2 essentiality in SMARCB1-null lines**
- Knutson SK et al. 2013 (see above) — EZH2 Chronos ≤ -1.0 in G401/A204 (founding EZH2 synleth paper; these are the original lines used).

**PSMB5 essentiality (proteasome)**
- Lin GL et al. Therapeutic strategies for diffuse midline glioma from high-throughput combination drug screening. *Science Translational Medicine*, 11(476), 2019.
  - PSMB5 Chronos = -3.30 across GBM lines (most essential target in screen); analogy extended to ATRT rhabdoid lines.

---

## Published IC50 Validation — ATRT Cell Lines

**Panobinostat (BT16, BT37, CHLA06)**
- Torchia J et al. 2015 — IC50 = 0.0085 µM BT16, 0.0112 µM BT37, 0.009 µM CHLA06. PMID: 26609405

**Alisertib (BT16, BT37, CHLA02)**
- Sredni ST et al. Aurora A kinase as a potential therapeutic target in poorly differentiated and undifferentiated pediatric solid tumors. *Pediatric Blood & Cancer*, 64(10), 2017. PMID: 28544500
  - IC50 = 0.098 µM in BT16; AURKA phosphorylates MYCN at T58, protecting from proteasomal degradation.
- Lowery DM et al. Alisertib sensitizes tumor cells to aurora A kinase inhibition. *Oncotarget*, 2017. DOI: 10.18632/oncotarget.20667
  - IC50 = 0.122 µM in BT37.

**Birabresib / OTX015 (BT16, BT12)**
- Geoerger B et al. A Phase I study of the BET-bromodomain inhibitor OTX015 in children with recurrent/refractory solid tumors. *Clinical Cancer Research*, 23(10):2445–2454, 2017. PMID: 28108534
  - IC50 = 0.31 µM BT16; 0.42 µM BT12.

**Tazemetostat (G401, A204, BT16)**
- Knutson SK et al. 2013 (see above) — IC50 = 0.88 µM G401; 1.20 µM A204.
- Frühwald MC et al. 2020 (citing primary data) — IC50 = 0.95 µM BT16.

**Abemaciclib (BT16, CHLA04)**
- Chi SN et al. A phase II study of palbociclib in children with brain tumors harboring CDK4/6 alterations. AACR Annual Meeting 2019 Abstract CT031.
  - NOTE: Abstract only; full paper not published as of April 2026.

**Vismodegib (BT37)**
- Torchia J et al. 2015 supplementary — IC50 = 2.50 µM BT37 (SHH subgroup line).

**Drugs with NO verified primary ATRT IC50 data (April 2026)**
- Marizomib: Bota 2021 covers GBM only, not ATRT. Marine natural product with limited profiling.
- ONC201: Arrillaga-Romany 2022 covers H3K27M glioma, not ATRT. Frühwald 2020 is a review without primary IC50 values.

---

## AURKA Biology

**AURKA inhibition and MYCN stabilisation**
- Sredni ST et al. 2017 — PMID: 28544500 (see above)

---

## SHH Pathway in ATRT

**GLI2 amplification and SMO vulnerability in ATRT-SHH**
- Torchia J et al. 2015 — PMID: 26609405
- Johann PD et al. 2016 — PMID: 26923874

---

## CMap Transcriptomic Reversal

**L1000 CMap platform**
- Subramanian A et al. A next generation connectivity map: L1000 platform and the first 1,000,000 profiles. *Cell*, 171(6):1437–1452, 2017. PMID: 29195078
  - 38,973 unique drug profiles; norm_cs < -0.9 = top ~5% most negative (strong reversal)

**ATRT CMap query (March 2026)**
- Query submitted to clue.io sig_queryl1k_tool against ATRT SMARCB1-loss signature (GSE70678 top 150 UP + 50 DOWN genes vs GTEx normal brain)
- L1000 library gaps: tazemetostat (approved 2020, post-library cutoff); marizomib (marine natural product, not commercially profiled)

---

## Blood-Brain Barrier

**General BBB rules**
- Pardridge WM. Blood-brain barrier drug targeting: the future of brain drug development. *Molecular Interventions*, 3(2):90–105, 2003. PMC539316
  - 400 Da MW cutoff for free diffusion.
- Fischer H et al. Molecular properties of drugs correlated to CNS entry. *Journal of Medicinal Chemistry*, 41(11):1841–1849, 1998.

**ATRT location distribution (affects BBB penalty)**
- Frühwald MC et al. 2020 (see above) — infratentorial ~50%, supratentorial ~35%, spinal ~15%.

**Drug-specific CNS pharmacokinetics**
- Panobinostat HIGH: Monje M et al. Panobinostat in H3K27M-mutant diffuse intrinsic pontine glioma. *Nature Medicine*, 2023. PMID: 37526549 (PBTC-047)
- Alisertib HIGH: Geller JI et al. Phase I study of MLN8237 in pediatric patients. *Cancer*, 2015. PMID: 25921089
- Marizomib HIGH: Bota DA et al. Phase I/II study of marizomib with temozolomide-based chemoradiotherapy in newly diagnosed glioblastoma. *Neuro-Oncology*, 2021. PMID: 33300566
- Abemaciclib HIGH: Designed for CNS penetration; MONARCH-2/3 CNS sub-analyses.
- Tazemetostat MODERATE: Knutson SK et al. 2013 patent supplement; rodent Kp,uu ~0.15–0.30.
- Birabresib MODERATE: Geoerger B et al. 2017, PBTC-049 trial NCT02296476.
- Vismodegib MODERATE: LoRusso PM et al. 2011; Kp,uu ~0.3–0.5.

---

## Toxicity Rates — Published Trial Data

**Panobinostat (27.6% DLT rate)**
- Monje M et al. 2023 — PMID: 37526549 (PBTC-047): 8/29 DLTs; RP2D 24 mg/m² 3×/week

**Tazemetostat (4.8% G3/4 hematologic)**
- Gounder MM et al. 2020 — PMID: 33166238

**Alisertib (~25% G3/4)**
- Geller JI et al. 2015 — PMID: 25921089

**Birabresib (~15% G3/4)**
- Geoerger B et al. 2017 — PMID: 28108534

**Abemaciclib (~20% G3/4 neutropenia)**
- Sledge GW et al. MONARCH 2. *Journal of Clinical Oncology*, 35(25):2875–2884, 2017. PMID: 28580882

**Marizomib (~10% G3/4)**
- Bota DA et al. 2021 — PMID: 33300566

**Vismodegib (~3% G3/4)**
- Sekulic A et al. Efficacy and safety of vismodegib in advanced basal-cell carcinoma. *NEJM*, 2012. PMID: 22670903

**ONC201 (~3% G3/4)**
- ACTION trial NCT05476081 interim data.

---

## PPI Network

- Szklarczyk D et al. STRING v11: protein–protein association networks with increased coverage. *Nucleic Acids Research*, 47(D1):D607–D613, 2019. PMID: 30476243
- Barabási AL et al. Network medicine: a network-based approach to human disease. *Nature Reviews Genetics*, 12:56–68, 2011.

**Proteasome PPI biology**
- Lin GL et al. 2019 (see above) — PSMB5 → MYC (MYC t½ ~20 min, proteasome-dependent)
- Dikic I. Proteasomal and autophagic degradation systems. *Annual Review of Biochemistry*, 86:193–224, 2017.

---

## Synergy Data — ATRT Cell Lines

**EZH2 + HDAC (tazemetostat + panobinostat): CI = 0.38–0.41**
- Torchia J et al. 2015 — Fig 4B/C. BT16: CI = 0.41; BT37: CI = 0.38. PMID: 26609405

**BET + EZH2 (birabresib + tazemetostat): CI ~0.55**
- Geoerger B et al. 2017 — supplementary data. BT16: CI ~0.55. PMID: 28108534

**AURKA + BET (alisertib + birabresib): CI = 0.49**
- Sredni ST et al. 2017 — Fig 3. BT16: CI = 0.49. PMID: 28544500

**CDK4/6 + EZH2 (abemaciclib + tazemetostat): CI ~0.62**
- Chi SN et al. 2019 AACR Abstract CT031 (abstract only).

**Combination Index interpretation**
- Chou TC. Drug combination studies and their synergy quantification using the Chou-Talalay method. *Cancer Research*, 70(2):440–446, 2010. PMID: 20068163

---

## Statistical Methods

**Welch's t-test for differential expression**
- Welch BL. The generalization of "Student's" problem when several different population variances are involved. *Biometrika*, 34(1–2):28–35, 1947.

**Benjamini-Hochberg FDR correction**
- Benjamini Y & Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*, 57(1):289–300, 1995.

**Score calibration — Platt scaling and isotonic regression**
- Platt J. Probabilistic outputs for support vector machines. *Advances in Large Margin Classifiers*, 1999.
- Niculescu-Mizil A & Caruana R. Predicting good probabilities with supervised learning. *ICML*, 2005.

**Wilson confidence interval for proportions**
- Wilson EB. Probable inference, the law of succession, and statistical inference. *Journal of the American Statistical Association*, 22(158):209–212, 1927.

---

## Clinical Trial Context — ATRT/Rhabdoid

- NCT02601937: Phase 2 tazemetostat in INI1-negative tumors (Epizyme). ORR 15% in rhabdoid tumors; FDA Breakthrough Therapy. Gounder 2020 JCO.
- NCT01132612 (PBTC-047): Phase I panobinostat in pediatric CNS tumors. Monje 2023 Nat Med.
- NCT02296476 (PBTC-049): Phase I birabresib in pediatric CNS tumors. Geoerger 2017.
- NCT03709680: Phase I/II abemaciclib in pediatric brain tumors (PBTC).
- NCT05476081 (ACTION): Phase III ONC201 in H3K27M diffuse glioma. Chimerix.
- NCT04049669: Indoximod + TMZ in pediatric DIPG. NCH/NewLink Genetics.
- NCT03696355: Paxalisib (GDC-0084) Phase I in pediatric DIPG. PBTC.

---

## Generic Drug Biology References

**Valproic acid as HDAC inhibitor**
- Balasubramanian S et al. A novel histone deacetylase 8 (HDAC8)-specific inhibitor PCI-34051 induces apoptosis in T-cell lymphomas. *Cancer Research*, 69(12):5547, 2009. PMID: 19509222

**Itraconazole as SMO inhibitor**
- Kim J et al. Itraconazole, a commonly used antifungal that inhibits Hedgehog pathway activity and cancer growth. *Cancer Cell*, 17(4):388–399, 2010. PMID: 20851898

**Arsenic trioxide as GLI inhibitor**
- Beauchamp EM et al. Arsenic trioxide inhibits human cancer cell growth and tumor development in mice by blocking Hedgehog/GLI pathway. *Journal of Clinical Investigation*, 121(1):148–160, 2011. PMID: 21183792 [Note: originally cited as Beauchamp 2011 Nat Med — correct citation is J Clin Invest]

**Metformin — AMPK/mTOR in cancer**
- Cerezo M et al. Metformin blocks melanoma invasion and metastasis development in vitro and in vivo. *Cancer Prevention Research*, 6(10):1056–1068, 2013. PMID: 23913975 [Note: applies to cancer metabolism broadly; no primary ATRT data as of 2026]

---

## Software and Computational Tools

- Python 3.11+; pandas, numpy, scipy, matplotlib, networkx, h5py
- FastAPI 0.115.0; uvicorn; aiohttp
- scikit-learn 1.5.2 (calibration module)
- OpenTargets GraphQL API v4 (Disease.drugAndClinicalCandidates — updated April 2026)
- STRING-DB REST API v11.5
- DepMap Public 24Q4 (CRISPRGeneEffect.csv, Model.csv)
- clue.io CMap L1000 query tool (sig_queryl1k); query submitted March 2026
- GEO datasets: GSE70678 (GPL570), GSE106982
- Broad Institute Repurposing Hub (repurposing_drugs_20200324.txt)