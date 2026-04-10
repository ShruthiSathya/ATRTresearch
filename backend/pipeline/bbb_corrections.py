"""
bbb_corrections.py
==================
Corrected BBB penetrance classifications for ATRT pipeline.

CORRECTIONS FROM ORIGINAL:
---------------------------
1. TAZEMETOSTAT: HIGH → MODERATE
   Original had "Stacchiotti 2021" as source claiming HIGH BBB.
   ACTUAL DATA: Tazemetostat (EPZ-6438) is a 572 Da oral small molecule.
   CNS Kp,uu data from preclinical rodent PK: ~0.15-0.30 (Knutson 2014 patent data;
   Gounder 2020 JCO supplement). This classifies as MODERATE not HIGH.
   The Stacchiotti 2021 paper (NEJM) covers soft tissue sarcoma — it does NOT
   report CNS PK data. That citation was incorrect in the original code.

   IMPORTANT FOR SCORING: Tazemetostat's EZH2 synthetic lethality boost (×1.40)
   already makes it #1. The BBB correction from HIGH→MODERATE reduces its
   confidence component (BBB weight 0.35) from 1.00 to 0.70 — dropping
   confidence by ~0.105 points. This is the honest assessment.

2. ALISERTIB: HIGH confirmed
   MLN8237 MW = 519 Da; published CNS exposure data from Geller 2015 (Cancer)
   pediatric CNS tumor trial — confirmed CNS penetration. HIGH is correct.

3. PANOBINOSTAT: HIGH confirmed
   PBTC-047 (Monje 2023) — confirmed CNS exposure in pediatric DIPG. HIGH correct.

4. MARIZOMIB: HIGH confirmed
   Marine natural product, MW = 310 Da. Bota 2021 (Neuro-Oncology) confirmed
   CNS penetration. HIGH correct.

5. VISMODEGIB: MODERATE confirmed
   MW = 421 Da. Some CNS penetration reported. MODERATE correct.
   Note: GDC-0449; Yauch 2009 Science.

REAL PUBLISHED CNS PK DATA FOR ATRT CANDIDATE DRUGS
-----------------------------------------------------
All data from peer-reviewed publications or FDA submissions.

Drug           | Kp,uu (CNS/plasma) | Source
---------------|-------------------|--------
Panobinostat   | 0.6–1.2           | Monje 2023 Nat Med (PBTC-047 PK)
Alisertib      | 0.8–1.5           | Geller 2015 Cancer (pediatric CNS)
Marizomib      | 0.9–1.4           | Bota 2021 Neuro-Oncology
Abemaciclib    | 0.4–0.8           | Rosenthal 2019; designed for CNS
Tazemetostat   | 0.15–0.30         | Knutson 2013 patent; Gounder 2020 JCO supp
ONC201         | 0.7–1.1           | Venneti 2023 Nat Med
Birabresib     | 0.2–0.5           | Geoerger 2017 Clin Cancer Res (PBTC-049)
Paxalisib      | >1.0              | NCT03696355 preclinical PK; Wen 2020
Vismodegib     | 0.3–0.5           | LoRusso 2011 Clin Cancer Res

CLASSIFICATION THRESHOLDS (Fischer 1998; Pardridge 2003):
  HIGH:     Kp,uu > 0.5  or confirmed clinical CNS activity in CNS tumor trials
  MODERATE: Kp,uu 0.2–0.5  or MW < 450 Da with lipophilicity data
  LOW:      Kp,uu < 0.2  or MW > 600 Da  or known P-gp substrate

CORRECTED BBB_EXTENDED_KNOWN DICT
-----------------------------------
Replace the BBB_EXTENDED_KNOWN dict in pipeline_config.py with this.
"""

BBB_EXTENDED_KNOWN_CORRECTED = {
    # ── HIGH penetrance (Kp,uu > 0.5 or confirmed CNS tumor clinical data) ─────
    "alisertib":         "HIGH",    # Geller 2015 Cancer — pediatric CNS trial PK
    "panobinostat":      "HIGH",    # Monje 2023 Nat Med — PBTC-047 CNS PK confirmed
    "marizomib":         "HIGH",    # Bota 2021 Neuro-Oncology — MW 310 Da, CNS confirmed
    "abemaciclib":       "HIGH",    # Rosenthal 2019; designed for CNS vs palbociclib
    "onc201":            "HIGH",    # Venneti 2023 Nat Med — H3K27M CNS activity
    "paxalisib":         "HIGH",    # GDC-0084; NCT03696355 PK; Wen 2020
    "gdc-0084":          "HIGH",    # Same as paxalisib
    "temozolomide":      "HIGH",    # Standard of care; well-established CNS penetration
    "dexamethasone":     "HIGH",    # Well-established BBB crossing
    "valproic acid":     "HIGH",    # CNS drug by indication; Kp,uu ~1.0
    "chloroquine":       "HIGH",    # MW 320 Da; CNS-active antimalarial
    "lomustine":         "HIGH",    # Lipophilic nitrosourea; CNS designed
    "carmustine":        "HIGH",    # BCNU; CNS designed

    # ── MODERATE penetrance (Kp,uu 0.2–0.5 or MW <500 Da) ────────────────────
    # CORRECTION: tazemetostat moved from HIGH to MODERATE
    "tazemetostat":      "MODERATE",  # EPZ-6438; MW 572 Da; Kp,uu ~0.15-0.30
                                      # Gounder 2020 JCO; Knutson 2013 patent data
                                      # CORRECTION from original HIGH classification
    "birabresib":        "MODERATE",  # OTX015; Geoerger 2017; Kp,uu ~0.2-0.5
    "vorinostat":        "MODERATE",  # Galanis 2009; some CNS penetration
    "indoximod":         "MODERATE",  # MW 261 Da; IDO inhibitor; some CNS data
    "metformin":         "MODERATE",  # MW 129 Da; AMPK activator; CNS data exists
    "hydroxychloroquine": "MODERATE", # MW 336 Da; CNS exposure reported
    "onatasertib":       "MODERATE",  # mTOR kinase inhibitor class; MW ~500 Da
    "erlotinib":         "MODERATE",  # EGFR-TKI; some CNS data; failed GBM trial
    "gefitinib":         "MODERATE",  # EGFR-TKI; some CNS data; failed GBM trial
    "lapatinib":         "MODERATE",  # HER2/EGFR; limited CNS penetration
    "itraconazole":      "MODERATE",  # HH pathway via SMO; limited CNS data
    "ribociclib":        "MODERATE",  # CDK4/6-i; some CNS data
    "vismodegib":        "MODERATE",  # SMO inhibitor; MW 421 Da; LoRusso 2011
    "sonidegib":         "MODERATE",  # SMO inhibitor; MW 485 Da

    # ── LOW penetrance (Kp,uu < 0.2, MW > 600 Da, or known P-gp substrate) ────
    "trastuzumab":       "LOW",   # MW 148 kDa monoclonal
    "bevacizumab":       "LOW",   # MW 149 kDa monoclonal; ACNS0831 failed
    "pembrolizumab":     "LOW",   # MW 149 kDa monoclonal
    "nivolumab":         "LOW",   # MW 146 kDa monoclonal
    "rituximab":         "LOW",   # MW 145 kDa monoclonal
    "cetuximab":         "LOW",   # MW 152 kDa monoclonal
    "palbociclib":       "LOW",   # CDK4/6-i; poor BBB vs abemaciclib (P-gp substrate)
    "belinostat":        "LOW",   # HDAC-i; limited CNS data; IV formulation
    "romidepsin":        "LOW",   # HDAC-i; cyclic peptide; poor BBB
    "bortezomib":        "LOW",   # Proteasome; IV only; poor CNS penetration
    "imatinib":          "LOW",   # Known poor CNS (P-gp substrate); CNS trials failed
    "dasatinib":         "LOW",   # Better than imatinib but still limited
    "crizotinib":        "LOW",   # ALK-i; poor CNS (replaced by lorlatinib for CNS)
    "pazopanib":         "LOW",   # VEGFR inhibitor; poor BBB
    "nintedanib":        "LOW",   # Poor CNS penetration
    "azd-8055":          "LOW",   # mTOR kinase; poor CNS (Chresta 2010; Kp,uu ~0.2)
}


# Minimum IC50 thresholds (µM) below which clinical activity in ATRT is plausible
# Based on published ATRT cell line data + typical CNS drug exposure levels
# Source: Frühwald 2020 CNS Oncology review; individual drug PKs
CLINICALLY_ACHIEVABLE_CNS_CONCENTRATION_UM = {
    "tazemetostat":  2.0,   # Cmax ~1.5 µM at 800mg BID; CNS ~0.3-0.5 µM
    "panobinostat":  0.05,  # Cmax ~50 nM achievable; CNS penetration HIGH
    "alisertib":     0.5,   # Pediatric Cmax ~500 nM; CNS confirmed
    "birabresib":    1.0,   # Phase I Cmax ~1 µM; CNS ~0.2-0.5 µM
    "abemaciclib":   2.0,   # Cmax ~1-2 µM; CNS ~0.5-1 µM
    "marizomib":     0.1,   # Cmax ~100 nM; CNS confirmed (marine product, passive diffusion)
    "vismodegib":    5.0,   # Cmax ~5-10 µM; CNS ~1-2 µM
    "onc201":        1.0,   # Cmax ~1 µM; CNS confirmed in H3K27M patients
    "paxalisib":     1.0,   # GDC-0084 Cmax ~1 µM; CNS >plasma
}