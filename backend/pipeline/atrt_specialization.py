"""
atrt_specialization.py
======================
ATRT-specific scoring, gene sets, and subgroup stratification.

ATRT BIOLOGY OVERVIEW
----------------------
Atypical Teratoid/Rhabdoid Tumor (ATRT) is a highly aggressive pediatric CNS tumour
driven almost exclusively by biallelic loss of SMARCB1 (INI1, ~95% of cases) or
SMARCA4 (~5%). Both encode subunits of the SWI/SNF (BAF) chromatin remodelling complex.

CRITICAL BIOLOGICAL INVERSION VS DIPG
---------------------------------------
In DIPG: H3K27M inhibits PRC2/EZH2 → EZH2 inhibitors are PENALISED (non-rational).
In ATRT: SMARCB1 loss REMOVES the antagonist of PRC2 → EZH2 becomes hyperactive
         and ESSENTIAL. EZH2 inhibitors (tazemetostat) are STRONGLY RATIONAL.
         This is a textbook synthetic lethality relationship.

Source: Knutson et al. 2013 PNAS; Wilson & Roberts 2011 Nat Rev Cancer.
Tazemetostat received FDA Breakthrough Therapy Designation for SMARCB1-deficient tumors.

THREE MOLECULAR SUBGROUPS (Torchia et al. 2015; Johann et al. 2016)
----------------------------------------------------------------------
ATRT-TYR  : Tyrosinase expression. Youngest patients (median age ~1yr). Infratentorial.
             Somewhat better prognosis. Driven by neural crest / melanocytic programs.
             Key vulnerability: HDAC + BET.

ATRT-SHH  : Sonic Hedgehog pathway activation. Mixed ages. Supratentorial predominance.
             GLI2 amplification common. Key vulnerability: SMO/GLI inhibitors + EZH2.

ATRT-MYC  : MYC amplification/upregulation. Worst prognosis. Supratentorial predominant.
             Oldest patients in ATRT. Key vulnerability: BET bromodomain (BRD4→MYC),
             Aurora kinase A (stabilises MYCN), CDK4/6, EZH2.

BBB CONSIDERATIONS
-------------------
Unlike DIPG (always brainstem), ATRT location varies:
  - Infratentorial/posterior fossa: ~50% (tighter BBB, similar to DIPG)
  - Supratentorial: ~35% (standard BBB)
  - Spinal: ~15% (no BBB relevance)
Without location data, use MODERATE BBB assumption (not DIPG's severe brainstem penalty).

DATA SOURCES REQUIRED
---------------------
1. RNA-seq cohort:
   GSE70678 — Torchia et al. 2015 Cell. 49 ATRT + normal samples. Subgroup labels included.
   Access: GEO (open). Platform: Affymetrix HuGene 1.0 ST (GPL6244).
   Key columns: subgroup annotation (TYR/SHH/MYC), SMARCB1 status.

2. Methylation-based subgroup:
   GSE106982 — Johann et al. 2016 Cancer Cell. 150 ATRT. Gold-standard subgroup calls.
   Access: GEO (open). Used for subgroup prevalence priors.

3. Additional RNA:
   CBTN ATRT cohort on Cavatica (OpenPedCan). Controlled access via DCC.
   ~30–40 ATRT samples. Same portal as PNOC/PBTA.
   URL: https://portal.kidsfirstdrc.org

4. scRNA-seq (limited):
   No H3K27M-equivalent scRNA atlas exists for ATRT as of 2026.
   Use GSE70678 bulk RNA as primary tissue expression source.
   Curated cell-line expression from published ATRT cell line studies as fallback.

5. DepMap cell lines — ATRT/rhabdoid:
   OncotreeSubtype filter: "MRT" (malignant rhabdoid tumor), "ATRT"
   Known lines: BT16, BT37, G401, A204, KP-MRT-NS, MON, BT12, CHLA02, CHLA04, CHLA06
   Many are SMARCB1-null and EZH2-dependent — critical for validating EZH2 scoring.

6. Published IC50 data: see published_ic50_atrt_validation.py

REFERENCES
-----------
Torchia J et al. (2015). Integrated (epi)-genomic analyses identify subgroup-specific
  therapeutic targets in CNS rhabdoid tumours. Cancer Cell, 30(6):891–908. PMID: 264...
  [Cell, 30(6):891-908, 2015. PMID 26609405]

Johann PD et al. (2016). Atypical teratoid/rhabdoid tumors are comprised of three
  epigenetic subgroups with distinct enhancer landscapes. Cancer Cell, 29(3):379–393.
  PMID 26923874.

Knutson SK et al. (2013). Durable tumor regression in genetically altered malignant
  rhabdoid tumors by inhibition of methyltransferase EZH2. PNAS 110(19):7922–7927.
  PMID 23620515.

Wilson BG & Roberts CW (2011). SWI/SNF nucleosome remodellers and cancer.
  Nature Reviews Cancer, 11(7):481–492. PMID 21654818.

Frühwald MC et al. (2020). ATRT—current biology, recent advances and emerging
  therapies. CNS Oncology, 9(2):CNS56. PMID 32432484.

Sredni ST et al. (2017). Aurora A kinase as a potential therapeutic target in
  poorly differentiated and undifferentiated pediatric solid tumors. Pediatric Blood
  & Cancer, 64(10). PMID 28544500.

Daigle SR et al. (2011). Selective killing of mixed lineage leukemia cells by a
  potent small-molecule DOT1L inhibitor. Cancer Cell, 20(1):53–65. PMID 21741596.
  [Note: Alisertib (MLN8237) AURKA data from Lowery et al. 2017 Oncotarget]

Chi SN et al. (2019). A phase II study of palbociclib (PD0332991) in children with
  brain tumors harboring CDK4/6 alterations. AACR 2019 Abstract.
"""

import logging
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# ATRT CORE GENE SETS
# ─────────────────────────────────────────────────────────────────────────────

ATRT_CORE_GENES: List[str] = [
    # SWI/SNF complex — primary driver
    "SMARCB1", "SMARCA4", "SMARCA2", "SMARCC1", "SMARCC2",
    "SMARCD1", "SMARCD2", "SMARCD3", "SMARCE1", "ARID1A",
    "ARID1B", "ARID2", "PBRM1",
    # PRC2 complex — synthetic lethality target
    "EZH2", "EED", "SUZ12", "RBBP4", "RBBP7",
    # BET bromodomain — H3K27ac reader, super-enhancer dependency
    "BRD4", "BRD2", "BRD3",
    # MYC — upregulated in ATRT-MYC subgroup, proteasome-sensitive
    "MYC", "MYCN", "MAX",
    # Aurora kinase — MYCN stabiliser, mitotic checkpoint
    "AURKA", "AURKB",
    # CDK4/6 — cell cycle
    "CDK4", "CDK6", "CCND1", "CCND2", "RB1", "CDKN2A",
    # HDAC
    "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC6",
    # mTOR / PI3K
    "MTOR", "PIK3CA", "AKT1", "PTEN",
    # SHH subgroup
    "GLI2", "GLI1", "SMO", "PTCH1", "SHH",
    # Apoptosis
    "BCL2", "BCL2L1", "MCL1", "TP53",
    # Proteasome — MYC degradation
    "PSMB5", "PSMB2", "PSMB1",
    # Stemness
    "SOX2", "SOX9", "LIN28A", "SALL4",
    # TYR subgroup markers
    "TYR", "DCT", "MITF",
    # DNA damage
    "PARP1", "ATM", "ATR",
    # Immune
    "CD274", "PDCD1",
]

# Subgroup-specific gene sets (Torchia 2015; Johann 2016)
ATRT_TYR_GENES: Set[str] = {
    "TYR", "DCT", "MITF", "SOX10", "PMEL",
    "HDAC1", "HDAC2", "BRD4",
    "SMARCB1",
}

ATRT_SHH_GENES: Set[str] = {
    "GLI2", "GLI1", "SMO", "PTCH1", "SHH", "HHIP",
    "EZH2", "BRD4", "CDK4",
    "SMARCB1",
}

ATRT_MYC_GENES: Set[str] = {
    "MYC", "MYCN", "AURKA", "CDK4", "CDK6",
    "BRD4", "BRD2", "EZH2",
    "PSMB5",  # MYC degraded by proteasome, t½ ~20 min
    "SMARCB1",
}

# Subgroup prevalence priors (Johann 2016, n=150)
ATRT_SUBGROUP_PREVALENCE = {
    "TYR": 0.36,
    "SHH": 0.37,
    "MYC": 0.27,
}

# Gene score weights for ATRT — higher = more essential/targetable
ATRT_GENE_SCORE_WEIGHTS: Dict[str, float] = {
    # Primary driver — all ATRT
    "SMARCB1": 1.00,  # Loss defines the disease
    "SMARCA4": 0.95,
    "EZH2":    1.00,  # SYNTHETIC LETHALITY with SMARCB1 loss — highest priority
    "EED":     0.90,
    "SUZ12":   0.88,
    # BET bromodomain
    "BRD4":    0.90,
    "BRD2":    0.82,
    "BRD3":    0.79,
    # HDAC
    "HDAC1":   0.85,
    "HDAC2":   0.82,
    "HDAC3":   0.80,
    # MYC/MYCN — especially ATRT-MYC subgroup
    "MYC":     0.85,
    "MYCN":    0.80,
    # Aurora kinase — MYCN stabiliser
    "AURKA":   0.82,
    "AURKB":   0.75,
    # CDK4/6
    "CDK4":    0.80,
    "CDK6":    0.75,
    # Proteasome
    "PSMB5":   0.72,
    "PSMB2":   0.68,
    "PSMB1":   0.65,
    # mTOR
    "MTOR":    0.68,
    "PIK3CA":  0.65,
    # SHH
    "GLI2":    0.70,
    "SMO":     0.65,
    # Apoptosis
    "BCL2":    0.60,
    "BCL2L1":  0.62,
    "MCL1":    0.65,
    # Stemness
    "LIN28A":  0.70,
    "SOX2":    0.72,
    "SALL4":   0.68,
    # PARP
    "PARP1":   0.62,
    # Tumour suppressors (lost)
    "CDKN2A":  0.20,
    "PTEN":    0.25,
    "RB1":     0.30,
    "TP53":    0.50,
}

# Pathway weights for ATRT composite scoring
ATRT_PATHWAY_WEIGHTS: Dict[str, float] = {
    "SWI/SNF complex":                   1.00,
    "SMARCB1 chromatin remodelling":     1.00,
    "EZH2 histone methyltransferase":    1.00,  # NOTE: BOOST not penalise
    "PRC2 complex":                      1.00,
    "H3K27 methylation":                 0.95,
    "Epigenetic regulation":             0.90,
    "BET bromodomain":                   0.90,
    "HDAC deacetylase activity":         0.88,
    "Histone deacetylation":             0.85,
    "MYC signaling":                     0.88,
    "MYCN signaling":                    0.85,
    "Aurora kinase signaling":           0.85,
    "CDK4/6 signaling":                  0.82,
    "Cell cycle regulation":             0.78,
    "SHH signaling":                     0.80,  # ATRT-SHH subgroup
    "Sonic Hedgehog pathway":            0.80,
    "PI3K-Akt signaling":                0.72,
    "mTOR signaling":                    0.70,
    "Proteasome inhibition":             0.75,
    "Apoptosis":                         0.70,
    "DNA damage response":               0.65,
    "PARP signaling":                    0.65,
    "T-cell checkpoint signaling":       0.55,
}


# ─────────────────────────────────────────────────────────────────────────────
# EZH2 BOOST (OPPOSITE OF DIPG PENALTY)
# ─────────────────────────────────────────────────────────────────────────────

# ATRT-specific EZH2 inhibitor config — a BOOST, not a penalty.
# SMARCB1 normally opposes PRC2; its loss makes EZH2 hyperactive and essential.
# Source: Knutson et al. 2013 PNAS. Tazemetostat FDA Breakthrough Therapy.
EZH2_INHIBITOR_ATRT = {
    "known_inhibitors": {
        "tazemetostat", "gsk126", "epz-6438", "epz6438",
        "unc1999", "unc-1999", "cpi-1205", "cpi1205",
        "ds-3201", "ds3201", "valemetostat",
    },
    "mechanism_keywords": [
        "ezh2 inhibitor", "prc2 inhibitor", "ezh2 inhibition",
        "histone methyltransferase inhibitor",
    ],
    # BOOST composite score — opposite of DIPG
    "composite_boost": 1.40,
    "rationale": (
        "EZH2 inhibitor BOOSTED in ATRT: SMARCB1/SMARCA4 loss removes the "
        "natural antagonist of PRC2, making EZH2 hyperactive and essential. "
        "Synthetic lethality relationship confirmed in Knutson et al. 2013 PNAS. "
        "Tazemetostat received FDA Breakthrough Therapy Designation for "
        "SMARCB1-deficient tumors. OPPOSITE of DIPG biology."
    ),
}

# AURKA inhibitor config — relevant for ATRT-MYC subgroup (MYCN stabilisation)
AURKA_INHIBITOR_ATRT = {
    "known_inhibitors": {
        "alisertib", "mln8237", "mln-8237",
        "barasertib", "at9283", "vx-680",
    },
    "mechanism_keywords": [
        "aurora kinase a inhibitor", "aurka inhibitor",
        "aurora a inhibitor",
    ],
    "composite_boost_myc_subgroup": 1.30,
    "composite_boost_other": 1.15,
    "rationale": (
        "AURKA inhibition destabilises MYCN in ATRT-MYC subgroup. "
        "AURKA phosphorylates MYCN at T58, protecting it from proteasomal degradation. "
        "Source: Sredni et al. 2017 Pediatric Blood & Cancer; "
        "Lowery et al. 2017 Oncotarget — alisertib IC50 ~100nM in BT16/BT37."
    ),
}

# SMO/GLI inhibitor config — relevant for ATRT-SHH subgroup
SMO_INHIBITOR_ATRT = {
    "known_inhibitors": {
        "vismodegib", "sonidegib", "glasdegib",
        "taladegib", "saridegib",
    },
    "mechanism_keywords": [
        "smoothened inhibitor", "smo inhibitor",
        "hedgehog pathway inhibitor", "gli inhibitor",
    ],
    "composite_boost_shh_subgroup": 1.25,
    "composite_boost_other": 1.00,  # Neutral in non-SHH
    "rationale": (
        "SMO/GLI inhibition is rationale in ATRT-SHH subgroup (~37% of ATRT). "
        "GLI2 amplification drives SHH pathway activation. "
        "Source: Torchia et al. 2015 Cancer Cell; Johann et al. 2016 Cancer Cell."
    ),
}


# ─────────────────────────────────────────────────────────────────────────────
# RESISTANCE BYPASS MAP (ATRT-specific)
# Based on SMARCB1-loss biology and published resistance mechanisms
# ─────────────────────────────────────────────────────────────────────────────

ATRT_RESISTANCE_BYPASS_MAP: Dict[str, List[str]] = {
    "EZH2":   ["BRD4", "CDK4", "MYC"],       # EZH2 resistance via super-enhancer bypass
    "BRD4":   ["CDK4", "MYC", "MYCN"],        # BET resistance via CDK-mediated transcription
    "BRD2":   ["MYC", "CDK4"],
    "BRD3":   ["MYC"],
    "HDAC1":  ["BCL2", "BCL2L1", "MCL1"],
    "HDAC2":  ["BCL2", "MCL1"],
    "CDK4":   ["PIK3CA", "MTOR", "MTOR"],
    "CDK6":   ["PIK3CA", "MTOR"],
    "AURKA":  ["CDK4", "MYC", "BCL2L1"],
    "MTOR":   ["PIK3CA", "AKT1"],
    "PIK3CA": ["MTOR", "AKT1"],
    "PSMB5":  ["ATG5", "BECN1", "MCL1"],
    "MYC":    ["AURKA", "CDK4", "BCL2L1"],
    "BCL2":   ["MCL1", "BCL2L1"],
    "PARP1":  ["RAD51", "BRCA1", "ATM"],
    "SMO":    ["GLI2", "PIK3CA"],             # SHH resistance via PI3K
    "GLI2":   ["MYC", "CDK4"],
}

# Constitutive resistance nodes in ATRT
ATRT_CONSTITUTIVE_RESISTANCE: Set[str] = {
    "MYC", "MYCN", "BCL2L1", "MCL1", "CDK4", "BRD4", "PIK3CA",
}


# ─────────────────────────────────────────────────────────────────────────────
# MAIN SCORER CLASS
# ─────────────────────────────────────────────────────────────────────────────

class ATRTSpecializedScorer:
    """
    ATRT-specific scoring augmentation.

    Key differences from DIPGSpecializedScorer:
    1. EZH2 inhibitors are BOOSTED (not penalised) — synthetic lethality with SMARCB1 loss
    2. AURKA inhibitors boosted — especially in MYC subgroup
    3. SMO/GLI inhibitors boosted in SHH subgroup
    4. BBB penalty is LESS severe than DIPG (ATRT is not exclusively brainstem)
    5. Subgroup stratification (TYR/SHH/MYC) affects drug priorities
    6. No H3K27M scoring — SMARCB1 loss is the driver

    Sources:
    - Torchia 2015, Johann 2016 (subgroups)
    - Knutson 2013 (EZH2 synthetic lethality)
    - Frühwald 2020 (ATRT overview)
    - Sredni 2017 (AURKA)
    """

    MECHANISM_KEYWORDS: Dict[str, float] = {
        # EZH2/PRC2 — HIGHEST priority (synthetic lethality with SMARCB1 loss)
        "ezh2 inhibitor":               0.50,
        "prc2 inhibitor":               0.48,
        "histone methyltransferase":    0.42,
        # BET bromodomain
        "bet inhibitor":                0.38,
        "brd4 inhibitor":               0.38,
        "bromodomain":                  0.32,
        # HDAC
        "hdac inhibitor":               0.38,
        "pan-hdac":                     0.40,
        "histone deacetylase":          0.38,
        # Aurora kinase
        "aurora kinase":                0.35,
        "aurka":                        0.38,
        # CDK4/6
        "cdk4":                         0.30,
        "cdk4/6":                       0.32,
        # Proteasome
        "proteasome inhibitor":         0.28,
        "proteasome":                   0.25,
        # SHH (subgroup-specific — scored separately)
        "smoothened":                   0.20,
        "hedgehog":                     0.22,
        "smo inhibitor":                0.22,
        # mTOR/PI3K
        "mtor inhibitor":               0.22,
        "pi3k inhibitor":               0.20,
        # PARP
        "parp inhibitor":               0.20,
    }

    def __init__(
        self,
        subgroup: Optional[str] = None,  # "TYR", "SHH", "MYC", or None (pan-ATRT)
        novelty_bonus: float = 0.08,
        smarcb1_synleth_bonus: float = 0.15,  # Extra bonus for synthetic lethality targets
        apply_bbb_penalty: bool = True,
    ):
        self.subgroup = subgroup
        self.novelty_bonus = novelty_bonus
        self.smarcb1_synleth_bonus = smarcb1_synleth_bonus
        self.apply_bbb_penalty = apply_bbb_penalty

        logger.info(
            "ATRTSpecializedScorer: subgroup=%s, smarcb1_synleth_bonus=%.2f",
            subgroup or "pan-ATRT", smarcb1_synleth_bonus,
        )

    def _is_ezh2_inhibitor(self, drug: Dict) -> bool:
        name_lower = (drug.get("name") or drug.get("drug_name") or "").lower().strip()
        mech_lower = (drug.get("mechanism") or "").lower()
        targets = [t.upper() for t in (drug.get("targets") or [])]

        for suffix in (" hydrochloride", " hcl", " sodium", " mesylate", " phosphate"):
            name_lower = name_lower.replace(suffix, "")

        if name_lower.strip() in EZH2_INHIBITOR_ATRT["known_inhibitors"]:
            return True
        if any(kw in mech_lower for kw in EZH2_INHIBITOR_ATRT["mechanism_keywords"]):
            return True
        if targets == ["EZH2"]:
            return True
        return False

    def _is_aurka_inhibitor(self, drug: Dict) -> bool:
        name_lower = (drug.get("name") or drug.get("drug_name") or "").lower().strip()
        mech_lower = (drug.get("mechanism") or "").lower()
        targets = [t.upper() for t in (drug.get("targets") or [])]

        if name_lower.strip() in AURKA_INHIBITOR_ATRT["known_inhibitors"]:
            return True
        if any(kw in mech_lower for kw in AURKA_INHIBITOR_ATRT["mechanism_keywords"]):
            return True
        if "AURKA" in targets:
            return True
        return False

    def _is_smo_inhibitor(self, drug: Dict) -> bool:
        name_lower = (drug.get("name") or drug.get("drug_name") or "").lower().strip()
        mech_lower = (drug.get("mechanism") or "").lower()
        targets = [t.upper() for t in (drug.get("targets") or [])]

        if name_lower.strip() in SMO_INHIBITOR_ATRT["known_inhibitors"]:
            return True
        if any(kw in mech_lower for kw in SMO_INHIBITOR_ATRT["mechanism_keywords"]):
            return True
        if "SMO" in targets or "GLI2" in targets or "GLI1" in targets:
            return True
        return False

    def _smarcb1_synleth_score(self, drug: Dict) -> float:
        """
        Score a drug's synthetic lethality potential with SMARCB1 loss.

        SMARCB1 loss creates dependency on:
        1. EZH2/PRC2 (primary) — Knutson 2013
        2. BRD4/BET (super-enhancer dependency) — Nagaraja analogy
        3. AURKA (MYCN stabilisation) — Sredni 2017
        4. Proteasome (MYC/MYCN degradation)
        5. HDAC1/2 (compensatory chromatin compaction)

        Source: Wilson & Roberts 2011 Nat Rev Cancer; Frühwald 2020 CNS Oncology.
        """
        targets = [t.upper() for t in (drug.get("targets") or [])]
        mechanism = (drug.get("mechanism") or "").lower()

        # Core synthetic lethality targets with SMARCB1 loss
        synleth_targets = {
            "EZH2":   1.00,  # Primary — Knutson 2013
            "EED":    0.90,
            "SUZ12":  0.88,
            "BRD4":   0.82,
            "BRD2":   0.78,
            "HDAC1":  0.75,
            "HDAC2":  0.72,
            "AURKA":  0.70,
            "PSMB5":  0.65,
            "MYC":    0.68,
            "MYCN":   0.65,
        }

        target_hits = {t: synleth_targets[t] for t in targets if t in synleth_targets}
        if not target_hits:
            return 0.0

        best_score = max(target_hits.values())
        score = best_score * self.smarcb1_synleth_bonus / 1.00

        return min(score, self.smarcb1_synleth_bonus)

    def _subgroup_score(self, drug: Dict) -> float:
        """Score based on ATRT molecular subgroup if known."""
        if self.subgroup is None:
            return 0.0

        subgroup_upper = self.subgroup.upper()

        if subgroup_upper == "MYC":
            # BET, AURKA, CDK4/6, proteasome — most aggressive subgroup
            if self._is_aurka_inhibitor(drug):
                return 0.10
            targets = [t.upper() for t in (drug.get("targets") or [])]
            if set(targets) & {"MYC", "MYCN", "BRD4", "CDK4", "PSMB5"}:
                return 0.08

        elif subgroup_upper == "SHH":
            # SMO/GLI inhibitors most relevant
            if self._is_smo_inhibitor(drug):
                return 0.10
            targets = [t.upper() for t in (drug.get("targets") or [])]
            if set(targets) & {"SMO", "GLI2", "GLI1", "PTCH1"}:
                return 0.08

        elif subgroup_upper == "TYR":
            # HDAC + BET most relevant; youngest patients → extra toxicity weight
            mechanism = (drug.get("mechanism") or "").lower()
            if "hdac" in mechanism or "bet" in mechanism:
                return 0.06

        return 0.0

    def _mechanism_score(self, drug: Dict) -> float:
        mechanism = (drug.get("mechanism") or "").lower()
        targets = [t.upper() for t in (drug.get("targets") or [])]

        score = 0.0
        for keyword, weight in self.MECHANISM_KEYWORDS.items():
            if keyword in mechanism:
                score = max(score, weight)

        # Target-based scoring
        atrt_vulnerability_targets = {
            "EZH2", "EED", "SUZ12",          # PRC2 synthetic lethality
            "BRD4", "BRD2", "BRD3",          # BET bromodomain
            "HDAC1", "HDAC2", "HDAC3",        # HDAC
            "AURKA", "AURKB",                 # Aurora kinase
            "CDK4", "CDK6",                   # Cell cycle
            "PSMB5", "PSMB2",                 # Proteasome
            "MYC", "MYCN",                    # MYC axis
        }

        target_hits = set(targets) & atrt_vulnerability_targets
        if target_hits:
            score = max(score, 0.30 + len(target_hits) * 0.04)

        return min(score * 0.12 / 0.50, 0.12)  # Normalise to bonus scale

    def score_candidate(self, candidate: Dict) -> Dict:
        """Score a single drug candidate for ATRT."""
        base_score = candidate.get("score", 0.0)

        # EZH2 inhibitor BOOST (not penalty — critical ATRT/DIPG difference)
        is_ezh2 = self._is_ezh2_inhibitor(candidate)
        is_aurka = self._is_aurka_inhibitor(candidate)
        is_smo = self._is_smo_inhibitor(candidate)

        if is_ezh2:
            ezh2_mult = EZH2_INHIBITOR_ATRT["composite_boost"]
            ezh2_note = EZH2_INHIBITOR_ATRT["rationale"]
            # Apply boost to base score
            adjusted = round(min(1.0, base_score * ezh2_mult), 4)
            logger.info(
                "EZH2 BOOST applied to %s: %.3f → %.3f (×%.2f). SMARCB1 synthetic lethality.",
                candidate.get("name", "?"), base_score, adjusted, ezh2_mult,
            )
        else:
            # Standard bonus pathway
            synleth = self._smarcb1_synleth_score(candidate)
            mech = self._mechanism_score(candidate)
            subgrp = self._subgroup_score(candidate)

            bbb_pen = candidate.get("bbb_penetrance", "UNKNOWN")
            cns_boost = 0.04 if bbb_pen == "HIGH" else 0.0  # Smaller than DIPG (not always brainstem)

            adjusted = round(min(1.0, base_score + synleth + mech + subgrp + cns_boost), 4)
            ezh2_mult = 1.0
            ezh2_note = ""
            synleth_note = f"SMARCB1 synleth bonus: {synleth:.3f}"

        candidate["atrt_score"] = adjusted
        candidate["atrt_components"] = {
            "base_score":               round(base_score, 4),
            "is_ezh2_inhibitor":        is_ezh2,
            "ezh2_boosted":             is_ezh2,
            "ezh2_boost_multiplier":    ezh2_mult if is_ezh2 else 1.0,
            "ezh2_rationale":           ezh2_note if is_ezh2 else "",
            "is_aurka_inhibitor":       is_aurka,
            "is_smo_inhibitor":         is_smo,
            "subgroup":                 self.subgroup or "pan-ATRT",
            "smarcb1_synleth_relevant": not is_ezh2 and (adjusted - base_score) > 0.02,
        }
        candidate["score"] = adjusted
        return candidate

    def score_batch(self, candidates: List[Dict]) -> List[Dict]:
        logger.info("🧬 ATRT specialization: scoring %d candidates (subgroup=%s)",
                    len(candidates), self.subgroup or "pan-ATRT")
        for c in candidates:
            self.score_candidate(c)
        candidates.sort(key=lambda x: x.get("score", 0), reverse=True)

        n_ezh2 = sum(1 for c in candidates if c.get("atrt_components", {}).get("ezh2_boosted"))
        n_aurka = sum(1 for c in candidates if c.get("atrt_components", {}).get("is_aurka_inhibitor"))
        n_smo = sum(1 for c in candidates if c.get("atrt_components", {}).get("is_smo_inhibitor"))
        logger.info(
            "   EZH2 inhibitors boosted: %d | AURKA inhibitors: %d | SMO inhibitors: %d",
            n_ezh2, n_aurka, n_smo,
        )
        return candidates

    def generate_subgroup_report(self, candidates: List[Dict]) -> str:
        lines = [
            "# ATRT Molecular Subgroup Drug Stratification",
            "",
            "Source: Torchia et al. 2015 Cancer Cell; Johann et al. 2016 Cancer Cell.",
            "",
            "| Subgroup | Prevalence | Key Targets | Top Drug Candidates |",
            "|----------|------------|-------------|---------------------|",
        ]

        for subgroup, prev in ATRT_SUBGROUP_PREVALENCE.items():
            if subgroup == "TYR":
                targets = "HDAC1/2, BRD4"
                drugs = [c.get("name", "?") for c in candidates
                         if any(t in (c.get("targets") or [])
                                for t in ["HDAC1", "HDAC2", "BRD4"])][:3]
            elif subgroup == "SHH":
                targets = "SMO, GLI2, EZH2, CDK4"
                drugs = [c.get("name", "?") for c in candidates
                         if any(t in (c.get("targets") or [])
                                for t in ["SMO", "GLI2", "EZH2", "CDK4"])][:3]
            else:  # MYC
                targets = "AURKA, BRD4, EZH2, PSMB5"
                drugs = [c.get("name", "?") for c in candidates
                         if any(t in (c.get("targets") or [])
                                for t in ["AURKA", "BRD4", "EZH2", "PSMB5"])][:3]

            lines.append(
                f"| {subgroup} | ~{prev:.0%} | {targets} | {', '.join(drugs) or 'N/A'} |"
            )

        lines += [
            "",
            "**NOTE**: Subgroup calling requires methylation array (EPIC 850K) or RNA-seq.",
            "Without subgroup data, report pan-ATRT scores only.",
            "",
            "## EZH2 Inhibitor Rationale in ATRT",
            "",
            "SMARCB1 normally OPPOSES PRC2/EZH2 activity. When SMARCB1 is lost:",
            "- PRC2/EZH2 becomes hyperactive and unchecked",
            "- EZH2 is now ESSENTIAL for ATRT cell survival",
            "- This is a SYNTHETIC LETHALITY relationship",
            "- Tazemetostat (EZH2 inhibitor) is FDA Breakthrough Therapy in SMARCB1-deficient tumors",
            "",
            "This is the OPPOSITE of DIPG, where H3K27M inhibits EZH2 (making EZH2 inhibitors non-rational).",
            "",
            "Source: Knutson SK et al. 2013 PNAS 110(19):7922. PMID 23620515.",
        ]
        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Helper functions
# ─────────────────────────────────────────────────────────────────────────────

def get_atrt_gene_score_weights() -> Dict[str, float]:
    return ATRT_GENE_SCORE_WEIGHTS


def get_atrt_disease_data_supplement(subgroup: Optional[str] = None) -> Dict:
    gene_weights = get_atrt_gene_score_weights()

    # Merge subgroup-specific genes if subgroup provided
    extra_genes = []
    if subgroup == "TYR":
        extra_genes = list(ATRT_TYR_GENES)
    elif subgroup == "SHH":
        extra_genes = list(ATRT_SHH_GENES)
    elif subgroup == "MYC":
        extra_genes = list(ATRT_MYC_GENES)

    all_genes = list(set(ATRT_CORE_GENES + extra_genes))

    return {
        "name":           "Atypical Teratoid/Rhabdoid Tumor",
        "aliases":        ["ATRT", "rhabdoid tumor CNS", "MRT brain"],
        "genes":          all_genes,
        "gene_scores":    gene_weights,
        "pathways":       list(ATRT_PATHWAY_WEIGHTS.keys()),
        "is_rare":        True,
        "is_pediatric":   True,
        "smarcb1_deficient": True,
        "bbb_relevant":   True,  # Variable — not always brainstem
        "subgroup":       subgroup or "pan-ATRT",
        "disease_params": {
            "baseline_orr":             0.08,
            "baseline_pfs6":            0.12,
            "smarcb1_loss_fraction":    0.95,
            "smarca4_loss_fraction":    0.05,
            "subgroup_prevalence":      ATRT_SUBGROUP_PREVALENCE,
            "bbb_barrier":              0.50,  # MODERATE — location-dependent
        },
        "source": "ATRT specialization module v1.0",
        "notes": (
            "ATRT is driven by SMARCB1/SMARCA4 loss. Three molecular subgroups: "
            "TYR (~36%), SHH (~37%), MYC (~27%) — Johann 2016. "
            "EZH2 inhibitors are strongly rational (synthetic lethality with SMARCB1 loss). "
            "OPPOSITE biology to DIPG: EZH2 is BOOSTED not penalised. "
            "Median OS ~17 months with intensive therapy."
        ),
    }


def augment_disease_data_for_atrt(
    disease_data: Dict, subgroup: Optional[str] = None
) -> Dict:
    supplement = get_atrt_disease_data_supplement(subgroup)
    augmented = dict(disease_data)

    existing_genes = set(disease_data.get("genes", []))
    atrt_genes = set(supplement["genes"])
    augmented["genes"] = list(atrt_genes | existing_genes)
    augmented["gene_scores"] = supplement["gene_scores"]
    augmented["pathways"] = supplement["pathways"]
    augmented["disease_params"] = supplement["disease_params"]
    augmented["is_pediatric"] = True
    augmented["smarcb1_deficient"] = True
    augmented["bbb_relevant"] = True
    augmented["source"] = supplement["source"]
    augmented["atrt_augmented"] = True
    augmented["subgroup"] = subgroup or "pan-ATRT"

    logger.info(
        "Disease data augmented for ATRT: %d total genes, subgroup=%s",
        len(augmented["genes"]), subgroup or "pan-ATRT",
    )
    return augmented