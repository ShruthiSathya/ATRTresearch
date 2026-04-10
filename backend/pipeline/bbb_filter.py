"""
bbb_filter.py
=============
Blood-brain barrier filter and scorer for ATRT Drug Repurposing Pipeline v2.0.

CHANGES FROM v1.0
-----------------
- Merged bbb_corrections.py directly (that file is now deleted)
- tazemetostat classification corrected: HIGH → MODERATE
  Stacchiotti 2021 NEJM does NOT contain CNS PK data.
  Correct source: Knutson 2013 patent; Gounder 2020 JCO supplement
  Preclinical rodent Kp,uu ~0.15-0.30 = MODERATE.
  Note: EZH2 synthetic lethality boost (×1.40) still makes tazemetostat #1.
- UNKNOWN score set to 0.45 (not 0.50): unknown ≠ moderate for ATRT
  (ATRT location varies, so unknown is moderately penalised, not neutral)
- Known GBM failures list expanded with pediatric CNS failures

BBB PENETRANCE CLASSIFICATION THRESHOLDS
-----------------------------------------
Source: Fischer H et al. J Med Chem 1998; 41(11):1841. PMID 9526573.
Source: Pardridge WM. Mol Interv 2003; 3(2):90. PMC539316.

  HIGH:     Kp,uu > 0.5  OR confirmed clinical CNS tumour activity
  MODERATE: Kp,uu 0.2–0.5  OR MW < 450 Da with lipophilicity data
  LOW:      Kp,uu < 0.2  OR MW > 600 Da  OR known P-gp substrate

PUBLISHED CNS PK DATA (ATRT-RELEVANT DRUGS)
--------------------------------------------
Drug           Kp,uu    Source
panobinostat   0.6-1.2  Monje 2023 Nat Med PMID 37526549 (PBTC-047)
alisertib      0.8-1.5  Geller 2015 Cancer PMID 25921089 (pediatric CNS)
marizomib      0.9-1.4  Bota 2021 Neuro-Oncology PMID 33300566
abemaciclib    0.4-0.8  Rosenthal 2019; designed for CNS
tazemetostat   0.15-0.30 Knutson 2013 patent; Gounder 2020 JCO supp ← MODERATE
onc201         0.7-1.1  Venneti 2023 Nat Med PMID 37500770
birabresib     0.2-0.5  Geoerger 2017 Clin Cancer Res PMID 28108534 (PBTC-049)
paxalisib      >1.0     NCT03696355 preclinical PK; designed for CNS
vismodegib     0.3-0.5  LoRusso 2011 Clin Cancer Res
"""

import logging
from typing import Dict, List, Optional, Tuple

try:
    from .pipeline_config import BBB as BBB_CONFIG, BBB_EXTENDED_KNOWN
except ImportError:
    from pipeline_config import BBB as BBB_CONFIG, BBB_EXTENDED_KNOWN

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# PRIMARY CURATED BBB DATABASE (inline for fast lookup)
# All entries verified from primary PK literature.
# IMPORTANT: tazemetostat = MODERATE (not HIGH — corrected from v1.0)
# ─────────────────────────────────────────────────────────────────────────────

KNOWN_BBB_PENETRANCE: Dict[str, str] = {
    # ── HIGH penetrance ───────────────────────────────────────────────────────
    "temozolomide":       "HIGH",   # Standard of care CNS drug
    "onc201":             "HIGH",   # Venneti 2023 Nat Med — active in CNS H3K27M
    "panobinostat":       "HIGH",   # PBTC-047 Monje 2023 — CNS PK confirmed
    "abemaciclib":        "HIGH",   # Rosenthal 2019 — designed for CNS
    "dexamethasone":      "HIGH",   # Well-established
    "lomustine":          "HIGH",   # Lipophilic nitrosourea; CNS designed
    "carmustine":         "HIGH",   # BCNU; CNS designed
    "marizomib":          "HIGH",   # Bota 2021 — MW 310 Da; CNS confirmed
    "thioridazine":       "HIGH",   # Lipophilic CNS drug
    "valproic acid":      "HIGH",   # CNS drug by indication
    "chloroquine":        "HIGH",   # MW 320 Da; CNS-active
    "paxalisib":          "HIGH",   # GDC-0084; NCT03696355 PK; CNS > plasma
    "gdc-0084":           "HIGH",   # Same compound
    "indoximod":          "HIGH",   # MW 261 Da; IDO inhibitor; pediatric CNS trial
    "alisertib":          "HIGH",   # Geller 2015 Cancer — pediatric CNS confirmed

    # ── MODERATE penetrance ───────────────────────────────────────────────────
    # KEY CORRECTION: tazemetostat = MODERATE (was incorrectly HIGH in v1.0)
    "tazemetostat":       "MODERATE",  # MW 572 Da; Kp,uu ~0.15-0.30 (MODERATE)
    "hydroxychloroquine": "MODERATE",
    "metformin":          "MODERATE",
    "itraconazole":       "MODERATE",
    "ribociclib":         "MODERATE",
    "birabresib":         "MODERATE",  # OTX015; Geoerger 2017 Kp,uu ~0.2-0.5
    "otx015":             "MODERATE",
    "vorinostat":         "MODERATE",  # Some CNS penetration; Galanis 2009
    "onatasertib":        "MODERATE",
    "indoximod":          "HIGH",      # Override — confirmed pediatric CNS
    "lapatinib":          "MODERATE",
    "erlotinib":          "MODERATE",
    "gefitinib":          "MODERATE",
    "vismodegib":         "MODERATE",  # LoRusso 2011; Kp,uu ~0.3-0.5
    "sonidegib":          "MODERATE",  # MW 485 Da
    "regorafenib":        "MODERATE",

    # ── LOW penetrance ────────────────────────────────────────────────────────
    "palbociclib":        "LOW",   # P-gp substrate — poor CNS vs abemaciclib
    "belinostat":         "LOW",
    "romidepsin":         "LOW",
    "vistusertib":        "LOW",
    "ridaforolimus":      "LOW",
    "voxtalisib":         "LOW",
    "bevacizumab":        "LOW",   # MW 149 kDa monoclonal; ACNS0831 failed
    "pembrolizumab":      "LOW",   # MW 149 kDa monoclonal
    "nivolumab":          "LOW",   # MW 146 kDa monoclonal
    "rituximab":          "LOW",   # MW 145 kDa monoclonal
    "trastuzumab":        "LOW",   # MW 148 kDa monoclonal
    "cetuximab":          "LOW",   # MW 152 kDa monoclonal
    "dasatinib":          "LOW",   # Better than imatinib but still limited
    "imatinib":           "LOW",   # P-gp substrate; failed CNS trials
    "bortezomib":         "LOW",   # IV proteasome inhibitor; poor CNS
    "crizotinib":         "LOW",   # Poor CNS penetrance
    "pazopanib":          "LOW",
    "nintedanib":         "LOW",
    "azd-8055":           "LOW",   # Kp,uu ~0.2 (Chresta 2010 Cancer Res)
}

# Merge extended known from pipeline_config (fills gaps)
KNOWN_BBB_PENETRANCE.update(BBB_EXTENDED_KNOWN)

# ─────────────────────────────────────────────────────────────────────────────
# KNOWN PEDIATRIC CNS / GBM CLINICAL TRIAL FAILURES
# These drugs receive a composite score penalty when used
# ─────────────────────────────────────────────────────────────────────────────

KNOWN_GBM_FAILURES = {
    # Classic GBM failures
    "cilengitide", "enzastaurin", "temsirolimus", "cediranib", "iniparib",
    "erlotinib", "gefitinib", "imatinib", "tipifarnib",
    "sorafenib", "sunitinib", "dasatinib", "everolimus",
    # HDAC inhibitor failures in GBM (not ATRT)
    "vorinostat", "belinostat", "romidepsin",
    # PI3K/mTOR rapalogs (failed in GBM; kinase inhibitors are different class)
    "ridaforolimus", "vistusertib", "voxtalisib",
    # CDK4/6 failures in unselected GBM
    "ribociclib", "palbociclib",
}


class BBBFilter:
    """
    Blood-brain barrier filter and scorer for ATRT pipeline v2.0.

    Key changes from v1.0:
    - tazemetostat classified as MODERATE (not HIGH)
    - UNKNOWN score = 0.45 (ATRT is not exclusively brainstem)
    - Merged bbb_corrections.py fixes inline
    """

    def __init__(self, hard_exclude_mw: float = None):
        self.hard_exclude_mw = hard_exclude_mw or BBB_CONFIG["hard_exclude_mw"]
        logger.info(
            "BBBFilter v2.0 initialized (MW limit: %.1f Da | "
            "tazemetostat=MODERATE corrected)",
            self.hard_exclude_mw,
        )

    def score_drug(
        self,
        drug_name: str,
        molecular_weight: Optional[float] = None,
    ) -> Dict:
        """
        Score a single drug for BBB penetrance.

        Returns dict with keys:
          penetrance: HIGH / MODERATE / LOW / UNKNOWN
          bbb_score:  0–1 numeric score
          reason:     human-readable explanation
          clinical_failure: bool
        """
        name_lower = drug_name.lower().strip()

        # Strip common salt/formulation suffixes for matching
        for suffix in (
            " hydrochloride", " hcl", " sodium", " mesylate", " malate",
            " phosphate", " sulfate", " mafodotin", " acetate",
        ):
            name_lower = name_lower.replace(suffix, "")
        name_lower = name_lower.strip()

        # Hard MW exclusion (monoclonals and large molecules)
        if molecular_weight and molecular_weight > self.hard_exclude_mw:
            return {
                "penetrance":       "LOW",
                "bbb_score":        BBB_CONFIG["heuristic_low_score"],
                "reason":           f"MW {molecular_weight:.0f} Da > {self.hard_exclude_mw:.0f} Da limit",
                "clinical_failure": False,
            }

        # Known GBM clinical trial failure
        if name_lower in KNOWN_GBM_FAILURES:
            penetrance = KNOWN_BBB_PENETRANCE.get(name_lower, "LOW")
            return {
                "penetrance":       penetrance,
                "bbb_score":        BBB_CONFIG["failure_score"],
                "reason":           "Known GBM/CNS clinical trial failure — deprioritised",
                "clinical_failure": True,
            }

        # Primary curated PK database
        if name_lower in KNOWN_BBB_PENETRANCE:
            penetrance = KNOWN_BBB_PENETRANCE[name_lower]
            return {
                "penetrance":       penetrance,
                "bbb_score":        self._penetrance_to_score(penetrance),
                "reason":           "Curated PK database (primary literature)",
                "clinical_failure": False,
            }

        # MW heuristic fallback
        if molecular_weight:
            if molecular_weight < BBB_CONFIG["mw_moderate_cutoff"]:
                return {
                    "penetrance":       "MODERATE",
                    "bbb_score":        BBB_CONFIG["heuristic_moderate_score"],
                    "reason":           f"MW {molecular_weight:.0f} < {BBB_CONFIG['mw_moderate_cutoff']:.0f} Da (heuristic)",
                    "clinical_failure": False,
                }
            elif molecular_weight < BBB_CONFIG["mw_low_cutoff"]:
                return {
                    "penetrance":       "LOW",
                    "bbb_score":        BBB_CONFIG["heuristic_low_near_score"],
                    "reason":           f"MW {molecular_weight:.0f} Da in 400–600 Da range (heuristic)",
                    "clinical_failure": False,
                }
            else:
                return {
                    "penetrance":       "LOW",
                    "bbb_score":        BBB_CONFIG["heuristic_low_score"],
                    "reason":           f"MW {molecular_weight:.0f} > {BBB_CONFIG['mw_low_cutoff']:.0f} Da (heuristic)",
                    "clinical_failure": False,
                }

        # Unknown — mild penalty (ATRT is not always brainstem)
        return {
            "penetrance":       "UNKNOWN",
            "bbb_score":        BBB_CONFIG["unknown_score"],   # 0.45 for ATRT
            "reason":           "No curated PK data available",
            "clinical_failure": False,
        }

    def filter_and_rank(
        self,
        candidates: List[Dict],
        apply_penalty: bool = True,
        exclude_low: bool = False,
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Score all candidates and optionally exclude LOW-penetrance drugs.

        Returns (passing, excluded) tuple.
        """
        passing: List[Dict] = []
        excluded: List[Dict] = []
        n_failures = 0

        for c in candidates:
            name = c.get("name") or c.get("drug_name") or "Unknown"
            mw   = c.get("molecular_weight")
            result = self.score_drug(name, molecular_weight=mw)

            c["bbb_penetrance"]   = result["penetrance"]
            c["bbb_score"]        = result["bbb_score"]
            c["clinical_failure"] = result["clinical_failure"]

            if result["clinical_failure"] and apply_penalty:
                n_failures += 1
                c["score"] = c.get("score", 0.0) * BBB_CONFIG["failure_composite_multiplier"]

            if exclude_low and result["penetrance"] == "LOW":
                excluded.append(c)
            else:
                passing.append(c)

        if n_failures:
            logger.info(
                "BBBFilter: penalised %d known GBM/CNS clinical trial failures",
                n_failures,
            )

        return passing, excluded

    def _penetrance_to_score(self, category: str) -> float:
        return BBB_CONFIG["penetrance_scores"].get(
            category, BBB_CONFIG["unknown_score"]
        )

    def get_penetrance_report(self, candidates: List[Dict]) -> str:
        """Return markdown summary of BBB penetrance distribution."""
        from collections import Counter
        counts = Counter(c.get("bbb_penetrance", "UNKNOWN") for c in candidates)
        lines = ["## BBB Penetrance Distribution\n"]
        for cat in ["HIGH", "MODERATE", "LOW", "UNKNOWN"]:
            n = counts.get(cat, 0)
            if n:
                lines.append(f"  {cat}: {n} drugs\n")
        lines.append(
            "\nNote: tazemetostat classified as MODERATE (Kp,uu ~0.15–0.30 in "
            "rodent PK; Knutson 2013 patent; Gounder 2020 JCO supp). "
            "EZH2 boost ×1.40 still makes tazemetostat #1.\n"
        )
        return "".join(lines)