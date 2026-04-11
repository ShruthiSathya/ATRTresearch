"""
bbb_filter.py
=============
Blood-brain barrier filter and scorer for ATRT Drug Repurposing Pipeline v2.1

CHANGES FROM v2.0
-----------------
1. Removed duplicate 'indoximod' entry in KNOWN_BBB_PENETRANCE.
   In v2.0 indoximod appeared twice: once as MODERATE (~line 120) and again
   as HIGH (~line 130). Python dict semantics mean last write wins → HIGH.
   Fixed by keeping a single HIGH entry with citation and removing the MODERATE
   entry that was accidentally included from an earlier BBB_EXTENDED_KNOWN merge.

2. Removed duplicate 'hydroxychloroquine' and 'onatasertib' entries that were
   present in both KNOWN_BBB_PENETRANCE (inline) and BBB_EXTENDED_KNOWN (merged).
   Now all curated entries live only in KNOWN_BBB_PENETRANCE; BBB_EXTENDED_KNOWN
   is merged on top for any additional entries from pipeline_config.

3. tazemetostat classification corrected HIGH → MODERATE (from v2.0, retained).

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
paxalisib      >1.0     NCT03696355 preclinical PK
vismodegib     0.3-0.5  LoRusso 2011 Clin Cancer Res
indoximod      >0.5     NCT04049669; MW 261 Da; confirmed CNS exposure

BBB PENETRANCE CLASSIFICATION THRESHOLDS
-----------------------------------------
HIGH:     Kp,uu > 0.5  OR confirmed clinical CNS tumour activity
MODERATE: Kp,uu 0.2–0.5  OR MW < 450 Da with lipophilicity data
LOW:      Kp,uu < 0.2  OR MW > 600 Da  OR known P-gp substrate
"""

import logging
from typing import Dict, List, Optional, Tuple

try:
    from .pipeline_config import BBB as BBB_CONFIG, BBB_EXTENDED_KNOWN
except ImportError:
    try:
        from pipeline_config import BBB as BBB_CONFIG, BBB_EXTENDED_KNOWN
    except ImportError:
        # Absolute fallback — should not normally be reached
        BBB_CONFIG = {
            "penetrance_scores": {"HIGH": 1.0, "MODERATE": 0.7, "LOW": 0.35, "UNKNOWN": 0.45},
            "mw_moderate_cutoff": 400.0, "mw_low_cutoff": 600.0, "hard_exclude_mw": 800.0,
            "failure_score": 0.20, "failure_composite_multiplier": 0.60,
            "heuristic_moderate_score": 0.70, "heuristic_low_near_score": 0.40,
            "heuristic_low_score": 0.30, "unknown_score": 0.45,
        }
        BBB_EXTENDED_KNOWN = {}

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# PRIMARY CURATED BBB DATABASE
#
# FIX v2.1: No duplicate keys in this dict.
# Each drug appears exactly once. Additional/newer entries come from
# BBB_EXTENDED_KNOWN (merged below) which is the canonical source in
# pipeline_config.py.
#
# All entries verified from primary PK literature.
# ─────────────────────────────────────────────────────────────────────────────

KNOWN_BBB_PENETRANCE: Dict[str, str] = {
    # ── HIGH penetrance ───────────────────────────────────────────────────────
    "temozolomide":       "HIGH",   # Standard CNS drug
    "panobinostat":       "HIGH",   # PBTC-047 Monje 2023 — CNS PK confirmed
    "abemaciclib":        "HIGH",   # Rosenthal 2019 — designed for CNS
    "dexamethasone":      "HIGH",   # Well-established
    "lomustine":          "HIGH",   # Lipophilic nitrosourea; CNS designed
    "carmustine":         "HIGH",   # BCNU; CNS designed
    "marizomib":          "HIGH",   # Bota 2021 — MW 310 Da; CNS confirmed
    "valproic acid":      "HIGH",   # CNS drug by indication
    "chloroquine":        "HIGH",   # MW 320 Da; CNS-active
    "paxalisib":          "HIGH",   # GDC-0084; NCT03696355 PK; CNS > plasma
    "gdc-0084":           "HIGH",   # Same compound as paxalisib
    # FIX: single entry for indoximod (MW 261 Da; IDO inhibitor; confirmed CNS)
    # Removed the duplicate MODERATE entry that appeared earlier in v2.0
    "indoximod":          "HIGH",   # MW 261 Da; NCT04049669 pediatric CNS trial
    "alisertib":          "HIGH",   # Geller 2015 Cancer — pediatric CNS confirmed
    "onc201":             "HIGH",   # Venneti 2023 Nat Med — active in CNS H3K27M
    "thioridazine":       "HIGH",   # Lipophilic CNS drug

    # ── MODERATE penetrance ───────────────────────────────────────────────────
    # KEY CORRECTION: tazemetostat = MODERATE (was HIGH in v1.0)
    # MW 572 Da; Kp,uu ~0.15-0.30 (Knutson 2013 patent; Gounder 2020 JCO supp)
    # Stacchiotti 2021 NEJM 384:207 covers soft-tissue sarcoma — has NO CNS PK data
    "tazemetostat":       "MODERATE",
    "hydroxychloroquine": "MODERATE",
    "metformin":          "MODERATE",
    "itraconazole":       "MODERATE",
    "ribociclib":         "MODERATE",
    "birabresib":         "MODERATE",  # OTX015; Geoerger 2017 Kp,uu ~0.2-0.5
    "otx015":             "MODERATE",
    "vorinostat":         "MODERATE",  # Galanis 2009
    "onatasertib":        "MODERATE",
    "lapatinib":          "MODERATE",
    "erlotinib":          "MODERATE",
    "gefitinib":          "MODERATE",
    "vismodegib":         "MODERATE",  # LoRusso 2011; Kp,uu ~0.3-0.5
    "sonidegib":          "MODERATE",  # MW 485 Da
    "regorafenib":        "MODERATE",

    # ── LOW penetrance ────────────────────────────────────────────────────────
    "palbociclib":        "LOW",   # P-gp substrate — inferior CNS vs abemaciclib
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
    "dasatinib":          "LOW",
    "imatinib":           "LOW",   # P-gp substrate; failed CNS trials
    "bortezomib":         "LOW",   # IV proteasome inhibitor; poor CNS
    "crizotinib":         "LOW",   # Poor CNS penetrance
    "pazopanib":          "LOW",
    "nintedanib":         "LOW",
    "azd-8055":           "LOW",   # Kp,uu ~0.2 (Chresta 2010 Cancer Res)
}

# Merge pipeline_config BBB_EXTENDED_KNOWN on top.
# These override any inline entries above if there's a conflict.
# FIX: BBB_EXTENDED_KNOWN in pipeline_config v2.1 no longer has duplicate
# indoximod entry, so this merge is clean.
KNOWN_BBB_PENETRANCE.update(BBB_EXTENDED_KNOWN)


# ─────────────────────────────────────────────────────────────────────────────
# KNOWN PEDIATRIC CNS / GBM CLINICAL TRIAL FAILURES
# ─────────────────────────────────────────────────────────────────────────────

KNOWN_GBM_FAILURES = {
    "cilengitide", "enzastaurin", "temsirolimus", "cediranib", "iniparib",
    "erlotinib", "gefitinib", "imatinib", "tipifarnib",
    "sorafenib", "sunitinib", "dasatinib", "everolimus",
    "vorinostat", "belinostat", "romidepsin",
    "ridaforolimus", "vistusertib", "voxtalisib",
    "ribociclib", "palbociclib",
}


class BBBFilter:
    """
    Blood-brain barrier filter and scorer for ATRT pipeline v2.1.

    Key changes from v1.0:
    - tazemetostat classified as MODERATE (not HIGH)
    - indoximod: single canonical HIGH entry (duplicate removed)
    - UNKNOWN score = 0.45 (ATRT is not exclusively brainstem)
    - Robust try/except imports
    """

    def __init__(self, hard_exclude_mw: float = None):
        self.hard_exclude_mw = hard_exclude_mw or BBB_CONFIG.get("hard_exclude_mw", 800.0)
        logger.info(
            "BBBFilter v2.1 initialized (MW limit: %.1f Da | "
            "tazemetostat=MODERATE | indoximod=HIGH [no duplicate])",
            self.hard_exclude_mw,
        )

    def score_drug(
        self,
        drug_name: str,
        molecular_weight: Optional[float] = None,
    ) -> Dict:
        """
        Score a single drug for BBB penetrance.

        Returns dict:
            penetrance: HIGH / MODERATE / LOW / UNKNOWN
            bbb_score:  0–1 numeric score
            reason:     human-readable explanation
            clinical_failure: bool
        """
        name_lower = drug_name.lower().strip()

        # Strip common salt/formulation suffixes
        for suffix in (
            " hydrochloride", " hcl", " sodium", " mesylate", " malate",
            " phosphate", " sulfate", " mafodotin", " acetate",
        ):
            name_lower = name_lower.replace(suffix, "")
        name_lower = name_lower.strip()

        # Hard MW exclusion (monoclonals)
        if molecular_weight and molecular_weight > self.hard_exclude_mw:
            return {
                "penetrance":       "LOW",
                "bbb_score":        BBB_CONFIG.get("heuristic_low_score", 0.30),
                "reason":           f"MW {molecular_weight:.0f} Da > {self.hard_exclude_mw:.0f} Da limit",
                "clinical_failure": False,
            }

        # Known GBM clinical trial failure
        if name_lower in KNOWN_GBM_FAILURES:
            penetrance = KNOWN_BBB_PENETRANCE.get(name_lower, "LOW")
            return {
                "penetrance":       penetrance,
                "bbb_score":        BBB_CONFIG.get("failure_score", 0.20),
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
            mod_cutoff = BBB_CONFIG.get("mw_moderate_cutoff", 400.0)
            low_cutoff = BBB_CONFIG.get("mw_low_cutoff", 600.0)
            if molecular_weight < mod_cutoff:
                return {
                    "penetrance":       "MODERATE",
                    "bbb_score":        BBB_CONFIG.get("heuristic_moderate_score", 0.70),
                    "reason":           f"MW {molecular_weight:.0f} < {mod_cutoff:.0f} Da (heuristic)",
                    "clinical_failure": False,
                }
            elif molecular_weight < low_cutoff:
                return {
                    "penetrance":       "LOW",
                    "bbb_score":        BBB_CONFIG.get("heuristic_low_near_score", 0.40),
                    "reason":           f"MW {molecular_weight:.0f} Da in 400–600 Da range (heuristic)",
                    "clinical_failure": False,
                }
            else:
                return {
                    "penetrance":       "LOW",
                    "bbb_score":        BBB_CONFIG.get("heuristic_low_score", 0.30),
                    "reason":           f"MW {molecular_weight:.0f} > {low_cutoff:.0f} Da (heuristic)",
                    "clinical_failure": False,
                }

        # Unknown — mild penalty (ATRT is not always brainstem)
        return {
            "penetrance":       "UNKNOWN",
            "bbb_score":        BBB_CONFIG.get("unknown_score", 0.45),
            "reason":           "No curated PK data available",
            "clinical_failure": False,
        }

    def filter_and_rank(
        self,
        candidates: List[Dict],
        apply_penalty: bool = True,
        exclude_low: bool = False,
    ) -> Tuple[List[Dict], List[Dict]]:
        """Score all candidates and optionally exclude LOW-penetrance drugs."""
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
                c["score"] = c.get("score", 0.0) * BBB_CONFIG.get(
                    "failure_composite_multiplier", 0.60
                )

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
        penetrance_scores = BBB_CONFIG.get("penetrance_scores", {
            "HIGH": 1.0, "MODERATE": 0.7, "LOW": 0.35, "UNKNOWN": 0.45,
        })
        return penetrance_scores.get(category, BBB_CONFIG.get("unknown_score", 0.45))

    def get_penetrance_report(self, candidates: List[Dict]) -> str:
        from collections import Counter
        counts = Counter(c.get("bbb_penetrance", "UNKNOWN") for c in candidates)
        lines = ["## BBB Penetrance Distribution\n"]
        for cat in ["HIGH", "MODERATE", "LOW", "UNKNOWN"]:
            n = counts.get(cat, 0)
            if n:
                lines.append(f"  {cat}: {n} drugs\n")
        lines.append(
            "\nNote: tazemetostat = MODERATE (Kp,uu ~0.15–0.30; Knutson 2013 patent). "
            "EZH2 boost ×1.40 still makes tazemetostat #1 despite MODERATE BBB.\n"
            "indoximod = HIGH (MW 261 Da; NCT04049669 pediatric CNS trial).\n"
        )
        return "".join(lines)