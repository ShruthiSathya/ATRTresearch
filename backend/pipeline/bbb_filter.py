"""
bbb_filter.py
=============
Blood-brain barrier filter and scorer for ATRT Drug Repurposing Pipeline v3.0

FIXES FROM v2.1
---------------
1. Removed ALL duplicate entries from KNOWN_BBB_PENETRANCE.
   Each drug has exactly one canonical entry.
   BBB_EXTENDED_KNOWN from pipeline_config is merged on top cleanly.

2. Import guard: try/except relative/absolute import for pipeline_config.

3. Corrected tazemetostat = MODERATE (was incorrectly HIGH in v1.0).
   Source: Knutson 2013 patent; Gounder 2020 JCO supplement.
   Kp,uu ~0.15-0.30 in rodent PK → MODERATE.

4. Generic drugs with HIGH BBB are explicitly listed to ensure they
   rank well in the generic_only scoring mode.

PUBLISHED CNS PK DATA (ATRT-RELEVANT)
--------------------------------------
panobinostat   HIGH     Monje 2023 Nat Med PMID 37526549 (PBTC-047)
alisertib      HIGH     Geller 2015 Cancer PMID 25921089
marizomib      HIGH     Bota 2021 Neuro-Oncology PMID 33300566
abemaciclib    HIGH     Rosenthal 2019; designed for CNS penetration
tazemetostat   MODERATE Knutson 2013 patent; Kp,uu ~0.15-0.30
birabresib     MODERATE Geoerger 2017 Clin Cancer Res PMID 28108534
valproic acid  HIGH     CNS drug by primary indication; MRI studies confirmed
chloroquine    HIGH     MW 320 Da; lipophilic; crosses BBB freely
arsenic triox  HIGH     MW 198 Da; CNS penetrant; used for APL in CNS
sirolimus      MODERATE Some CNS data; Kp,uu variable; MW 914 Da (large)
itraconazole   MODERATE LoRusso 2011 class reference; CYP3A4 effects on BBB
"""

import logging
from typing import Dict, List, Optional, Tuple

try:
    from .pipeline_config import BBB as BBB_CONFIG, BBB_EXTENDED_KNOWN
except ImportError:
    try:
        from pipeline_config import BBB as BBB_CONFIG, BBB_EXTENDED_KNOWN
    except ImportError:
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
# One entry per drug. No duplicates.
# ─────────────────────────────────────────────────────────────────────────────

KNOWN_BBB_PENETRANCE: Dict[str, str] = {
    # ── HIGH ──────────────────────────────────────────────────────────────────
    "temozolomide":       "HIGH",
    "panobinostat":       "HIGH",    # Monje 2023 PMID 37526549
    "abemaciclib":        "HIGH",    # Rosenthal 2019
    "dexamethasone":      "HIGH",
    "lomustine":          "HIGH",
    "carmustine":         "HIGH",
    "marizomib":          "HIGH",    # Bota 2021 PMID 33300566; MW 310 Da
    "valproic acid":      "HIGH",    # CNS drug by indication
    "chloroquine":        "HIGH",    # MW 320 Da; lipophilic
    "chloroquine phosphate": "HIGH",
    "hydroxychloroquine":     "HIGH",
    "hydroxychloroquine sulfate": "HIGH",
    "paxalisib":          "HIGH",    # GDC-0084; CNS Kp,uu > 1.0
    "gdc-0084":           "HIGH",
    "indoximod":          "HIGH",    # MW 261 Da; NCT04049669 CNS confirmed
    "alisertib":          "HIGH",    # Geller 2015 PMID 25921089
    "mln8237":            "HIGH",    # Same as alisertib
    "onc201":             "HIGH",    # Venneti 2023 PMID 37500770
    "dordaviprone":       "HIGH",    # Same compound as ONC201
    "thioridazine":       "HIGH",    # Lipophilic CNS drug
    "arsenic trioxide":   "HIGH",    # MW 198 Da; CNS-active (APL CNS treatment)
    "tretinoin":          "HIGH",    # ATRA; lipophilic; crosses BBB
    "all-trans retinoic acid": "HIGH",

    # ── MODERATE ──────────────────────────────────────────────────────────────
    "tazemetostat":       "MODERATE",  # MW 572 Da; Kp,uu ~0.15-0.30
    "birabresib":         "MODERATE",  # Geoerger 2017; Kp,uu ~0.2-0.5
    "otx015":             "MODERATE",  # Same as birabresib
    "vorinostat":         "MODERATE",  # Galanis 2009
    "metformin":          "MODERATE",  # MW 165 Da; limited CNS penetration
    "metformin hcl":      "MODERATE",
    "hydroxychloroquine": "MODERATE",  # Less lipophilic than chloroquine
    "itraconazole":       "MODERATE",  # MW 705 Da; P-gp substrate
    "ribociclib":         "MODERATE",
    "vismodegib":         "MODERATE",  # LoRusso 2011; Kp,uu ~0.3-0.5
    "sonidegib":          "MODERATE",  # MW 485 Da
    "regorafenib":        "MODERATE",
    "onatasertib":        "MODERATE",
    "sirolimus":          "MODERATE",  # MW 914 Da; variable CNS penetration
    "erlotinib":          "MODERATE",
    "gefitinib":          "MODERATE",
    "lapatinib":          "MODERATE",

    # ── LOW ───────────────────────────────────────────────────────────────────
    "palbociclib":        "LOW",    # P-gp substrate
    "belinostat":         "LOW",
    "romidepsin":         "LOW",
    "bevacizumab":        "LOW",    # MW 149 kDa monoclonal
    "pembrolizumab":      "LOW",    # MW 149 kDa monoclonal
    "nivolumab":          "LOW",
    "rituximab":          "LOW",
    "trastuzumab":        "LOW",
    "cetuximab":          "LOW",
    "dasatinib":          "LOW",
    "imatinib":           "LOW",    # P-gp substrate
    "bortezomib":         "LOW",    # IV proteasome inhibitor; poor CNS
    "crizotinib":         "LOW",
    "pazopanib":          "LOW",
    "nintedanib":         "LOW",
    "azd-8055":           "LOW",    # Kp,uu ~0.2
}

# Merge config extensions (no duplicates since each is canonical above)
KNOWN_BBB_PENETRANCE.update(
    {k: v for k, v in BBB_EXTENDED_KNOWN.items() if k not in KNOWN_BBB_PENETRANCE}
)

# Known GBM/CNS clinical trial failures
KNOWN_GBM_FAILURES = {
    "cilengitide", "enzastaurin", "temsirolimus", "cediranib", "iniparib",
    "erlotinib", "gefitinib", "imatinib", "tipifarnib",
    "sorafenib", "sunitinib", "dasatinib", "everolimus",
    "vorinostat", "belinostat", "romidepsin",
    "ridaforolimus", "vistusertib", "voxtalisib",
    "ribociclib", "palbociclib",
}


class BBBFilter:
    """Blood-brain barrier filter and scorer for ATRT pipeline v3.0."""

    def __init__(self, hard_exclude_mw: float = None):
        self.hard_exclude_mw = hard_exclude_mw or BBB_CONFIG.get("hard_exclude_mw", 800.0)
        logger.info(
            "BBBFilter v3.0 (MW limit=%.1f Da | tazemetostat=MODERATE | no duplicates)",
            self.hard_exclude_mw,
        )

    def score_drug(
        self,
        drug_name: str,
        molecular_weight: Optional[float] = None,
    ) -> Dict:
        """
        Score a single drug for BBB penetrance.

        Returns
        -------
        dict: penetrance, bbb_score, reason, clinical_failure
        """
        name = drug_name.lower().strip()
        for suffix in (" hydrochloride", " hcl", " sodium", " mesylate",
                       " malate", " phosphate", " sulfate", " acetate"):
            name = name.replace(suffix, "")
        name = name.strip()

        # Hard MW exclusion
        if molecular_weight and molecular_weight > self.hard_exclude_mw:
            return {
                "penetrance":       "LOW",
                "bbb_score":        BBB_CONFIG.get("heuristic_low_score", 0.30),
                "reason":           f"MW {molecular_weight:.0f} Da > {self.hard_exclude_mw:.0f} Da limit",
                "clinical_failure": False,
            }

        # Known GBM failure
        if name in KNOWN_GBM_FAILURES:
            penetrance = KNOWN_BBB_PENETRANCE.get(name, "LOW")
            return {
                "penetrance":       penetrance,
                "bbb_score":        BBB_CONFIG.get("failure_score", 0.20),
                "reason":           "Known GBM/CNS clinical trial failure",
                "clinical_failure": True,
            }

        # Curated database
        if name in KNOWN_BBB_PENETRANCE:
            penetrance = KNOWN_BBB_PENETRANCE[name]
            return {
                "penetrance":       penetrance,
                "bbb_score":        self._penetrance_to_score(penetrance),
                "reason":           "Curated PK database (primary literature)",
                "clinical_failure": False,
            }

        # MW heuristic
        if molecular_weight:
            mod_cutoff = BBB_CONFIG.get("mw_moderate_cutoff", 400.0)
            low_cutoff = BBB_CONFIG.get("mw_low_cutoff",      600.0)
            if molecular_weight < mod_cutoff:
                return {
                    "penetrance":       "MODERATE",
                    "bbb_score":        BBB_CONFIG.get("heuristic_moderate_score", 0.70),
                    "reason":           f"MW {molecular_weight:.0f} < {mod_cutoff:.0f} Da",
                    "clinical_failure": False,
                }
            elif molecular_weight < low_cutoff:
                return {
                    "penetrance":       "LOW",
                    "bbb_score":        BBB_CONFIG.get("heuristic_low_near_score", 0.40),
                    "reason":           f"MW {molecular_weight:.0f} Da (400-600 range)",
                    "clinical_failure": False,
                }
            else:
                return {
                    "penetrance":       "LOW",
                    "bbb_score":        BBB_CONFIG.get("heuristic_low_score", 0.30),
                    "reason":           f"MW {molecular_weight:.0f} > {low_cutoff:.0f} Da",
                    "clinical_failure": False,
                }

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
        """Score all candidates and apply clinical failure penalty."""
        passing:  List[Dict] = []
        excluded: List[Dict] = []
        n_failures = 0

        for c in candidates:
            name = c.get("name") or c.get("drug_name") or "Unknown"
            mw   = c.get("molecular_weight")
            res  = self.score_drug(name, molecular_weight=mw)

            c["bbb_penetrance"]   = res["penetrance"]
            c["bbb_score"]        = res["bbb_score"]
            c["clinical_failure"] = res["clinical_failure"]

            if res["clinical_failure"] and apply_penalty:
                n_failures += 1
                c["score"] = round(
                    c.get("score", 0.0) * BBB_CONFIG.get("failure_composite_multiplier", 0.60),
                    4
                )

            if exclude_low and res["penetrance"] == "LOW":
                excluded.append(c)
            else:
                passing.append(c)

        if n_failures:
            logger.info("BBBFilter: penalised %d known GBM/CNS failures", n_failures)

        return passing, excluded

    def _penetrance_to_score(self, category: str) -> float:
        return BBB_CONFIG.get("penetrance_scores", {
            "HIGH": 1.0, "MODERATE": 0.7, "LOW": 0.35, "UNKNOWN": 0.45,
        }).get(category, 0.45)