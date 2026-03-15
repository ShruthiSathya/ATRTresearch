import logging
from typing import Dict, List, Tuple, Optional

from .pipeline_config import BBB as BBB_CONFIG, BBB_EXTENDED_KNOWN

logger = logging.getLogger(__name__)


KNOWN_BBB_PENETRANCE: Dict[str, str] = {
    # ── HIGH ──────────────────────────────────────────────────────────────────
    "temozolomide":   "HIGH",
    "onc201":         "HIGH",
    "panobinostat":   "HIGH",
    "abemaciclib":    "HIGH",
    "dexamethasone":  "HIGH",
    "lomustine":      "HIGH",
    "carmustine":     "HIGH",
    "marizomib":      "HIGH",
    "thioridazine":   "HIGH",
    "valproic acid":  "HIGH",
    "chloroquine":    "HIGH",
    "paxalisib":      "HIGH",
    "gdc-0084":       "HIGH",

    # ── MODERATE ──────────────────────────────────────────────────────────────
    "hydroxychloroquine": "MODERATE",
    "metformin":      "MODERATE",
    "itraconazole":   "MODERATE",
    "ribociclib":     "MODERATE",
    "birabresib":     "MODERATE",
    "vorinostat":     "MODERATE",
    "onatasertib":    "MODERATE",
    "indoximod":      "MODERATE",
    "lapatinib":      "MODERATE",
    "erlotinib":      "MODERATE",
    "gefitinib":      "MODERATE",

    # ── LOW ───────────────────────────────────────────────────────────────────
    "palbociclib":    "LOW",
    "belinostat":     "LOW",
    "romidepsin":     "LOW",
    "vistusertib":    "LOW",
    "ridaforolimus":  "LOW",
    "voxtalisib":     "LOW",
    "bevacizumab":    "LOW",
    "pembrolizumab":  "LOW",
    "nivolumab":      "LOW",
    "rituximab":      "LOW",
    "trastuzumab":    "LOW",
    "cetuximab":      "LOW",
    "dasatinib":      "LOW",
    "imatinib":       "LOW",
    "bortezomib":     "LOW",
}

# FIX: Merge extended known BBB data from pipeline_config
# This fills in UNKNOWN gaps for drugs ranked in top-20 (AZD-8055, crizotinib, etc.)
KNOWN_BBB_PENETRANCE.update(BBB_EXTENDED_KNOWN)


KNOWN_GBM_FAILURES = {
    "cilengitide", "enzastaurin", "temsirolimus", "cediranib", "iniparib",
    "vorinostat", "erlotinib", "gefitinib", "imatinib", "tipifarnib",
    "sorafenib", "sunitinib", "dasatinib", "everolimus", "ribociclib",
    "palbociclib", "belinostat", "romidepsin", "ridaforolimus",
    "vistusertib", "voxtalisib",
}


class BBBFilter:
    """
    Blood-brain barrier filter and scorer.
    v5.5: UNKNOWN score changed from 0.5 → 0.4 (unknown ≠ moderate).
          Extended known list merged from pipeline_config.BBB_EXTENDED_KNOWN.
    """

    def __init__(self, hard_exclude_mw: float = None):
        self.hard_exclude_mw = hard_exclude_mw or BBB_CONFIG["hard_exclude_mw"]
        logger.info("BBB Filter initialized (MW limit: %.1f Da)", self.hard_exclude_mw)

    def score_drug(self, drug_name: str,
                   molecular_weight: Optional[float] = None) -> Dict:
        name_lower = drug_name.lower().strip()

        # Strip common salt suffixes for matching
        for suffix in (" hydrochloride", " hcl", " sodium", " mesylate",
                       " malate", " phosphate", " sulfate", " mafodotin"):
            name_lower = name_lower.replace(suffix, "")
        name_lower = name_lower.strip()

        # Hard MW exclusion
        if molecular_weight and molecular_weight > self.hard_exclude_mw:
            return {
                "penetrance": "LOW",
                "bbb_score":  BBB_CONFIG["heuristic_low_score"],
                "reason":     f"MW {molecular_weight:.0f} Da > {self.hard_exclude_mw:.0f} Da limit",
                "clinical_failure": False,
            }

        # Known GBM clinical trial failure
        if name_lower in KNOWN_GBM_FAILURES:
            penetrance = KNOWN_BBB_PENETRANCE.get(name_lower, "LOW")
            return {
                "penetrance": penetrance,
                "bbb_score":  BBB_CONFIG["failure_score"],
                "reason":     "Known GBM/CNS clinical trial failure — deprioritised",
                "clinical_failure": True,
            }

        # Curated PK database (primary + extended)
        if name_lower in KNOWN_BBB_PENETRANCE:
            penetrance = KNOWN_BBB_PENETRANCE[name_lower]
            return {
                "penetrance": penetrance,
                "bbb_score":  self._penetrance_to_score(penetrance),
                "reason":     "Curated PK database",
                "clinical_failure": False,
            }

        # MW heuristic fallback
        if molecular_weight:
            if molecular_weight < BBB_CONFIG["mw_moderate_cutoff"]:
                return {"penetrance": "MODERATE",
                        "bbb_score": BBB_CONFIG["heuristic_moderate_score"],
                        "reason": f"MW {molecular_weight:.0f} < {BBB_CONFIG['mw_moderate_cutoff']:.0f} Da (heuristic)",
                        "clinical_failure": False}
            elif molecular_weight < BBB_CONFIG["mw_low_cutoff"]:
                return {"penetrance": "LOW",
                        "bbb_score": BBB_CONFIG["heuristic_low_near_score"],
                        "reason": f"MW {molecular_weight:.0f} in 400-600 Da range (heuristic)",
                        "clinical_failure": False}
            else:
                return {"penetrance": "LOW",
                        "bbb_score": BBB_CONFIG["heuristic_low_score"],
                        "reason": f"MW {molecular_weight:.0f} > {BBB_CONFIG['mw_low_cutoff']:.0f} Da (heuristic)",
                        "clinical_failure": False}

        # FIX: UNKNOWN score is now 0.4 (from config), not 0.5
        # Unknown BBB ≠ MODERATE — lack of data should be penalised mildly.
        return {
            "penetrance": "UNKNOWN",
            "bbb_score":  BBB_CONFIG["unknown_score"],
            "reason":     "No PK data available",
            "clinical_failure": False,
        }

    def filter_and_rank(self, candidates: List[Dict],
                        apply_penalty: bool = True,
                        exclude_low: bool = False) -> Tuple[List[Dict], List[Dict]]:
        passing, excluded = [], []
        n_failures = 0

        for c in candidates:
            name   = c.get("name") or c.get("drug_name") or "Unknown"
            mw     = c.get("molecular_weight")
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
            logger.info("BBB filter: penalised %d known GBM clinical trial failures", n_failures)

        return passing, excluded

    def _penetrance_to_score(self, category: str) -> float:
        return BBB_CONFIG["penetrance_scores"].get(category, BBB_CONFIG["unknown_score"])