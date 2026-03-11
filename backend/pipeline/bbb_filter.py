import logging
from typing import Dict, List, Tuple, Optional

from .pipeline_config import BBB as BBB_CONFIG

logger = logging.getLogger(__name__)


KNOWN_BBB_PENETRANCE: Dict[str, str] = {
    # ── HIGH ──────────────────────────────────────────────────────────────────
    "temozolomide":   "HIGH",
    "onc201":         "HIGH",
    "panobinostat":   "HIGH",    # PBTC-047 confirmed CNS exposure (Monje 2023)
    "abemaciclib":    "HIGH",    # Designed for CNS penetrance vs palbociclib
    "dexamethasone":  "HIGH",
    "lomustine":      "HIGH",
    "carmustine":     "HIGH",
    "marizomib":      "HIGH",    # Engineered for CNS — crosses BBB unlike bortezomib
    "thioridazine":   "HIGH",
    "valproic acid":  "HIGH",
    "chloroquine":    "HIGH",
    # Paxalisib: designed as brain-penetrant PI3K/mTOR inhibitor (Genentech).
    # FDA Orphan Drug for DIPG (2020). Pediatric Phase I (NCT03696355) and
    # Phase II (NCT05009992). CNS penetrance by design and preclinical PK.
    "paxalisib":      "HIGH",
    "gdc-0084":       "HIGH",

    # ── MODERATE ──────────────────────────────────────────────────────────────
    "hydroxychloroquine": "MODERATE",
    "metformin":      "MODERATE",
    "itraconazole":   "MODERATE",
    "ribociclib":     "MODERATE",
    "birabresib":     "MODERATE",  # PBTC-049 confirmed brain exposure (Geoerger 2017)
    "vorinostat":     "MODERATE",
    "onatasertib":    "MODERATE",  # mTOR kinase; class analogue AZD-8055 shows CNS exposure
    "indoximod":      "MODERATE",  # IDO inhibitor; pediatric DIPG trial ongoing NCT04049669
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
    All numeric thresholds read from pipeline_config.BBB — no magic numbers.
    """

    def __init__(self, hard_exclude_mw: float = None):
        self.hard_exclude_mw = hard_exclude_mw or BBB_CONFIG["hard_exclude_mw"]
        logger.info("BBB Filter initialized (MW limit: %.1f Da)", self.hard_exclude_mw)

    def score_drug(self, drug_name: str,
                   molecular_weight: Optional[float] = None) -> Dict:
        name_lower = drug_name.lower().strip()

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

        # Curated PK database
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