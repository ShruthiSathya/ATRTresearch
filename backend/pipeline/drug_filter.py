"""
drug_filter.py
==============
Safety filter for drug candidates. Removes contraindicated drugs before
returning results to the user.

FIX: Removed the broken duplicate ProductionPipeline class that imported
non-existent modules (CMAPQueryEngine, TMEScorer, TrialOutcomeCalibrator,
DrugDiseaseGCN, etc.) and conflicted with discovery_pipeline.ProductionPipeline.
"""

import logging
from typing import Dict, List, Tuple, Optional

logger = logging.getLogger(__name__)


# Absolute contraindications: drug class → reasons
# Sources: FDA prescribing information, published safety reviews
ABSOLUTE_CONTRAINDICATIONS = {
    # Live vaccines — immunosuppressed pediatric patients
    "live vaccine": "Contraindicated in immunocompromised pediatric CNS tumor patients",
    # Strong CYP3A4 inhibitors with narrow-TI oncology drugs
}

# Drug-disease relative contraindications for pediatric CNS tumors
# Sources: COG supportive care guidelines, SIOPE ATRT guidelines
RELATIVE_CONTRAINDICATIONS: Dict[str, Dict] = {
    # Drugs with strong evidence against use in pediatric ATRT/DIPG context
    "dexamethasone": {
        "reason": "Chronic steroid use impairs immune response and growth; reserve for acute edema only",
        "severity": "relative",
        "exception": "acute CNS edema management",
    },
    "bevacizumab": {
        "reason": "Poor BBB penetrance (MW 149 kDa); no survival benefit in pediatric CNS tumors (ACNS0831)",
        "severity": "relative",
        "source": "Fangusaro 2021 Neuro-Oncology ACNS0831",
    },
    "temsirolimus": {
        "reason": "Known GBM/pediatric CNS trial failure (NCT00087815); mTOR rapalog class failed",
        "severity": "relative",
        "source": "Wick 2011 JCO",
    },
    "everolimus": {
        "reason": "mTOR rapalog; failed in adult GBM Phase 2 (Fouladi 2007); use kinase inhibitor class instead",
        "severity": "relative",
    },
    "erlotinib": {
        "reason": "Known CNS trial failure; EGFR not primary driver in ATRT",
        "severity": "relative",
        "source": "van den Bent 2009",
    },
    "gefitinib": {
        "reason": "Known CNS trial failure; EGFR not primary driver in ATRT",
        "severity": "relative",
    },
}

# Drugs known to have failed in pediatric CNS trials — should be deprioritized
# but not necessarily removed (may still be useful in combinations)
KNOWN_PEDIATRIC_CNS_FAILURES = {
    "temsirolimus", "everolimus", "erlotinib", "gefitinib",
    "imatinib", "dasatinib", "enzastaurin", "cediranib",
    "cilengitide", "tipifarnib", "vorinostat",  # vorinostat failed in GBM; ATRT data limited
    "bortezomib",  # poor BBB; IV formulation
}


class DrugSafetyFilter:
    """
    Filters drug candidates for safety contraindications in pediatric CNS tumors.

    Two pass modes:
      remove_absolute: removes drugs with absolute contraindications
      remove_relative: removes drugs with relative contraindications
                       (e.g. known trial failures, poor BBB for CNS-exclusive drug)

    Usage
    -----
        filter = DrugSafetyFilter()
        safe, filtered_out = await filter.filter_candidates(
            candidates, disease_name="atrt",
            remove_absolute=True, remove_relative=True
        )
    """

    def __init__(self):
        logger.info("DrugSafetyFilter initialized")

    async def filter_candidates(
        self,
        candidates: List[Dict],
        disease_name: str,
        remove_absolute: bool = True,
        remove_relative: bool = True,
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Filter drug candidates for safety contraindications.

        Parameters
        ----------
        candidates : list of drug dicts (must have 'drug_name' or 'name' key)
        disease_name : str — used to determine disease-specific contraindications
        remove_absolute : bool — remove absolutely contraindicated drugs
        remove_relative : bool — remove relatively contraindicated drugs

        Returns
        -------
        (safe_candidates, filtered_out)
        """
        safe: List[Dict] = []
        filtered_out: List[Dict] = []

        for candidate in candidates:
            name = (
                candidate.get("drug_name")
                or candidate.get("name")
                or ""
            ).lower().strip()

            contraindication = self._check_contraindication(
                name, disease_name, remove_absolute, remove_relative
            )

            if contraindication:
                candidate["contraindication"] = contraindication
                filtered_out.append(candidate)
                logger.debug(
                    "Filtered out %s: %s (%s)",
                    name,
                    contraindication["reason"],
                    contraindication["severity"],
                )
            else:
                safe.append(candidate)

        logger.info(
            "Safety filter: %d safe | %d filtered out (absolute=%s, relative=%s)",
            len(safe), len(filtered_out), remove_absolute, remove_relative,
        )
        return safe, filtered_out

    def _check_contraindication(
        self,
        drug_name: str,
        disease_name: str,
        remove_absolute: bool,
        remove_relative: bool,
    ) -> Optional[Dict]:
        """Return contraindication dict if drug should be filtered, else None."""

        # Check absolute contraindications
        if remove_absolute:
            for pattern, reason in ABSOLUTE_CONTRAINDICATIONS.items():
                if pattern in drug_name:
                    return {
                        "reason": reason,
                        "severity": "absolute",
                    }

        # Check relative contraindications
        if remove_relative:
            if drug_name in RELATIVE_CONTRAINDICATIONS:
                ci = RELATIVE_CONTRAINDICATIONS[drug_name]
                return {
                    "reason": ci["reason"],
                    "severity": ci.get("severity", "relative"),
                    "source": ci.get("source", ""),
                }

        return None

    def get_known_failures_report(self) -> List[Dict]:
        """Return list of known pediatric CNS trial failures for annotation."""
        return [
            {"drug_name": d, "reason": "Known pediatric CNS trial failure — deprioritized"}
            for d in KNOWN_PEDIATRIC_CNS_FAILURES
        ]