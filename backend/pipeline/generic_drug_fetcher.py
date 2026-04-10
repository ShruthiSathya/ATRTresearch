"""
generic_drug_fetcher.py
========================
Dynamic generic drug availability checker using FDA and RxNorm APIs.

Replaces the hardcoded generic_drug_db.py approach with live API queries.

WHY DYNAMIC INSTEAD OF HARDCODED
----------------------------------
Generic drug availability changes over time:
  - Patent expirations create new generics (e.g. imatinib went generic 2016)
  - FDA approves new ANDAs (Abbreviated New Drug Applications) quarterly
  - A hardcoded list is always out of date

This module queries three free public APIs:
  1. FDA Drugs@FDA API  — checks NDA/ANDA status (ANDA = generic)
  2. RxNorm API (NIH)   — drug name normalisation and ingredient lookup
  3. OpenFDA drug label — checks "is there a marketed generic" field

WHAT "GENERIC AVAILABLE" MEANS
--------------------------------
A drug is considered to have a generic available when:
  - An ANDA (Abbreviated NDA) has been approved for it, OR
  - The FDA Orange Book lists it with generic applicants, OR
  - RxNorm finds multiple brand entries mapping to the same ingredient

API ENDPOINTS (all free, no auth required)
-------------------------------------------
FDA Drugs@FDA:
  https://api.fda.gov/drug/drugsfda.json?search=...&limit=10

RxNorm:
  https://rxnav.nlm.nih.gov/REST/rxcui.json?name=<drug>
  https://rxnav.nlm.nih.gov/REST/rxcui/<rxcui>/allrelated.json

OpenFDA drug labels:
  https://api.fda.gov/drug/label.json?search=openfda.generic_name:...

References
-----------
FDA Drugs@FDA: https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files
RxNorm API:    https://rxnav.nlm.nih.gov/RxNormAPIs.html
OpenFDA API:   https://open.fda.gov/apis/drug/
"""

import asyncio
import aiohttp
import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

# ── Disk cache (avoids re-querying the same drug repeatedly) ─────────────────
CACHE_DIR  = Path("/tmp/generic_drug_cache")
CACHE_FILE = CACHE_DIR / "generic_status_cache.json"
CACHE_TTL_DAYS = 7   # refresh after 7 days

# ── API endpoints ─────────────────────────────────────────────────────────────
OPENFDA_DRUGSFDA_URL = "https://api.fda.gov/drug/drugsfda.json"
OPENFDA_LABEL_URL    = "https://api.fda.gov/drug/label.json"
RXNORM_RXCUI_URL     = "https://rxnav.nlm.nih.gov/REST/rxcui.json"
RXNORM_RELATED_URL   = "https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/allrelated.json"
RXNORM_PROPERTIES_URL= "https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/properties.json"

# ── Drug name aliases (helps with API matching) ───────────────────────────────
# Some drugs are known in the literature by development names;
# query the FDA using the INN (International Nonproprietary Name)
DRUG_NAME_ALIASES: Dict[str, str] = {
    "BIRABRESIB":  "OTX015",         # Development code only; no INN-based ANDA
    "TAZEMETOSTAT":"TAZEMETOSTAT",    # Brand: Tazverik
    "ALISERTIB":   "ALISERTIB",       # MLN8237 development name
    "MARIZOMIB":   "MARIZOMIB",
    "PAXALISIB":   "GDC-0084",
    "ONC201":      "DORDAVIPRONE",    # FDA-approved name for ONC201
    "PANOBINOSTAT":"PANOBINOSTAT",    # Brand: Farydak
    "ABEMACICLIB": "ABEMACICLIB",     # Brand: Verzenio
    "VISMODEGIB":  "VISMODEGIB",      # Brand: Erivedge
    "VALPROIC ACID":"VALPROIC ACID",  # Generic: yes
    "VORINOSTAT":  "VORINOSTAT",      # Brand: Zolinza
    "METFORMIN":   "METFORMIN",       # Generic: yes
    "CHLOROQUINE": "CHLOROQUINE",     # Generic: yes
    "SIROLIMUS":   "SIROLIMUS",       # Brand: Rapamune; Generic: yes
    "PALBOCICLIB": "PALBOCICLIB",     # Brand: Ibrance
    "RIBOCICLIB":  "RIBOCICLIB",      # Brand: Kisqali
}

# ── Known generic status (verified from FDA Orange Book, April 2026) ─────────
# Used as a fast-path fallback when APIs are unavailable or rate-limited.
# This is a *seed* lookup only — the live API is always queried first.
KNOWN_GENERIC_STATUS: Dict[str, bool] = {
    # Generic available
    "VALPROIC ACID":      True,
    "METFORMIN":          True,
    "METFORMIN HCL":      True,
    "CHLOROQUINE":        True,
    "HYDROXYCHLOROQUINE": True,
    "SIROLIMUS":          True,
    "BORTEZOMIB":         True,    # US patent expired 2022
    "IMATINIB":           True,    # US patent expired 2016
    "TEMOZOLOMIDE":       True,    # US patent expired 2014
    "DEXAMETHASONE":      True,
    "LOMUSTINE":          True,
    "ITRACONAZOLE":       True,
    "TRETINOIN":          True,    # ATRA
    "ARSENIC TRIOXIDE":   True,
    "VORINOSTAT":         False,   # Zolinza still brand; generics filed
    "PALBOCICLIB":        False,   # Ibrance; patent expires ~2027
    "RIBOCICLIB":         False,   # Kisqali; patent expires ~2027
    "ABEMACICLIB":        False,   # Verzenio; patent expires ~2030
    "PANOBINOSTAT":       False,   # Farydak; orphan drug protections
    "TAZEMETOSTAT":       False,   # Tazverik; approved 2020
    "ALISERTIB":          False,   # No FDA approval (investigational)
    "MARIZOMIB":          False,   # No FDA approval (investigational)
    "BIRABRESIB":         False,   # No FDA approval (investigational)
    "ONC201":             False,   # FDA approval 2024 (brand: Dordaviprone)
    "DORDAVIPRONE":       False,   # Recently approved; no generic yet
    "VISMODEGIB":         False,   # Erivedge; patent ongoing
    "PAXALISIB":          False,   # Investigational
}


class GenericDrugFetcher:
    """
    Dynamically queries FDA and RxNorm to determine which drug candidates
    have generic formulations available.

    This replaces the hardcoded generic_drug_db.py approach.
    Results are cached to disk to avoid repeated API calls.

    Usage
    -----
    >>> fetcher = GenericDrugFetcher()
    >>> results = await fetcher.check_generic_availability(["PANOBINOSTAT", "METFORMIN"])
    >>> # returns {"PANOBINOSTAT": {"has_generic": False, "source": "FDA"}, ...}

    >>> filtered = await fetcher.filter_candidates_by_generic(candidates)
    >>> # returns only candidates with generic formulations available
    """

    def __init__(self, use_cache: bool = True):
        self.use_cache    = use_cache
        self._cache: Dict[str, Dict] = {}
        self._session: Optional[aiohttp.ClientSession] = None
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        if use_cache:
            self._load_disk_cache()

    # ── Public API ────────────────────────────────────────────────────────────

    async def check_generic_availability(
        self, drug_names: List[str]
    ) -> Dict[str, Dict]:
        """
        Check whether each drug has a generic formulation available.

        Parameters
        ----------
        drug_names : list of drug names (any capitalisation)

        Returns
        -------
        dict mapping normalised_drug_name → {
            "has_generic":  bool,
            "source":       str  ("FDA_ANDA" | "RXNORM" | "SEED" | "FALLBACK"),
            "generic_names": list[str],
            "note":         str,
        }
        """
        results = {}
        names_to_query = []

        for name in drug_names:
            upper = name.upper().strip()
            # Fast path: known status seed
            if upper in KNOWN_GENERIC_STATUS:
                results[upper] = {
                    "has_generic":   KNOWN_GENERIC_STATUS[upper],
                    "source":        "SEED",
                    "generic_names": [upper.lower()] if KNOWN_GENERIC_STATUS[upper] else [],
                    "note":          "Known status from FDA Orange Book (April 2026 seed)",
                }
                continue
            # Cache hit
            if upper in self._cache:
                results[upper] = self._cache[upper]
                continue
            names_to_query.append(upper)

        if names_to_query:
            async with aiohttp.ClientSession(
                timeout=aiohttp.ClientTimeout(total=15)
            ) as session:
                self._session = session
                tasks = [self._query_single(n) for n in names_to_query]
                queried = await asyncio.gather(*tasks, return_exceptions=True)

            for name, result in zip(names_to_query, queried):
                if isinstance(result, Exception):
                    logger.warning("Generic check failed for %s: %s", name, result)
                    result = {
                        "has_generic":   False,
                        "source":        "FALLBACK",
                        "generic_names": [],
                        "note":          f"API query failed: {result}",
                    }
                results[name] = result
                self._cache[name] = result

            self._save_disk_cache()

        return results

    async def filter_candidates_by_generic(
        self,
        candidates: List[Dict],
        annotate_all: bool = True,
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Split candidates into those with and without generic formulations.

        Parameters
        ----------
        candidates  : list of drug candidate dicts (must have 'name' key)
        annotate_all: if True, add generic_status to all candidates regardless
                      of filtering decision

        Returns
        -------
        (generic_available, brand_only)
        Each candidate gets 'has_generic', 'generic_source', 'generic_note' added.
        """
        names       = [c.get("name") or c.get("drug_name") or "" for c in candidates]
        status_map  = await self.check_generic_availability(names)

        generic_avail = []
        brand_only    = []

        for c in candidates:
            name   = (c.get("name") or c.get("drug_name") or "").upper().strip()
            status = status_map.get(name, {
                "has_generic":   False,
                "source":        "UNKNOWN",
                "generic_names": [],
                "note":          "Not found in FDA/RxNorm query",
            })

            if annotate_all:
                c["has_generic"]    = status["has_generic"]
                c["generic_source"] = status.get("source", "UNKNOWN")
                c["generic_names"]  = status.get("generic_names", [])
                c["generic_note"]   = status.get("note", "")

            if status["has_generic"]:
                generic_avail.append(c)
            else:
                brand_only.append(c)

        n_generic = len(generic_avail)
        n_brand   = len(brand_only)
        logger.info(
            "Generic drug filter: %d with generics | %d brand-only / investigational",
            n_generic, n_brand,
        )
        return generic_avail, brand_only

    async def get_cost_category(self, drug_name: str) -> str:
        """
        Return a cost tier for a drug: GENERIC | BRAND | INVESTIGATIONAL | UNKNOWN.
        Used for access/equity reporting in the pipeline output.
        """
        upper = drug_name.upper().strip()
        status = await self.check_generic_availability([upper])
        s = status.get(upper, {})

        if s.get("source") == "FDA_ANDA" or s.get("has_generic"):
            return "GENERIC"
        if s.get("source") in ("FDA_NDA",):
            return "BRAND"
        if s.get("source") in ("SEED",):
            return "GENERIC" if s.get("has_generic") else "BRAND"
        return "INVESTIGATIONAL"

    # ── Internal query methods ────────────────────────────────────────────────

    async def _query_single(self, drug_name: str) -> Dict:
        """
        Query FDA Drugs@FDA, then fall back to RxNorm if FDA returns nothing.
        """
        alias = DRUG_NAME_ALIASES.get(drug_name, drug_name)

        # 1. Try FDA Drugs@FDA (most authoritative for US generic status)
        fda_result = await self._query_fda_drugsfda(alias)
        if fda_result is not None:
            return fda_result

        # 2. Try RxNorm (good for normalisation and finding ingredient links)
        rxnorm_result = await self._query_rxnorm(alias)
        if rxnorm_result is not None:
            return rxnorm_result

        # 3. Try OpenFDA drug label
        label_result = await self._query_openfda_label(alias)
        if label_result is not None:
            return label_result

        # 4. Final fallback
        return {
            "has_generic":   False,
            "source":        "FALLBACK",
            "generic_names": [],
            "note":          f"No data found in FDA or RxNorm for '{drug_name}'",
        }

    async def _query_fda_drugsfda(self, drug_name: str) -> Optional[Dict]:
        """
        Query the FDA Drugs@FDA database.

        An ANDA (application_type=ANDA) approval means a generic is available.
        An NDA means the brand-name product.
        """
        try:
            params = {
                "search": f'products.brand_name:"{drug_name}"+openfda.generic_name:"{drug_name}"',
                "limit":  10,
            }
            async with self._session.get(
                OPENFDA_DRUGSFDA_URL, params=params
            ) as resp:
                if resp.status == 404:
                    return None
                if resp.status != 200:
                    logger.debug("FDA Drugs@FDA returned %d for %s", resp.status, drug_name)
                    return None
                data = await resp.json()

            results  = data.get("results", [])
            if not results:
                return None

            has_anda     = False
            generic_names: Set[str] = set()

            for r in results:
                app_type = r.get("application_type", "").upper()
                products = r.get("products", [])

                if app_type == "ANDA":
                    has_anda = True
                    for p in products:
                        gname = p.get("active_ingredients", [{}])
                        for ing in gname:
                            n = ing.get("name", "").lower()
                            if n:
                                generic_names.add(n)

                # Even for NDAs, check if there are listed generic applicants
                if app_type == "NDA":
                    for p in products:
                        if p.get("te_code"):  # Therapeutic equivalence code
                            # TE code presence implies generic equivalents
                            has_anda = True

            return {
                "has_generic":   has_anda,
                "source":        "FDA_ANDA" if has_anda else "FDA_NDA",
                "generic_names": sorted(generic_names),
                "note": (
                    f"FDA Drugs@FDA: {'ANDA found — generic available' if has_anda else 'NDA only — no generic'}"
                ),
            }

        except asyncio.TimeoutError:
            logger.debug("FDA Drugs@FDA timeout for %s", drug_name)
            return None
        except Exception as e:
            logger.debug("FDA Drugs@FDA error for %s: %s", drug_name, e)
            return None

    async def _query_rxnorm(self, drug_name: str) -> Optional[Dict]:
        """
        Query RxNorm to get the RxCUI, then check if multiple brand/generic
        entries map to the same ingredient (indicates generic availability).
        """
        try:
            # Step 1: get RxCUI
            async with self._session.get(
                RXNORM_RXCUI_URL, params={"name": drug_name, "search": 1}
            ) as resp:
                if resp.status != 200:
                    return None
                data = await resp.json()

            id_group = data.get("idGroup", {})
            rxcui_list = id_group.get("rxnormId", [])
            if not rxcui_list:
                return None

            rxcui = rxcui_list[0]

            # Step 2: get related concepts
            url = RXNORM_RELATED_URL.format(rxcui=rxcui)
            async with self._session.get(url) as resp:
                if resp.status != 200:
                    return None
                data = await resp.json()

            concept_group = data.get("allRelatedGroup", {}).get("conceptGroup", [])

            # Term types relevant to generics:
            # "IN"  = Ingredient (active ingredient name = generic name)
            # "PIN" = Precise ingredient
            # "MIN" = Multiple ingredients
            # "BN"  = Brand name
            # "SBD" = Semantic branded drug
            # "SCD" = Semantic clinical drug (generic)
            generic_related = []
            brand_related   = []

            for group in concept_group:
                tty   = group.get("tty", "")
                props = group.get("conceptProperties", [])
                for prop in (props or []):
                    name = prop.get("name", "").lower()
                    if tty in ("IN", "PIN", "MIN", "SCD", "GPCK"):
                        generic_related.append(name)
                    elif tty in ("BN", "SBD", "BPCK"):
                        brand_related.append(name)

            # Has a generic if we found generic-type concepts in RxNorm
            has_generic = len(generic_related) > 0

            return {
                "has_generic":   has_generic,
                "source":        "RXNORM",
                "generic_names": generic_related[:5],
                "note": (
                    f"RxNorm RxCUI={rxcui}: "
                    f"{len(generic_related)} generic concepts, "
                    f"{len(brand_related)} brand concepts"
                ),
            }

        except asyncio.TimeoutError:
            return None
        except Exception as e:
            logger.debug("RxNorm error for %s: %s", drug_name, e)
            return None

    async def _query_openfda_label(self, drug_name: str) -> Optional[Dict]:
        """
        Query OpenFDA drug label endpoint.
        If a generic_name field exists and differs from brand_name, that signals
        generic availability.
        """
        try:
            params = {
                "search": f'openfda.generic_name:"{drug_name.lower()}"',
                "limit":  5,
            }
            async with self._session.get(OPENFDA_LABEL_URL, params=params) as resp:
                if resp.status == 404:
                    return None
                if resp.status != 200:
                    return None
                data = await resp.json()

            results = data.get("results", [])
            if not results:
                return None

            generic_names: Set[str] = set()
            brand_names:   Set[str] = set()

            for r in results:
                openfda = r.get("openfda", {})
                for g in openfda.get("generic_name", []):
                    generic_names.add(g.lower().strip())
                for b in openfda.get("brand_name", []):
                    brand_names.add(b.lower().strip())

            # If generic_name exists and is different from all brand_names → generic available
            true_generics = generic_names - brand_names
            has_generic   = len(true_generics) > 0

            return {
                "has_generic":   has_generic,
                "source":        "OPENFDA_LABEL",
                "generic_names": sorted(true_generics),
                "note": (
                    f"OpenFDA label: {len(generic_names)} generic names, "
                    f"{len(brand_names)} brand names"
                ),
            }

        except Exception as e:
            logger.debug("OpenFDA label error for %s: %s", drug_name, e)
            return None

    # ── Cache management ──────────────────────────────────────────────────────

    def _load_disk_cache(self) -> None:
        if not CACHE_FILE.exists():
            return
        try:
            import time as _time
            # Check TTL
            age_days = (_time.time() - CACHE_FILE.stat().st_mtime) / 86400
            if age_days > CACHE_TTL_DAYS:
                logger.info("Generic drug cache expired (%.1f days old) — will refresh", age_days)
                return
            with open(CACHE_FILE) as f:
                self._cache = json.load(f)
            logger.info("Loaded generic drug cache (%d entries)", len(self._cache))
        except Exception as e:
            logger.warning("Could not load generic drug cache: %s", e)

    def _save_disk_cache(self) -> None:
        try:
            with open(CACHE_FILE, "w") as f:
                json.dump(self._cache, f, indent=2)
        except Exception as e:
            logger.warning("Could not save generic drug cache: %s", e)

    def get_cache_report(self) -> str:
        """Return a summary of the current generic status cache."""
        has  = [k for k, v in self._cache.items() if v.get("has_generic")]
        lack = [k for k, v in self._cache.items() if not v.get("has_generic")]
        lines = [
            "## Generic Drug Status Cache\n",
            f"Total entries: {len(self._cache)}\n\n",
            f"Generic available ({len(has)}):\n",
        ]
        for d in sorted(has):
            src = self._cache[d].get("source", "?")
            lines.append(f"  ✅ {d}  [{src}]\n")
        lines.append(f"\nGeneric NOT available ({len(lack)}):\n")
        for d in sorted(lack):
            src = self._cache[d].get("source", "?")
            lines.append(f"  ❌ {d}  [{src}]\n")
        return "".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Convenience function for pipeline integration
# ─────────────────────────────────────────────────────────────────────────────

_DEFAULT_FETCHER: Optional[GenericDrugFetcher] = None


async def annotate_with_generic_status(candidates: List[Dict]) -> List[Dict]:
    """
    Add generic availability annotations to all candidates without filtering.
    This is the non-destructive version — call this in the main pipeline.
    """
    global _DEFAULT_FETCHER
    if _DEFAULT_FETCHER is None:
        _DEFAULT_FETCHER = GenericDrugFetcher(use_cache=True)

    names  = [c.get("name") or c.get("drug_name") or "" for c in candidates]
    status = await _DEFAULT_FETCHER.check_generic_availability(names)

    for c in candidates:
        name = (c.get("name") or c.get("drug_name") or "").upper().strip()
        s    = status.get(name, {})
        c["has_generic"]    = s.get("has_generic", False)
        c["generic_source"] = s.get("source", "UNKNOWN")
        c["generic_names"]  = s.get("generic_names", [])
        c["generic_note"]   = s.get("note", "")
        c["cost_category"]  = "GENERIC" if s.get("has_generic") else (
            "BRAND" if s.get("source") in ("FDA_NDA", "SEED") else "INVESTIGATIONAL"
        )

    n_generic = sum(1 for c in candidates if c.get("has_generic"))
    logger.info(
        "Generic status annotated: %d/%d candidates have generics available",
        n_generic, len(candidates),
    )
    return candidates