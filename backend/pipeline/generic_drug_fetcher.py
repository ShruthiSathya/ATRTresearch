"""
generic_drug_fetcher.py  (v3.0 — High-Throughput & Rate Limit Protection)
======================================================
Dynamic generic drug availability checker using FDA and RxNorm APIs.

FIXES FROM v2.0
---------------
1. MASSIVE BATCHING PROTECTION: Unbiased screens can push thousands of drugs
   at once. This script now uses asyncio chunking (batches of 50) and 
   mandatory sleep timers to prevent the FDA API from issuing IP bans (HTTP 429).
2. Expanded KNOWN_GENERIC_STATUS seed to prevent unnecessary network calls
   for the most common repurposed oncology/CNS drugs.
"""

import asyncio
import aiohttp
import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Set

logger = logging.getLogger(__name__)

CACHE_DIR   = Path("/tmp/generic_drug_cache")
CACHE_FILE  = CACHE_DIR / "generic_status_cache.json"
CACHE_TTL_DAYS = 7
MIN_CACHE_ENTRIES = 10  

OPENFDA_DRUGSFDA_URL = "https://api.fda.gov/drug/drugsfda.json"
RXNORM_RXCUI_URL     = "https://rxnav.nlm.nih.gov/REST/rxcui.json"
RXNORM_RELATED_URL   = "https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/allrelated.json"

# Expanded to save API calls during massive screens
KNOWN_GENERIC_STATUS: Dict[str, bool] = {
    "valproic acid": True, "metformin": True, "metformin hcl": True,
    "chloroquine": True, "hydroxychloroquine": True, "sirolimus": True,
    "itraconazole": True, "arsenic trioxide": True, "all-trans retinoic acid": True,
    "tretinoin": True, "temozolomide": True, "dexamethasone": True,
    "lomustine": True, "carmustine": True, "imatinib": True, "bortezomib": True,
    "azacitidine": True, "decitabine": True, "dasatinib": True, "atorvastatin": True,
    "simvastatin": True, "olaparib": False, "venetoclax": False, "vismodegib": False,
    "sonidegib": False, "everolimus": True, "abemaciclib": False, "ribociclib": False,
    "palbociclib": False, "panobinostat": False, "tazemetostat": False,
    "alisertib": False, "birabresib": False, "onc201": False, "paxalisib": False
}

_DRUG_ALIASES: Dict[str, str] = {
    "otx-015": "birabresib", "epz-6438": "tazemetostat", "mln8237": "alisertib",
    "tic-10": "onc201", "gdc-0084": "paxalisib", "rapamycin": "sirolimus",
    "atra": "all-trans retinoic acid", "valproate": "valproic acid",
    "glucophage": "metformin", "plaquenil": "hydroxychloroquine",
    "sporanox": "itraconazole", "temodar": "temozolomide", "velcade": "bortezomib",
}

def _normalise_name(name: str) -> str:
    n = name.lower().strip()
    for suffix in (" hydrochloride", " hcl", " sodium", " phosphate", " sulfate", " mesylate", " acetate", " tartrate"):
        n = n.replace(suffix, "")
    return _DRUG_ALIASES.get(n.strip(), n.strip())

def _seed_lookup(drug_name: str) -> Optional[bool]:
    norm = _normalise_name(drug_name)
    if norm in KNOWN_GENERIC_STATUS: return KNOWN_GENERIC_STATUS[norm]
    raw = drug_name.lower().strip()
    if raw in KNOWN_GENERIC_STATUS: return KNOWN_GENERIC_STATUS[raw]
    return None

class GenericDrugFetcher:
    def __init__(self, use_cache: bool = True):
        self.use_cache = use_cache
        self._cache: Dict[str, Dict] = {}
        self._session: Optional[aiohttp.ClientSession] = None
        self._network_available: Optional[bool] = None  
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        if use_cache: self._load_disk_cache()

    async def check_generic_availability(self, drug_names: List[str]) -> Dict[str, Dict]:
        results: Dict[str, Dict] = {}
        names_needing_live: List[str] = []

        for name in drug_names:
            upper = name.upper().strip()
            if upper in self._cache:
                results[upper] = self._cache[upper]
                continue

            seed_val = _seed_lookup(name)
            if seed_val is not None:
                results[upper] = {
                    "has_generic": seed_val, "source": "SEED",
                    "generic_names": [_normalise_name(name)] if seed_val else [],
                    "note": "FDA Orange Book seed"
                }
                continue

            names_needing_live.append(name)

        if names_needing_live and self._network_available is not False:
            try:
                live_results = await self._batch_query_live(names_needing_live)
                results.update(live_results)
                self._network_available = True
            except Exception as e:
                self._network_available = False
                logger.debug("Live API unavailable: %s. Using fallback.", e)
                for name in names_needing_live:
                    results[name.upper().strip()] = {
                        "has_generic": False, "source": "FALLBACK",
                        "generic_names": [], "note": "Network unavailable"
                    }
        elif names_needing_live:
            for name in names_needing_live:
                results[name.upper().strip()] = {
                    "has_generic": False, "source": "FALLBACK",
                    "generic_names": [], "note": "Network unavailable"
                }

        for name, res in results.items():
            self._cache[name] = res
        self._save_disk_cache()

        return results

    async def _batch_query_live(self, names: List[str]) -> Dict[str, Dict]:
        """Safely chunks requests for massive unbiased screens to prevent IP bans."""
        results: Dict[str, Dict] = {}
        timeout = aiohttp.ClientTimeout(total=15, connect=5)
        connector = aiohttp.TCPConnector(limit=10) # Max 10 concurrent TCP connections
        
        CHUNK_SIZE = 50 
        DELAY_BETWEEN_CHUNKS = 2.0 # Wait 2 seconds between chunks to appease FDA API

        async with aiohttp.ClientSession(timeout=timeout, connector=connector) as session:
            self._session = session
            
            for i in range(0, len(names), CHUNK_SIZE):
                chunk = names[i:i + CHUNK_SIZE]
                if len(names) > 100:
                    logger.info("Processing generic API check chunk %d to %d...", i, min(i+CHUNK_SIZE, len(names)))
                
                tasks = [self._query_single(n) for n in chunk]
                outcomes = await asyncio.gather(*tasks, return_exceptions=True)
                
                for name, outcome in zip(chunk, outcomes):
                    upper = name.upper().strip()
                    if isinstance(outcome, Exception):
                        seed_val = _seed_lookup(name)
                        results[upper] = {
                            "has_generic": seed_val if seed_val is not None else False,
                            "source": "SEED" if seed_val is not None else "FALLBACK",
                            "generic_names": [], "note": f"API error: {outcome}"
                        }
                    else:
                        results[upper] = outcome
                
                if i + CHUNK_SIZE < len(names):
                    await asyncio.sleep(DELAY_BETWEEN_CHUNKS)

            self._session = None
        return results

    async def _query_single(self, drug_name: str) -> Dict:
        seed_val = _seed_lookup(drug_name)
        if seed_val is not None:
            return {"has_generic": seed_val, "source": "SEED", "generic_names": [], "note": "Known status"}

        try:
            fda = await self._query_fda(drug_name)
            if fda is not None: return fda
        except Exception: pass

        try:
            rxn = await self._query_rxnorm(drug_name)
            if rxn is not None: return rxn
        except Exception: pass

        return {"has_generic": False, "source": "FALLBACK", "generic_names": [], "note": f"No data for {drug_name}"}

    async def _query_fda(self, drug_name: str) -> Optional[Dict]:
        if self._session is None: return None
        try:
            params = {"search": f'products.brand_name:"{drug_name}"+OR+openfda.generic_name:"{drug_name}"', "limit": 5}
            async with self._session.get(OPENFDA_DRUGSFDA_URL, params=params) as resp:
                if resp.status != 200: return None
                data = await resp.json()

            results = data.get("results", [])
            has_anda = any(r.get("application_type", "").upper() == "ANDA" for r in results)
            generic_names: Set[str] = set()
            for r in results:
                if r.get("application_type", "").upper() == "ANDA":
                    for p in r.get("products", []):
                        for ing in p.get("active_ingredients", []):
                            n = ing.get("name", "").lower()
                            if n: generic_names.add(n)

            return {
                "has_generic": has_anda,
                "source": "FDA_ANDA" if has_anda else "FDA_NDA",
                "generic_names": sorted(generic_names),
                "note": f"FDA Drugs@FDA: {'ANDA found' if has_anda else 'NDA only'}"
            }
        except (aiohttp.ClientError, asyncio.TimeoutError): return None

    async def _query_rxnorm(self, drug_name: str) -> Optional[Dict]:
        if self._session is None: return None
        try:
            async with self._session.get(RXNORM_RXCUI_URL, params={"name": drug_name, "search": 1}) as resp:
                if resp.status != 200: return None
                data = await resp.json()

            rxcui_list = data.get("idGroup", {}).get("rxnormId", [])
            if not rxcui_list: return None

            rxcui = rxcui_list[0]
            async with self._session.get(RXNORM_RELATED_URL.format(rxcui=rxcui)) as resp:
                if resp.status != 200: return None
                data = await resp.json()

            groups = data.get("allRelatedGroup", {}).get("conceptGroup", [])
            generics = []
            for g in groups:
                if g.get("tty") in ("IN", "PIN", "SCD", "GPCK"):
                    for p in (g.get("conceptProperties") or []):
                        n = p.get("name", "").lower()
                        if n: generics.append(n)

            return {
                "has_generic": len(generics) > 0,
                "source": "RXNORM",
                "generic_names": generics[:5],
                "note": f"RxNorm RxCUI={rxcui}"
            }
        except (aiohttp.ClientError, asyncio.TimeoutError): return None

    def _load_disk_cache(self) -> None:
        if not CACHE_FILE.exists(): return
        try:
            with open(CACHE_FILE) as f: data = json.load(f)
            if len(data) >= MIN_CACHE_ENTRIES: self._cache = data
        except Exception: pass

    def _save_disk_cache(self) -> None:
        if len(self._cache) < 2: return
        try:
            CACHE_DIR.mkdir(parents=True, exist_ok=True)
            with open(CACHE_FILE, "w") as f: json.dump(self._cache, f, indent=2)
        except Exception: pass

_DEFAULT_FETCHER: Optional[GenericDrugFetcher] = None

async def annotate_with_generic_status(candidates: List[Dict]) -> List[Dict]:
    """Annotates massively scaled candidate lists safely."""
    global _DEFAULT_FETCHER
    if _DEFAULT_FETCHER is None:
        _DEFAULT_FETCHER = GenericDrugFetcher(use_cache=True)

    names = [c.get("name") or c.get("drug_name") or "" for c in candidates]

    try:
        status = await _DEFAULT_FETCHER.check_generic_availability(names)
    except Exception as e:
        logger.debug("annotate_with_generic_status failed: %s. Using seed only.", e)
        status = {n.upper().strip(): {"has_generic": _seed_lookup(n) or False, "source": "SEED"} for n in names}

    for c in candidates:
        name = (c.get("name") or c.get("drug_name") or "").upper().strip()
        s = status.get(name, {})
        c["has_generic"] = s.get("has_generic", False)
        c["generic_source"] = s.get("source", "UNKNOWN")
        c["generic_names"] = s.get("generic_names", [])
        c["generic_note"] = s.get("note", "")
        c["cost_category"] = "GENERIC" if s.get("has_generic") else "BRAND" if s.get("source") in ("FDA_NDA", "SEED") else "INVESTIGATIONAL"

    n_generic = sum(1 for c in candidates if c.get("has_generic"))
    logger.info("Screening Annotated: %d/%d candidates have generic formulations available.", n_generic, len(candidates))
    return candidates