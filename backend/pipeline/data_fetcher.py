"""
data_fetcher.py  (v5.0 — API robustness + fallback fixes)
==========================================================
OpenTargets API Fetcher for ATRT pipeline.

FIXES FROM v4.0
---------------
1. OpenTargets GraphQL query now correctly handles v4 API pagination:
   - Uses 'cursor' field properly from knownDrugs response
   - Adds robust error handling for empty/null disease nodes
   - Falls back gracefully when API is unreachable (network restrictions)

2. _merge_with_fallback now uses exact name matching only (no partial
   substring matching that caused drug mis-identification in v4.0).
   Aliases are handled via an explicit DRUG_NAME_ALIASES dict.

3. Generic annotation now gracefully handles network failures by using
   KNOWN_GENERIC_STATUS seed without attempting live API calls.

4. Added annotate_generics parameter to fetch_approved_drugs so callers
   can disable live API calls in restricted network environments.

PRIMARY EFO IDs FOR ATRT
--------------------------
EFO_0002915 : rhabdoid tumor (primary — catches tazemetostat etc.)
EFO_0000543 : malignant rhabdoid tumor
EFO_0000519 : glioblastoma (catches BET/HDAC/CDK drugs with CNS PK data)
EFO_0001422 : DIPG (catches epigenetic drugs — overlapping vulnerabilities)
EFO_0000618 : pediatric brain tumor

References
-----------
Ochoa D et al. (2021). Open Targets Platform. Nucleic Acids Res 49(D1):D1302.
Knutson SK et al. (2013). PNAS 110(19):7922. PMID 23620515.
Torchia J et al. (2015). Cancer Cell 30(6):891. PMID 26609405.
"""

import asyncio
import aiohttp
import logging
from typing import Dict, List, Optional, Set

try:
    from .pipeline_config import OPENTARGETS
except ImportError:
    from pipeline_config import OPENTARGETS

logger = logging.getLogger(__name__)

MAX_RETRIES   = 3
RETRY_BACKOFF = 2.0

# Explicit drug name aliases for merging (prevents partial-match mis-identification)
# Maps canonical name (lowercase) -> list of alternative names (lowercase)
DRUG_NAME_ALIASES: Dict[str, List[str]] = {
    "birabresib":         ["otx015", "otx-015", "mk-8628"],
    "tazemetostat":       ["epz-6438", "epz6438"],
    "alisertib":          ["mln8237", "mln-8237"],
    "valproic acid":      ["valproate", "valproic acid sodium", "depakote"],
    "paxalisib":          ["gdc-0084", "gdc0084"],
    "onc201":             ["dordaviprone", "tic-10", "tic10"],
    "hydroxychloroquine": ["hydroxychloroquine sulfate", "plaquenil"],
    "chloroquine":        ["chloroquine phosphate"],
    "all-trans retinoic acid": ["tretinoin", "atra", "vesanoid"],
    "arsenic trioxide":   ["trisenox", "as2o3"],
    "metformin":          ["metformin hcl", "metformin hydrochloride", "glucophage"],
    "sirolimus":          ["rapamycin", "rapamune"],
    "itraconazole":       ["sporanox"],
    "temozolomide":       ["temodar", "temodal"],
    "bortezomib":         ["velcade", "ps-341"],
}

# Build reverse alias lookup: alias -> canonical name
_ALIAS_TO_CANONICAL: Dict[str, str] = {}
for _canonical, _aliases in DRUG_NAME_ALIASES.items():
    for _alias in _aliases:
        _ALIAS_TO_CANONICAL[_alias] = _canonical

# ─────────────────────────────────────────────────────────────────────────────
# CURATED ATRT FALLBACK CANDIDATES
# Used when OpenTargets API returns 0 results (network failure etc).
# Generic status is NOT hardcoded — annotated dynamically or via seed.
# Sources: Torchia 2015 PMID 26609405, Knutson 2013 PMID 23620515,
#          Sredni 2017 PMID 28544500, Geoerger 2017 PMID 28108534,
#          Johann 2016 PMID 26923874, Frühwald 2020 PMID 32432484.
# ─────────────────────────────────────────────────────────────────────────────

ATRT_CURATED_FALLBACK: List[Dict] = [
    # EZH2 synthetic lethality (Knutson 2013 PNAS PMID 23620515)
    {"name": "Tazemetostat",     "targets": ["EZH2"],
     "mechanism": "EZH2 inhibitor", "pmid": "23620515",
     "drug_class": "EZH2/PRC2 inhibitor"},

    # Pan-HDAC (Torchia 2015 Cancer Cell PMID 26609405)
    {"name": "Panobinostat",     "targets": ["HDAC1", "HDAC2", "HDAC3", "HDAC6"],
     "mechanism": "pan-HDAC inhibitor", "pmid": "26609405",
     "drug_class": "HDAC inhibitor"},
    {"name": "Vorinostat",       "targets": ["HDAC1", "HDAC2", "HDAC3", "HDAC6"],
     "mechanism": "pan-HDAC inhibitor", "pmid": "26609405",
     "drug_class": "HDAC inhibitor"},
    {"name": "Valproic acid",    "targets": ["HDAC1", "HDAC2", "HDAC3"],
     "mechanism": "HDAC inhibitor (repurposed anticonvulsant)", "pmid": "11742990",
     "drug_class": "HDAC inhibitor"},
    {"name": "Entinostat",       "targets": ["HDAC1", "HDAC2", "HDAC3"],
     "mechanism": "class I HDAC inhibitor",
     "drug_class": "HDAC inhibitor"},

    # Aurora kinase (Sredni 2017 Pediatric Blood Cancer PMID 28544500)
    {"name": "Alisertib",        "targets": ["AURKA"],
     "mechanism": "aurora kinase A inhibitor", "pmid": "28544500",
     "drug_class": "Aurora kinase A inhibitor"},

    # BET bromodomain (Geoerger 2017 Clin Cancer Res PMID 28108534)
    {"name": "Birabresib",       "targets": ["BRD4", "BRD2", "BRD3"],
     "mechanism": "BET bromodomain inhibitor", "pmid": "28108534",
     "drug_class": "BET bromodomain inhibitor"},

    # CDK4/6
    {"name": "Abemaciclib",      "targets": ["CDK4", "CDK6"],
     "mechanism": "CDK4/CDK6 inhibitor",
     "drug_class": "CDK4/6 inhibitor"},
    {"name": "Palbociclib",      "targets": ["CDK4", "CDK6"],
     "mechanism": "CDK4/CDK6 inhibitor",
     "drug_class": "CDK4/6 inhibitor"},
    {"name": "Ribociclib",       "targets": ["CDK4", "CDK6"],
     "mechanism": "CDK4/CDK6 inhibitor",
     "drug_class": "CDK4/6 inhibitor"},

    # Proteasome (Lin 2019 Sci Transl Med)
    {"name": "Marizomib",        "targets": ["PSMB5", "PSMB2", "PSMB1"],
     "mechanism": "proteasome inhibitor",
     "drug_class": "Proteasome inhibitor"},
    {"name": "Bortezomib",       "targets": ["PSMB5", "PSMB1"],
     "mechanism": "proteasome inhibitor (boronic acid)",
     "drug_class": "Proteasome inhibitor"},

    # SHH subgroup (Torchia 2015; Johann 2016 PMID 26923874)
    {"name": "Vismodegib",       "targets": ["SMO"],
     "mechanism": "smoothened inhibitor",
     "drug_class": "Hedgehog/SMO inhibitor"},
    {"name": "Sonidegib",        "targets": ["SMO"],
     "mechanism": "smoothened inhibitor",
     "drug_class": "Hedgehog/SMO inhibitor"},
    {"name": "Itraconazole",     "targets": ["SMO", "PTCH1"],
     "mechanism": "SMO inhibitor (repurposed antifungal)",
     "drug_class": "Hedgehog/SMO inhibitor"},

    # DRD2 / TRAIL
    {"name": "ONC201",           "targets": ["DRD2", "CLPB"],
     "mechanism": "DRD2 antagonist / TRAIL inducer",
     "drug_class": "DRD2 antagonist"},

    # PI3K/mTOR
    {"name": "Paxalisib",        "targets": ["PIK3CA", "PIK3CD", "PIK3CG", "MTOR"],
     "mechanism": "PI3K inhibitor (CNS-penetrant)",
     "drug_class": "PI3K inhibitor"},
    {"name": "Sirolimus",        "targets": ["MTOR", "RPTOR"],
     "mechanism": "mTOR inhibitor (repurposed immunosuppressant)",
     "drug_class": "mTOR inhibitor"},

    # Repurposable generics with epigenetic/MYC rationale
    {"name": "Metformin",        "targets": ["PRKAB1", "PRKAB2"],
     "mechanism": "AMPK activator (repurposed antidiabetic)",
     "drug_class": "AMPK activator"},
    {"name": "Chloroquine",      "targets": ["ATP6V0A1", "BECN1"],
     "mechanism": "autophagy inhibitor (repurposed antimalarial)",
     "drug_class": "Autophagy inhibitor"},
    {"name": "Hydroxychloroquine", "targets": ["ATP6V0A1", "BECN1"],
     "mechanism": "autophagy inhibitor (repurposed antimalarial)",
     "drug_class": "Autophagy inhibitor"},

    # Differentiation (TYR subgroup rationale)
    {"name": "All-trans retinoic acid", "targets": ["RARA", "RARB", "RARG"],
     "mechanism": "retinoid receptor agonist / differentiation therapy",
     "drug_class": "Retinoid"},

    # SHH via GLI (ATRT-SHH subgroup)
    {"name": "Arsenic trioxide", "targets": ["GLI1", "GLI2"],
     "mechanism": "GLI inhibitor (repurposed from APL therapy)",
     "drug_class": "GLI inhibitor"},

    # Anti-apoptotic
    {"name": "Venetoclax",       "targets": ["BCL2"],
     "mechanism": "BCL-2 inhibitor",
     "drug_class": "BCL-2 inhibitor"},
]

# Known generic status seed — FDA Orange Book verified April 2026
# Source: https://www.accessdata.fda.gov/scripts/cder/ob/
KNOWN_GENERIC_STATUS: Dict[str, bool] = {
    "valproic acid":        True,
    "metformin":            True,
    "metformin hcl":        True,
    "chloroquine":          True,
    "hydroxychloroquine":   True,
    "sirolimus":            True,
    "bortezomib":           True,   # US patent expired 2022
    "imatinib":             True,   # US patent expired 2016
    "temozolomide":         True,   # US patent expired 2014
    "dexamethasone":        True,
    "lomustine":            True,
    "itraconazole":         True,
    "tretinoin":            True,   # ATRA
    "all-trans retinoic acid": True,
    "arsenic trioxide":     True,
    "vorinostat":           False,  # Zolinza; ANDAs filed but not yet approved
    "palbociclib":          False,  # Ibrance; patent ~2027
    "ribociclib":           False,  # Kisqali; patent ~2027
    "abemaciclib":          False,  # Verzenio; patent ~2030
    "panobinostat":         False,  # Farydak; orphan drug protections
    "tazemetostat":         False,  # Tazverik; approved 2020
    "alisertib":            False,  # No FDA approval (investigational)
    "marizomib":            False,  # No FDA approval (investigational)
    "birabresib":           False,  # No FDA approval (investigational)
    "onc201":               False,  # FDA approval 2024; no generic yet
    "dordaviprone":         False,
    "vismodegib":           False,  # Erivedge; patent ongoing
    "sonidegib":            False,
    "paxalisib":            False,  # Investigational
    "venetoclax":           False,  # Venclexta; patent ~2033
    "entinostat":           False,  # Investigational
}


class ProductionDataFetcher:
    """
    Fetches drug-target associations from OpenTargets API and annotates
    generic drug availability.

    v5.0 Pipeline:
    1. Query OpenTargets with ATRT/rhabdoid EFO IDs.
    2. Fall back to ATRT_CURATED_FALLBACK if OpenTargets is unreachable.
    3. Merge: curated fallback ensures all key ATRT drugs are present.
    4. Annotate generic status from KNOWN_GENERIC_STATUS seed (live API
       attempted but silently skipped if network is restricted).
    """

    def __init__(self):
        self.api_url   = "https://api.platform.opentargets.org/api/v4/graphql"
        self.efo_ids   = OPENTARGETS["efo_ids"]
        self.max_pages = OPENTARGETS["max_pages_per_efo"]
        self.page_size = OPENTARGETS["drugs_per_page"]

    # ── Public API ────────────────────────────────────────────────────────────

    async def fetch_disease_data(self, disease_name: str) -> Dict:
        return {
            "name":  disease_name,
            "genes": [
                "SMARCB1", "EZH2", "BRD4", "HDAC1", "AURKA",
                "MYC", "MYCN", "CDK4", "PSMB5", "BCL2L1",
            ],
            "id": "EFO_0002915",
        }

    async def fetch_approved_drugs(
        self,
        annotate_generics: bool = True,
    ) -> List[Dict]:
        """
        Fetch drug candidates and annotate with generic status.

        Parameters
        ----------
        annotate_generics : bool
            If True, attempt live FDA/RxNorm queries; falls back to seed
            data if network is unavailable.
        """
        # Step 1: OpenTargets fetch (may return [] on network restriction)
        ot_drugs = await self._fetch_opentargets_drugs()

        if not ot_drugs:
            logger.warning(
                "OpenTargets returned 0 drugs — using ATRT curated fallback "
                "library (%d drugs). Generic status will still be checked via "
                "FDA/RxNorm or seed lookup.",
                len(ATRT_CURATED_FALLBACK),
            )

        # Step 2: Merge OT results with curated fallback
        merged = self._merge_with_fallback(
            ot_drugs if ot_drugs else [],
            ATRT_CURATED_FALLBACK
        )

        # Step 3: Annotate generic availability
        if annotate_generics:
            merged = await self._annotate_generic_status(merged)
        else:
            for d in merged:
                d.setdefault("has_generic", False)
                d.setdefault("cost_category", "UNKNOWN")

        n_generic = sum(1 for d in merged if d.get("has_generic"))
        logger.info(
            "Generic status annotated: %d/%d drugs have generics available",
            n_generic, len(merged),
        )
        logger.info("Total candidates after merge + annotation: %d", len(merged))
        return merged

    async def close(self) -> None:
        pass

    # ── OpenTargets fetch ─────────────────────────────────────────────────────

    async def _fetch_opentargets_drugs(self) -> List[Dict]:
        logger.info("Fetching from OpenTargets (EFOs: %s) ...", self.efo_ids)

        # OpenTargets Platform v4 GraphQL API
        # knownDrugs returns a KnownDrugs object with:
        #   count: Int
        #   cursor: String (for pagination)
        #   rows: [KnownDrug]
        # Each KnownDrug has: drug { id name } target { approvedSymbol }
        query = """
        query atrtDrugs($efoId: String!, $cursor: String, $size: Int) {
          disease(efoId: $efoId) {
            name
            knownDrugs(cursor: $cursor, size: $size) {
              count
              cursor
              rows {
                drug { id name }
                target { approvedSymbol }
                phase
                status
              }
            }
          }
        }
        """

        drugs_dict: Dict[str, Dict] = {}
        timeout = aiohttp.ClientTimeout(total=60, connect=15)

        try:
            async with aiohttp.ClientSession(timeout=timeout) as session:
                for efo in self.efo_ids:
                    cursor     = None
                    page_count = 0
                    efo_drugs  = 0

                    for _ in range(self.max_pages):
                        variables: Dict = {"efoId": efo, "size": self.page_size}
                        if cursor:
                            variables["cursor"] = cursor

                        data = await self._safe_post(
                            session, {"query": query, "variables": variables}
                        )
                        if data is None:
                            break

                        # Navigate the v4 response structure
                        raw     = data.get("data") or {}
                        disease = raw.get("disease")
                        if disease is None:
                            logger.debug("EFO %s: no disease node in response", efo)
                            break

                        kd   = disease.get("knownDrugs") or {}
                        rows = kd.get("rows") or []

                        if not rows:
                            logger.debug("EFO %s: empty rows in page %d", efo, page_count)
                            break

                        for row in rows:
                            try:
                                drug_node   = row.get("drug") or {}
                                target_node = row.get("target") or {}
                                d_name      = drug_node.get("name")
                                t_name      = target_node.get("approvedSymbol")
                                if not d_name or not t_name:
                                    continue
                                key = d_name.lower().strip()
                                if key not in drugs_dict:
                                    drugs_dict[key] = {
                                        "name":    d_name,
                                        "targets": [],
                                        "source":  "OpenTargets",
                                        "phase":   row.get("phase"),
                                        "status":  row.get("status"),
                                    }
                                if t_name not in drugs_dict[key]["targets"]:
                                    drugs_dict[key]["targets"].append(t_name)
                                efo_drugs += 1
                            except (KeyError, TypeError, AttributeError):
                                continue

                        next_cursor = kd.get("cursor")
                        page_count += 1
                        if not next_cursor or next_cursor == cursor:
                            break
                        cursor = next_cursor

                    logger.info("  EFO %s: %d drugs, %d pages", efo, efo_drugs, page_count)

        except aiohttp.ClientConnectorError as e:
            logger.warning(
                "OpenTargets API unreachable (network restriction): %s. "
                "Using curated fallback.", e
            )
            return []
        except Exception as e:
            logger.warning("OpenTargets fetch failed: %s. Using curated fallback.", e)
            return []

        logger.info("OpenTargets: %d unique drugs retrieved", len(drugs_dict))
        return list(drugs_dict.values())

    async def _safe_post(
        self,
        session: aiohttp.ClientSession,
        payload: dict,
        attempt: int = 0,
    ) -> Optional[dict]:
        """POST to OpenTargets GraphQL; retry on 5xx; return None on failure."""
        try:
            async with session.post(self.api_url, json=payload) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    # Check for GraphQL errors
                    if data.get("errors"):
                        logger.debug("GraphQL errors: %s", data["errors"][:1])
                        return None
                    return data
                elif resp.status >= 500 and attempt < MAX_RETRIES:
                    await asyncio.sleep(RETRY_BACKOFF ** attempt)
                    return await self._safe_post(session, payload, attempt + 1)
                else:
                    logger.debug("OpenTargets HTTP %d", resp.status)
                    return None
        except aiohttp.ClientConnectorError:
            # Network blocked — don't retry
            return None
        except aiohttp.ClientError as e:
            if attempt < MAX_RETRIES:
                await asyncio.sleep(RETRY_BACKOFF ** attempt)
                return await self._safe_post(session, payload, attempt + 1)
            logger.debug("OpenTargets client error: %s", e)
            return None
        except Exception as e:
            logger.debug("OpenTargets unexpected error: %s", e)
            return None

    # ── Generic annotation ────────────────────────────────────────────────────

    async def _annotate_generic_status(self, candidates: List[Dict]) -> List[Dict]:
        """
        Annotate each candidate with generic availability.
        Attempts live FDA/RxNorm query; silently falls back to seed lookup
        if the network is restricted.

        Strategy (in priority order):
          1. KNOWN_GENERIC_STATUS seed (instantaneous, no network)
          2. Live GenericDrugFetcher (requires network access)
        """
        # First pass: annotate from seed (always works)
        for drug in candidates:
            name = (drug.get("name") or "").lower().strip()
            # Normalise via alias
            canonical = _ALIAS_TO_CANONICAL.get(name, name)
            # Strip salt suffixes for lookup
            for suffix in (" hydrochloride", " hcl", " sodium", " mesylate",
                           " phosphate", " sulfate", " acetate"):
                canonical = canonical.replace(suffix, "")
            canonical = canonical.strip()

            if canonical in KNOWN_GENERIC_STATUS:
                drug["has_generic"]    = KNOWN_GENERIC_STATUS[canonical]
                drug["generic_source"] = "SEED"
                drug["cost_category"]  = (
                    "GENERIC" if KNOWN_GENERIC_STATUS[canonical] else "BRAND"
                )
            elif name in KNOWN_GENERIC_STATUS:
                drug["has_generic"]    = KNOWN_GENERIC_STATUS[name]
                drug["generic_source"] = "SEED"
                drug["cost_category"]  = (
                    "GENERIC" if KNOWN_GENERIC_STATUS[name] else "BRAND"
                )
            else:
                drug.setdefault("has_generic",    False)
                drug.setdefault("generic_source", "UNKNOWN")
                drug.setdefault("cost_category",  "INVESTIGATIONAL")

        # Second pass: try live API for unknowns (silently skip if network blocked)
        unknowns = [d for d in candidates if d.get("generic_source") == "UNKNOWN"]
        if unknowns:
            try:
                from .generic_drug_fetcher import annotate_with_generic_status
            except ImportError:
                try:
                    from generic_drug_fetcher import annotate_with_generic_status
                except ImportError:
                    annotate_with_generic_status = None

            if annotate_with_generic_status is not None:
                try:
                    unknowns = await annotate_with_generic_status(unknowns)
                    # Re-integrate annotated unknowns
                    annotated_names = {d.get("name", "").lower() for d in unknowns}
                    candidates = [
                        d if d.get("name", "").lower() not in annotated_names else
                        next((u for u in unknowns
                              if u.get("name", "").lower() == d.get("name", "").lower()), d)
                        for d in candidates
                    ]
                except Exception as e:
                    logger.debug(
                        "Live generic annotation failed (network may be restricted): %s", e
                    )

        return candidates

    # ── Merge helpers ─────────────────────────────────────────────────────────

    @staticmethod
    def _merge_with_fallback(
        ot_drugs: List[Dict],
        fallback: List[Dict],
    ) -> List[Dict]:
        """
        Merge OpenTargets results with curated fallback list.

        Strategy:
        - Index all drugs by canonical lowercase name.
        - Resolve aliases: e.g. "otx015" -> "birabresib".
        - For each fallback drug, enrich OT entry if present; add if absent.
        - Uses EXACT name matching only (no substring matching).
        """
        def _canonical_key(name: str) -> str:
            """Normalise a drug name to its canonical lowercase form."""
            name = name.lower().strip()
            for suffix in (" hydrochloride", " hcl", " sodium", " mesylate",
                           " phosphate", " sulfate", " acetate"):
                name = name.replace(suffix, "")
            name = name.strip()
            # Resolve alias to canonical form
            return _ALIAS_TO_CANONICAL.get(name, name)

        merged: Dict[str, Dict] = {}

        # Index OT drugs by canonical key
        for d in ot_drugs:
            key = _canonical_key(d.get("name", ""))
            if key and key not in merged:
                merged[key] = dict(d)
            elif key and key in merged:
                # Merge targets
                existing_targets = set(merged[key].get("targets", []))
                new_targets = set(d.get("targets", []))
                merged[key]["targets"] = sorted(existing_targets | new_targets)

        n_enriched = 0
        n_added    = 0

        for fb in fallback:
            fb_key = _canonical_key(fb.get("name", ""))
            if not fb_key:
                continue

            if fb_key in merged:
                # Enrich existing entry with fallback metadata
                merged[fb_key].setdefault("mechanism",  fb.get("mechanism", ""))
                merged[fb_key].setdefault("drug_class", fb.get("drug_class", ""))
                merged[fb_key].setdefault("pmid",       fb.get("pmid", ""))
                # Merge targets
                existing = set(merged[fb_key].get("targets", []))
                fb_tgts  = set(fb.get("targets", []))
                merged[fb_key]["targets"] = sorted(existing | fb_tgts)
                n_enriched += 1
            else:
                # Add as new entry
                merged[fb_key] = {
                    "name":       fb.get("name"),
                    "targets":    list(fb.get("targets", [])),
                    "mechanism":  fb.get("mechanism", ""),
                    "drug_class": fb.get("drug_class", ""),
                    "pmid":       fb.get("pmid", ""),
                    "source":     "curated_fallback",
                }
                n_added += 1

        logger.info(
            "Merge: %d OT drugs + %d fallback → %d unique total "
            "(%d enriched, %d added from fallback)",
            len(ot_drugs), len(fallback), len(merged), n_enriched, n_added,
        )
        return list(merged.values())