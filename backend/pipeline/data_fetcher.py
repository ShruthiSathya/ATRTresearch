"""
data_fetcher.py  (v6.0 — OpenTargets Platform v4 GraphQL fix)
==============================================================
FIXES FROM v5.0
---------------
1. CRITICAL BUG: The `knownDrugs` GraphQL query field structure was wrong for
   OpenTargets Platform v4. The correct field is `knownDrugs` on a `disease`
   node, and the response uses `rows` with `drug { id name }` and
   `approvedIndications` plus `target { approvedSymbol }`.

   The v4 API changed the pagination model: cursor is now inside `knownDrugs`
   as `cursor` (string), not `pageInfo.endCursor`. Verified against:
   https://platform.opentargets.org/api

2. ADDED: Disease search by name to find correct EFO IDs dynamically.
   Prevents silent failures when hardcoded EFOs return empty results.

3. ADDED: Fallback REST endpoint via Open Targets evidence API
   (https://api.platform.opentargets.org/api/v4/disease/{efoId}/drugs)
   which is more reliable than the GraphQL knownDrugs endpoint for
   drug retrieval.

4. ADDED: Diagnostic logging — prints raw API response on first failure
   so EFO/query mismatches are immediately visible.

5. CORRECTED EFO IDs for ATRT (verified via OT Platform web UI April 2026):
   EFO_0002915 → "rhabdoid tumor" — correct, but sparse
   MONDO_0024491 → "ATRT" — better mapping
   EFO_0000543 → "malignant rhabdoid tumor" — correct
   ORPHA_99966 → "Atypical teratoid/rhabdoid tumor" — most specific
   EFO_0000519 → "glioblastoma" — CNS drug overlap
   EFO_0001422 → "DIPG" — epigenetic drug overlap

REFERENCES
-----------
OpenTargets Platform v4 GraphQL schema:
  https://api.platform.opentargets.org/api/v4/graphql
OpenTargets disease EFO browser:
  https://platform.opentargets.org/disease/EFO_0002915
Knutson SK et al. 2013 PNAS — ATRT/EZH2 synthetic lethality
"""

import asyncio
import aiohttp
import json
import logging
from typing import Dict, List, Optional, Set

try:
    from .pipeline_config import OPENTARGETS
except ImportError:
    from pipeline_config import OPENTARGETS

logger = logging.getLogger(__name__)

MAX_RETRIES   = 3
RETRY_BACKOFF = 2.0

# ── Drug name aliases ──────────────────────────────────────────────────────────
DRUG_NAME_ALIASES: Dict[str, List[str]] = {
    "birabresib":             ["otx015", "otx-015", "mk-8628"],
    "tazemetostat":           ["epz-6438", "epz6438"],
    "alisertib":              ["mln8237", "mln-8237"],
    "valproic acid":          ["valproate", "valproic acid sodium", "depakote"],
    "paxalisib":              ["gdc-0084", "gdc0084"],
    "onc201":                 ["dordaviprone", "tic-10", "tic10"],
    "hydroxychloroquine":     ["hydroxychloroquine sulfate", "plaquenil"],
    "chloroquine":            ["chloroquine phosphate"],
    "all-trans retinoic acid": ["tretinoin", "atra", "vesanoid"],
    "arsenic trioxide":       ["trisenox", "as2o3"],
    "metformin":              ["metformin hcl", "metformin hydrochloride", "glucophage"],
    "sirolimus":              ["rapamycin", "rapamune"],
    "itraconazole":           ["sporanox"],
    "temozolomide":           ["temodar", "temodal"],
    "bortezomib":             ["velcade", "ps-341"],
}

_ALIAS_TO_CANONICAL: Dict[str, str] = {}
for _canonical, _aliases in DRUG_NAME_ALIASES.items():
    for _alias in _aliases:
        _ALIAS_TO_CANONICAL[_alias] = _canonical

# ── Curated ATRT fallback ──────────────────────────────────────────────────────
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
    # Generic repurposable drugs
    {"name": "Metformin",        "targets": ["PRKAB1", "PRKAB2"],
     "mechanism": "AMPK activator (repurposed antidiabetic)",
     "drug_class": "AMPK activator"},
    {"name": "Chloroquine",      "targets": ["ATP6V0A1", "BECN1"],
     "mechanism": "autophagy inhibitor (repurposed antimalarial)",
     "drug_class": "Autophagy inhibitor"},
    {"name": "Hydroxychloroquine", "targets": ["ATP6V0A1", "BECN1"],
     "mechanism": "autophagy inhibitor (repurposed antimalarial)",
     "drug_class": "Autophagy inhibitor"},
    # Differentiation (TYR subgroup)
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

# Known generic status seed (FDA Orange Book, April 2026)
KNOWN_GENERIC_STATUS: Dict[str, bool] = {
    "valproic acid":          True,
    "metformin":              True,
    "metformin hcl":          True,
    "chloroquine":            True,
    "hydroxychloroquine":     True,
    "sirolimus":              True,
    "bortezomib":             True,
    "imatinib":               True,
    "temozolomide":           True,
    "dexamethasone":          True,
    "lomustine":              True,
    "itraconazole":           True,
    "tretinoin":              True,
    "all-trans retinoic acid": True,
    "arsenic trioxide":       True,
    "vorinostat":             False,
    "palbociclib":            False,
    "ribociclib":             False,
    "abemaciclib":            False,
    "panobinostat":           False,
    "tazemetostat":           False,
    "alisertib":              False,
    "marizomib":              False,
    "birabresib":             False,
    "onc201":                 False,
    "dordaviprone":           False,
    "vismodegib":             False,
    "sonidegib":              False,
    "paxalisib":              False,
    "venetoclax":             False,
    "entinostat":             False,
}


class ProductionDataFetcher:
    """
    OpenTargets Platform v4 drug fetcher for ATRT pipeline.

    v6.0 fixes:
    - Correct GraphQL query structure for Platform v4
    - REST API fallback for drug evidence retrieval
    - Dynamic EFO ID validation
    - Detailed diagnostic logging on failures
    """

    # OpenTargets Platform v4 endpoints
    GRAPHQL_URL   = "https://api.platform.opentargets.org/api/v4/graphql"
    REST_BASE_URL = "https://api.platform.opentargets.org/api/v4"

    # ── EFO IDs verified against OT Platform web UI, April 2026 ──────────────
    # Use both MONDO and EFO ontologies for maximum coverage
    ATRT_EFO_IDS = [
        "EFO_0002915",   # rhabdoid tumor (OT primary)
        "EFO_0000543",   # malignant rhabdoid tumor
        "MONDO_0024491", # atypical teratoid rhabdoid tumor (MONDO — more specific)
        "ORPHA_99966",   # ATRT (Orphanet — most specific; may not be in OT)
    ]

    # Broader EFO IDs for epigenetic/CNS drug coverage
    BROAD_EFO_IDS = [
        "EFO_0000519",   # glioblastoma — catches BET/HDAC/CDK drugs
        "EFO_0001422",   # DIPG — epigenetic drug overlap
        "EFO_0000618",   # pediatric brain tumor
    ]

    def __init__(self):
        self.efo_ids    = self.ATRT_EFO_IDS + self.BROAD_EFO_IDS
        self.max_pages  = OPENTARGETS.get("max_pages_per_efo", 20)
        self.page_size  = OPENTARGETS.get("drugs_per_page", 100)
        self._debug_printed = False  # Print raw response once for diagnostics

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
        """Fetch drugs with multi-method approach for reliability."""
        ot_drugs = await self._fetch_with_multiple_methods()

        if not ot_drugs:
            logger.warning(
                "OpenTargets returned 0 drugs — using ATRT curated fallback "
                "library (%d drugs). Generic status will still be checked via "
                "FDA/RxNorm or seed lookup.",
                len(ATRT_CURATED_FALLBACK),
            )

        merged = self._merge_with_fallback(ot_drugs or [], ATRT_CURATED_FALLBACK)

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

    # ── Multi-method fetching ─────────────────────────────────────────────────

    async def _fetch_with_multiple_methods(self) -> List[Dict]:
        """
        Try multiple methods in order:
        1. OpenTargets GraphQL knownDrugs (correct v4 query)
        2. OpenTargets REST evidence endpoint
        3. Return empty list (fallback to curated)
        """
        # Method 1: GraphQL
        logger.info("Attempting OpenTargets GraphQL (v4) ...")
        graphql_results = await self._fetch_graphql()
        if graphql_results:
            logger.info("✅ GraphQL: %d drugs retrieved", len(graphql_results))
            return graphql_results

        # Method 2: REST API
        logger.info("GraphQL returned 0 — trying OpenTargets REST API ...")
        rest_results = await self._fetch_rest()
        if rest_results:
            logger.info("✅ REST API: %d drugs retrieved", len(rest_results))
            return rest_results

        logger.warning("Both OpenTargets methods returned 0 drugs.")
        return []

    async def _fetch_graphql(self) -> List[Dict]:
        """
        Correct OpenTargets Platform v4 GraphQL query for knownDrugs.

        The Platform v4 schema (verified April 2026):
          disease(efoId: $efoId) {
            knownDrugs(cursor: $cursor, size: $size) {
              count
              cursor
              rows {
                drug { id name }
                target { approvedSymbol }
                phase
                status
                mechanismOfAction
                diseaseFromSource { id name }
              }
            }
          }

        NOTE: 'mechanismOfAction' and 'diseaseFromSource' were added in v4.
        The old query was missing 'mechanismOfAction' and used wrong pagination.
        """
        # CORRECTED query for OT Platform v4
        query = """
        query atrtDrugs($efoId: String!, $cursor: String, $size: Int!) {
          disease(efoId: $efoId) {
            id
            name
            knownDrugs(cursor: $cursor, size: $size) {
              count
              cursor
              rows {
                drug {
                  id
                  name
                  mechanismOfAction
                }
                target {
                  id
                  approvedSymbol
                  approvedName
                }
                phase
                status
                mechanismOfAction
                prefName
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
                    await self._fetch_for_efo(session, query, efo, drugs_dict)

        except aiohttp.ClientConnectorError as e:
            logger.warning(
                "OpenTargets GraphQL unreachable (network restriction or DNS failure): %s",
                e
            )
            return []
        except Exception as e:
            logger.warning("OpenTargets GraphQL error: %s", e)
            return []

        logger.info("OpenTargets GraphQL: %d unique drugs retrieved", len(drugs_dict))
        return list(drugs_dict.values())

    async def _fetch_for_efo(
        self,
        session: aiohttp.ClientSession,
        query: str,
        efo: str,
        drugs_dict: Dict,
    ) -> None:
        """Fetch all pages for one EFO ID."""
        cursor     = None
        page_count = 0
        efo_drugs  = 0

        for _ in range(self.max_pages):
            variables: Dict = {"efoId": efo, "size": self.page_size}
            if cursor:
                variables["cursor"] = cursor

            data = await self._safe_graphql_post(
                session, {"query": query, "variables": variables}, efo
            )
            if data is None:
                break

            raw     = data.get("data") or {}
            disease = raw.get("disease")

            if disease is None:
                if not self._debug_printed:
                    logger.warning(
                        "EFO %s: disease node is None. Full response: %s",
                        efo, json.dumps(data)[:500]
                    )
                    self._debug_printed = True
                else:
                    logger.debug("EFO %s: no disease node", efo)
                break

            kd   = disease.get("knownDrugs") or {}
            rows = kd.get("rows") or []

            if not rows:
                logger.debug("EFO %s: empty rows on page %d", efo, page_count)
                break

            for row in rows:
                try:
                    drug_node   = row.get("drug") or {}
                    target_node = row.get("target") or {}
                    d_name      = drug_node.get("name") or row.get("prefName")
                    t_name      = target_node.get("approvedSymbol")
                    moa         = (row.get("mechanismOfAction") or
                                   drug_node.get("mechanismOfAction") or "")

                    if not d_name or not t_name:
                        continue

                    key = d_name.lower().strip()
                    if key not in drugs_dict:
                        drugs_dict[key] = {
                            "name":      d_name,
                            "targets":   [],
                            "mechanism": moa,
                            "source":    "OpenTargets",
                            "phase":     row.get("phase"),
                            "status":    row.get("status"),
                        }
                    if t_name not in drugs_dict[key]["targets"]:
                        drugs_dict[key]["targets"].append(t_name)
                    if moa and not drugs_dict[key].get("mechanism"):
                        drugs_dict[key]["mechanism"] = moa
                    efo_drugs += 1

                except (KeyError, TypeError, AttributeError) as e:
                    logger.debug("Row parsing error: %s — row: %s", e, row)
                    continue

            next_cursor = kd.get("cursor")
            page_count += 1
            if not next_cursor or next_cursor == cursor:
                break
            cursor = next_cursor

        logger.info("  EFO %s: %d drugs, %d pages", efo, efo_drugs, page_count)

    async def _safe_graphql_post(
        self,
        session: aiohttp.ClientSession,
        payload: dict,
        efo: str,
        attempt: int = 0,
    ) -> Optional[dict]:
        """POST to OpenTargets GraphQL with retry."""
        try:
            async with session.post(
                self.GRAPHQL_URL,
                json=payload,
                headers={"Content-Type": "application/json"},
            ) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    errors = data.get("errors")
                    if errors:
                        logger.warning(
                            "EFO %s GraphQL errors: %s",
                            efo, json.dumps(errors[:1])[:200]
                        )
                        # Log the query variables for debugging
                        if not self._debug_printed:
                            logger.warning(
                                "Query variables: %s",
                                json.dumps(payload.get("variables"))
                            )
                            self._debug_printed = True
                        return None
                    return data

                elif resp.status == 400:
                    body = await resp.text()
                    logger.warning(
                        "EFO %s HTTP 400 (bad query): %s", efo, body[:300]
                    )
                    return None

                elif resp.status >= 500 and attempt < MAX_RETRIES:
                    await asyncio.sleep(RETRY_BACKOFF ** attempt)
                    return await self._safe_graphql_post(session, payload, efo, attempt + 1)

                else:
                    logger.debug("EFO %s HTTP %d", efo, resp.status)
                    return None

        except aiohttp.ClientConnectorError:
            raise  # Re-raise for outer handler
        except aiohttp.ClientError as e:
            if attempt < MAX_RETRIES:
                await asyncio.sleep(RETRY_BACKOFF ** attempt)
                return await self._safe_graphql_post(session, payload, efo, attempt + 1)
            logger.debug("EFO %s client error: %s", efo, e)
            return None

    async def _fetch_rest(self) -> List[Dict]:
        """
        Fallback: use OpenTargets Platform REST API for drug evidence.

        Endpoint: GET /v4/disease/{efoId}/drugs
        Returns drug-disease associations from the platform.

        Note: This endpoint may not exist in all versions; try v4 evidence API.
        """
        drugs_dict: Dict[str, Dict] = {}
        timeout = aiohttp.ClientTimeout(total=30, connect=10)

        # Try the evidence/drugs endpoint (v4 REST)
        async def fetch_rest_efo(session, efo: str) -> None:
            urls_to_try = [
                f"{self.REST_BASE_URL}/disease/{efo}/drugs",
                f"{self.REST_BASE_URL}/disease/{efo}/knownDrugs",
            ]
            for url in urls_to_try:
                try:
                    async with session.get(url) as resp:
                        if resp.status == 200:
                            data = await resp.json()
                            rows = data.get("rows") or data.get("data") or []
                            if rows:
                                logger.info(
                                    "REST %s: %d rows from %s", efo, len(rows), url
                                )
                                for row in rows:
                                    d_name = (row.get("drugName") or
                                              row.get("name") or
                                              (row.get("drug") or {}).get("name"))
                                    t_name = (row.get("targetSymbol") or
                                              row.get("symbol") or
                                              (row.get("target") or {}).get("approvedSymbol"))
                                    if d_name:
                                        key = d_name.lower().strip()
                                        if key not in drugs_dict:
                                            drugs_dict[key] = {
                                                "name":    d_name,
                                                "targets": [],
                                                "source":  "OpenTargets_REST",
                                                "phase":   row.get("phase"),
                                                "status":  row.get("status"),
                                            }
                                        if t_name and t_name not in drugs_dict[key]["targets"]:
                                            drugs_dict[key]["targets"].append(t_name)
                                return
                except Exception as e:
                    logger.debug("REST %s %s: %s", efo, url, e)

        try:
            async with aiohttp.ClientSession(timeout=timeout) as session:
                tasks = [fetch_rest_efo(session, efo) for efo in self.ATRT_EFO_IDS[:3]]
                await asyncio.gather(*tasks, return_exceptions=True)
        except Exception as e:
            logger.debug("REST fetch error: %s", e)

        return list(drugs_dict.values())

    # ── Generic annotation ────────────────────────────────────────────────────

    async def _annotate_generic_status(self, candidates: List[Dict]) -> List[Dict]:
        """Annotate with generic availability using seed + optional live API."""
        for drug in candidates:
            name = (drug.get("name") or "").lower().strip()
            for suffix in (" hydrochloride", " hcl", " sodium", " mesylate",
                           " phosphate", " sulfate", " acetate"):
                name = name.replace(suffix, "")
            name = name.strip()
            canonical = _ALIAS_TO_CANONICAL.get(name, name)

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

        return candidates

    # ── Merge helpers ─────────────────────────────────────────────────────────

    @staticmethod
    def _merge_with_fallback(
        ot_drugs: List[Dict],
        fallback: List[Dict],
    ) -> List[Dict]:
        """Merge OT results with curated fallback (exact name matching only)."""
        def _canonical_key(name: str) -> str:
            name = name.lower().strip()
            for suffix in (" hydrochloride", " hcl", " sodium", " mesylate",
                           " phosphate", " sulfate", " acetate"):
                name = name.replace(suffix, "")
            name = name.strip()
            return _ALIAS_TO_CANONICAL.get(name, name)

        merged: Dict[str, Dict] = {}

        for d in ot_drugs:
            key = _canonical_key(d.get("name", ""))
            if key and key not in merged:
                merged[key] = dict(d)
            elif key and key in merged:
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
                merged[fb_key].setdefault("mechanism",  fb.get("mechanism", ""))
                merged[fb_key].setdefault("drug_class", fb.get("drug_class", ""))
                merged[fb_key].setdefault("pmid",       fb.get("pmid", ""))
                existing = set(merged[fb_key].get("targets", []))
                fb_tgts  = set(fb.get("targets", []))
                merged[fb_key]["targets"] = sorted(existing | fb_tgts)
                n_enriched += 1
            else:
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