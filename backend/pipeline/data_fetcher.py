"""
data_fetcher.py  (v4.0)
========================
OpenTargets API Fetcher for ATRT pipeline.

KEY CHANGES FROM v3.0
----------------------
1. Generic drug filtering is now DYNAMIC via generic_drug_fetcher.py.
   The old generic_drug_db.py had a hardcoded list that became stale.
   Now we query FDA Drugs@FDA + RxNorm at runtime and cache results locally.

2. The generic drug candidate list is no longer manually maintained.
   We query OpenFDA for any drug that matches ATRT-relevant EFO IDs,
   then dynamically annotate which ones have generics available.

3. A curated ATRT-relevant drug list is used as a guaranteed fallback
   when OpenTargets API is unavailable, but generic STATUS is always
   checked dynamically — not hardcoded.

PRIMARY EFO IDs FOR ATRT
--------------------------
EFO_0002915 : rhabdoid tumor (primary — catches tazemetostat etc.)
EFO_0000543 : malignant rhabdoid tumor (non-CNS; identical SMARCB1 biology)
EFO_0000519 : glioblastoma (catches BET/HDAC/CDK drugs with CNS PK data)
EFO_0001422 : DIPG (catches epigenetic drugs — overlapping vulnerabilities)
EFO_0000618 : pediatric brain tumor (broad catch for CNS drugs)

References
-----------
Ochoa D et al. (2021). Open Targets Platform. Nucleic Acids Res 49(D1):D1302.
FDA Drugs@FDA: https://api.fda.gov/drug/drugsfda.json
RxNorm API:    https://rxnav.nlm.nih.gov/REST/
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

# ─────────────────────────────────────────────────────────────────────────────
# CURATED ATRT FALLBACK CANDIDATES
#
# Used ONLY when OpenTargets API returns 0 results (network failure etc).
# Generic STATUS is NOT hardcoded here — it is always checked via
# generic_drug_fetcher.py at runtime.
#
# This list covers the drugs most validated in published ATRT biology.
# Sources: Torchia 2015, Knutson 2013, Sredni 2017, Geoerger 2017,
#          Frühwald 2020, Johann 2016.
# ─────────────────────────────────────────────────────────────────────────────

ATRT_CURATED_FALLBACK: List[Dict] = [
    # EZH2 synthetic lethality (Knutson 2013 PNAS PMID 23620515)
    {"name": "Tazemetostat",     "targets": ["EZH2"],
     "mechanism": "EZH2 inhibitor", "pmid": "23620515"},

    # Pan-HDAC (Torchia 2015 Cancer Cell PMID 26609405)
    {"name": "Panobinostat",     "targets": ["HDAC1", "HDAC2", "HDAC3", "HDAC6"],
     "mechanism": "pan-HDAC inhibitor", "pmid": "26609405"},
    {"name": "Vorinostat",       "targets": ["HDAC1", "HDAC2", "HDAC3", "HDAC6"],
     "mechanism": "pan-HDAC inhibitor", "pmid": "26609405"},
    {"name": "Valproic acid",    "targets": ["HDAC1", "HDAC2", "HDAC3"],
     "mechanism": "HDAC inhibitor (repurposed anticonvulsant)", "pmid": "11742990"},
    {"name": "Entinostat",       "targets": ["HDAC1", "HDAC2", "HDAC3"],
     "mechanism": "class I HDAC inhibitor"},

    # Aurora kinase (Sredni 2017 Pediatric Blood Cancer PMID 28544500)
    {"name": "Alisertib",        "targets": ["AURKA"],
     "mechanism": "aurora kinase A inhibitor", "pmid": "28544500"},

    # BET bromodomain (Geoerger 2017 Clin Cancer Res PMID 28108534)
    {"name": "Birabresib",       "targets": ["BRD4", "BRD2", "BRD3"],
     "mechanism": "BET bromodomain inhibitor", "pmid": "28108534"},

    # CDK4/6
    {"name": "Abemaciclib",      "targets": ["CDK4", "CDK6"],
     "mechanism": "CDK4/CDK6 inhibitor"},
    {"name": "Palbociclib",      "targets": ["CDK4", "CDK6"],
     "mechanism": "CDK4/CDK6 inhibitor"},
    {"name": "Ribociclib",       "targets": ["CDK4", "CDK6"],
     "mechanism": "CDK4/CDK6 inhibitor"},

    # Proteasome (Lin 2019 Sci Transl Med)
    {"name": "Marizomib",        "targets": ["PSMB5", "PSMB2", "PSMB1"],
     "mechanism": "proteasome inhibitor"},
    {"name": "Bortezomib",       "targets": ["PSMB5", "PSMB1"],
     "mechanism": "proteasome inhibitor (boronic acid)"},

    # SHH subgroup (Torchia 2015; Johann 2016)
    {"name": "Vismodegib",       "targets": ["SMO"],
     "mechanism": "smoothened inhibitor"},
    {"name": "Sonidegib",        "targets": ["SMO"],
     "mechanism": "smoothened inhibitor"},
    {"name": "Itraconazole",     "targets": ["SMO", "PTCH1"],
     "mechanism": "SMO inhibitor (repurposed antifungal)"},

    # DRD2 / TRAIL
    {"name": "ONC201",           "targets": ["DRD2", "CLPB"],
     "mechanism": "DRD2 antagonist / TRAIL inducer"},

    # PI3K/mTOR
    {"name": "Paxalisib",        "targets": ["PIK3CA", "PIK3CD", "PIK3CG", "MTOR"],
     "mechanism": "PI3K inhibitor (CNS-penetrant)"},
    {"name": "Sirolimus",        "targets": ["MTOR", "RPTOR"],
     "mechanism": "mTOR inhibitor (repurposed immunosuppressant)"},

    # Repurposable generics with epigenetic/MYC rationale
    {"name": "Metformin",        "targets": ["PRKAB1", "PRKAB2"],
     "mechanism": "AMPK activator (repurposed antidiabetic)"},
    {"name": "Chloroquine",      "targets": ["ATP6V0A1"],
     "mechanism": "autophagy inhibitor (repurposed antimalarial)"},
    {"name": "Hydroxychloroquine","targets": ["ATP6V0A1"],
     "mechanism": "autophagy inhibitor (repurposed antimalarial)"},

    # Differentiation (TYR subgroup rationale)
    {"name": "All-trans retinoic acid", "targets": ["RARA", "RARB", "RARG"],
     "mechanism": "retinoid receptor agonist / differentiation therapy"},

    # SHH via GLI (ATRT-SHH subgroup)
    {"name": "Arsenic trioxide", "targets": ["GLI1", "GLI2"],
     "mechanism": "GLI inhibitor (repurposed from APL therapy)"},

    # Anti-apoptotic
    {"name": "Venetoclax",       "targets": ["BCL2"],
     "mechanism": "BCL-2 inhibitor"},
]


class ProductionDataFetcher:
    """
    Fetches drug-target associations from OpenTargets API and annotates
    generic drug availability dynamically via FDA/RxNorm APIs.

    Pipeline
    --------
    1. Query OpenTargets with ATRT/rhabdoid EFO IDs (≥500 drugs typical).
    2. Fall back to ATRT_CURATED_FALLBACK if OpenTargets is unreachable.
    3. Deduplicate candidates by lowercase name.
    4. Annotate ALL candidates with generic availability via GenericDrugFetcher.
       This is a LIVE query — not hardcoded.

    The caller (discovery_pipeline.py) can then:
      - Filter to generic_only=True for access/cost analysis
      - Use has_generic flag for reporting
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
        Fetch drug candidates and annotate with dynamic generic status.

        Parameters
        ----------
        annotate_generics : bool
            If True, query FDA/RxNorm for each drug's generic availability.
            Set to False in unit tests or when network is unavailable.

        Returns
        -------
        List of candidate dicts, each with:
            name          : str
            targets       : List[str]   (HGNC gene symbols)
            has_generic   : bool        (dynamic FDA/RxNorm check)
            cost_category : str         (GENERIC | BRAND | INVESTIGATIONAL)
            generic_source: str         (FDA_ANDA | RXNORM | SEED | FALLBACK)
        """
        # Step 1: OpenTargets fetch
        ot_drugs = await self._fetch_opentargets_drugs()

        # Step 2: If OT failed, use curated fallback (still annotate generics live)
        if not ot_drugs:
            logger.warning(
                "OpenTargets returned 0 drugs — using ATRT curated fallback library (%d drugs). "
                "Generic status will still be checked via FDA/RxNorm.",
                len(ATRT_CURATED_FALLBACK),
            )
            ot_drugs = [dict(d) for d in ATRT_CURATED_FALLBACK]

        # Step 3: Merge curated fallback to ensure key ATRT drugs are always present
        merged = self._merge_with_fallback(ot_drugs, ATRT_CURATED_FALLBACK)

        # Step 4: Dynamically annotate generic availability
        if annotate_generics:
            try:
                from .generic_drug_fetcher import annotate_with_generic_status
                merged = await annotate_with_generic_status(merged)
                n_generic = sum(1 for d in merged if d.get("has_generic"))
                logger.info(
                    "Generic annotation complete: %d/%d drugs have generics available",
                    n_generic, len(merged),
                )
            except Exception as e:
                logger.warning(
                    "Generic drug annotation failed (%s) — all set to UNKNOWN. "
                    "Pipeline continues without generic filtering.",
                    e,
                )
                for d in merged:
                    d.setdefault("has_generic", False)
                    d.setdefault("cost_category", "UNKNOWN")
        else:
            for d in merged:
                d.setdefault("has_generic", False)
                d.setdefault("cost_category", "UNKNOWN")

        logger.info("Total candidates after merge + annotation: %d", len(merged))
        return merged

    async def close(self) -> None:
        pass

    # ── OpenTargets fetch ─────────────────────────────────────────────────────

    async def _fetch_opentargets_drugs(self) -> List[Dict]:
        logger.info("Fetching from OpenTargets (EFOs: %s) ...", self.efo_ids)

        query = """
        query atrtDrugs($efoId: String!, $cursor: String, $size: Int) {
          disease(efoId: $efoId) {
            knownDrugs(cursor: $cursor, size: $size) {
              cursor
              rows {
                drug { name }
                target { approvedSymbol }
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

                        raw     = data.get("data", {})
                        disease = raw.get("disease")
                        if disease is None:
                            break

                        kd   = disease.get("knownDrugs")
                        rows = (kd or {}).get("rows") or []
                        if not rows:
                            break

                        for row in rows:
                            try:
                                d_name = (row.get("drug") or {}).get("name")
                                t_name = (row.get("target") or {}).get("approvedSymbol")
                                if not d_name or not t_name:
                                    continue
                                key = d_name.lower()
                                if key not in drugs_dict:
                                    drugs_dict[key] = {
                                        "name":    d_name,
                                        "targets": [],
                                        "source":  "OpenTargets",
                                    }
                                drugs_dict[key]["targets"].append(t_name)
                                efo_drugs += 1
                            except (KeyError, TypeError):
                                continue

                        cursor = (kd or {}).get("cursor")
                        page_count += 1
                        if not cursor:
                            break

                    logger.info("  EFO %s: %d drugs, %d pages", efo, efo_drugs, page_count)

        except Exception as e:
            logger.warning("OpenTargets fetch failed: %s", e)
            return []

        # Deduplicate targets within each drug
        for d in drugs_dict.values():
            d["targets"] = sorted(set(d["targets"]))

        logger.info("OpenTargets: %d unique drugs retrieved", len(drugs_dict))
        return list(drugs_dict.values())

    async def _safe_post(
        self,
        session: aiohttp.ClientSession,
        payload: dict,
        attempt: int = 0,
    ) -> Optional[dict]:
        try:
            async with session.post(self.api_url, json=payload) as resp:
                if resp.status >= 500:
                    if attempt < MAX_RETRIES:
                        await asyncio.sleep(RETRY_BACKOFF ** attempt)
                        return await self._safe_post(session, payload, attempt + 1)
                    return None
                if resp.status != 200:
                    return None
                data = await resp.json()
                if data.get("data") is None:
                    return None
                return data
        except aiohttp.ClientError as e:
            if attempt < MAX_RETRIES:
                await asyncio.sleep(RETRY_BACKOFF ** attempt)
                return await self._safe_post(session, payload, attempt + 1)
            logger.error("OpenTargets request failed: %s", e)
            return None
        except Exception as e:
            logger.error("OpenTargets request error: %s", e)
            return None

    # ── Merge helpers ─────────────────────────────────────────────────────────

    @staticmethod
    def _merge_with_fallback(
        ot_drugs: List[Dict],
        fallback: List[Dict],
    ) -> List[Dict]:
        """
        Merge OpenTargets results with curated fallback list.

        Strategy:
        - Index OT drugs by lowercase name.
        - For each fallback drug, if it's already in OT results, enrich the OT
          entry with mechanism/pmid from the fallback.
        - If not in OT results, add the fallback drug directly.
        - This ensures key ATRT drugs (tazemetostat, alisertib etc.) are always
          present regardless of what OpenTargets returned.
        """
        merged: Dict[str, Dict] = {}
        for d in ot_drugs:
            key = d["name"].lower().strip()
            merged[key] = d

        n_added = 0
        n_enriched = 0

        for fb in fallback:
            fb_name = fb.get("name", "").lower().strip()

            # Try exact match
            if fb_name in merged:
                existing = merged[fb_name]
                existing.setdefault("mechanism", fb.get("mechanism", ""))
                existing.setdefault("pmid", fb.get("pmid", ""))
                # Merge targets
                ot_targets = set(existing.get("targets", []))
                fb_targets = set(fb.get("targets", []))
                existing["targets"] = sorted(ot_targets | fb_targets)
                n_enriched += 1
                continue

            # Try partial match (e.g. "all-trans retinoic acid" ↔ "tretinoin")
            found = False
            for ot_key in list(merged.keys()):
                if (fb_name and fb_name[:8] in ot_key) or (ot_key[:8] in fb_name):
                    existing = merged[ot_key]
                    existing.setdefault("mechanism", fb.get("mechanism", ""))
                    existing.setdefault("pmid", fb.get("pmid", ""))
                    fb_targets = set(fb.get("targets", []))
                    existing["targets"] = sorted(
                        set(existing.get("targets", [])) | fb_targets
                    )
                    n_enriched += 1
                    found = True
                    break

            if not found:
                # New drug not in OT results
                merged[fb_name] = {
                    "name":      fb.get("name"),
                    "targets":   fb.get("targets", []),
                    "mechanism": fb.get("mechanism", ""),
                    "pmid":      fb.get("pmid", ""),
                    "source":    "curated_fallback",
                }
                n_added += 1

        logger.info(
            "Merge: %d OT drugs + %d fallback → %d unique total "
            "(%d enriched, %d added from fallback)",
            len(ot_drugs), len(fallback), len(merged), n_enriched, n_added,
        )
        return list(merged.values())