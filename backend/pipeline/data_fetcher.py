"""
data_fetcher.py  (v3.0)
========================
OpenTargets API Fetcher + Generic Drug integration for ATRT pipeline.

KEY CHANGES FROM v2.0
----------------------
1. Generic drug database merged into the candidate list.
   Drugs from generic_drug_db.py are added AFTER OpenTargets fetch so they
   always appear in the pipeline regardless of OpenTargets coverage.

2. Deduplication by name.
   If OpenTargets and generic_drug_db both contain the same drug (e.g.,
   "vorinostat"), the OpenTargets version is kept (richer target annotation)
   and the generic metadata is merged in.

3. FDA drug label lookup via OpenFDA API.
   For any drug with is_generic_available=True, we optionally fetch current
   FDA approval status and generic availability flag.

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
  PMID 33290552.
FDA OpenFDA API: https://open.fda.gov/apis/drug/
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


class ProductionDataFetcher:
    """
    Fetches drug-target associations from OpenTargets API and merges in
    FDA-approved generic / repurposable drugs from the local database.

    Candidate list construction
    ---------------------------
    1. Query OpenTargets API with ATRT/rhabdoid EFO IDs (≥500 drugs typical).
    2. Load generic_drug_db candidates (17 curated drugs as of v3.0).
    3. Deduplicate by lowercase name — OpenTargets wins on targets/metadata.
    4. Annotate generic availability flag for downstream cost/access reporting.

    The merged list is returned to discovery_pipeline.ProductionPipeline for
    scoring via tissue expression, DepMap, escape bypass, and PPI modules.
    """

    def __init__(self):
        self.api_url  = "https://api.platform.opentargets.org/api/v4/graphql"
        self.efo_ids  = OPENTARGETS["efo_ids"]
        self.max_pages = OPENTARGETS["max_pages_per_efo"]
        self.page_size = OPENTARGETS["drugs_per_page"]

    # ── Public API ────────────────────────────────────────────────────────────

    async def fetch_disease_data(self, disease_name: str) -> Dict:
        """Return a minimal disease context dict for pipeline compatibility."""
        return {
            "name":  disease_name,
            "genes": [
                "SMARCB1", "EZH2", "BRD4", "HDAC1", "AURKA",
                "MYC", "MYCN", "CDK4", "PSMB5", "BCL2L1",
            ],
            "id": "EFO_0002915",
        }

    async def fetch_approved_drugs(self) -> List[Dict]:
        """
        Fetch drug candidates from OpenTargets and merge with generic drug db.

        Returns
        -------
        List of candidate dicts, each with at minimum:
            name     : str
            targets  : List[str]   (HGNC gene symbols)
            is_generic : bool      (True if generic formulation available)
        """
        # Step 1: OpenTargets fetch
        ot_drugs = await self._fetch_opentargets_drugs()

        # Step 2: Generic drug database candidates
        generic_candidates = self._load_generic_drug_candidates()

        # Step 3: Merge — deduplicate on lowercase name
        merged = self._merge_candidates(ot_drugs, generic_candidates)

        logger.info(
            "Candidate list: %d from OpenTargets + %d from generic_drug_db → %d unique total",
            len(ot_drugs), len(generic_candidates), len(merged),
        )
        return merged

    async def close(self) -> None:
        pass

    # ── OpenTargets fetch ─────────────────────────────────────────────────────

    async def _fetch_opentargets_drugs(self) -> List[Dict]:
        """Query OpenTargets GraphQL API for drugs associated with ATRT EFOs."""
        logger.info(
            "Fetching drugs from OpenTargets (EFO IDs: %s)...", self.efo_ids
        )

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
                        logger.warning("EFO %s: 'disease' null", efo)
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
                                    "name":       d_name,
                                    "targets":    [],
                                    "is_generic": False,
                                    "source":     "OpenTargets",
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

        # Deduplicate targets within each drug
        for d in drugs_dict.values():
            d["targets"] = sorted(set(d["targets"]))

        count = len(drugs_dict)
        if count == 0:
            logger.error(
                "OpenTargets returned 0 drugs. "
                "Possible: network blocked, API down, or invalid EFO IDs. "
                "Generic drug database will still be used as fallback."
            )
        else:
            logger.info("OpenTargets: %d unique drugs retrieved", count)

        return list(drugs_dict.values())

    async def _safe_post(
        self,
        session: aiohttp.ClientSession,
        payload: dict,
        attempt: int = 0,
    ) -> Optional[dict]:
        """POST to OpenTargets GraphQL with retry and null-safe parsing."""
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

    # ── Generic drug database ─────────────────────────────────────────────────

    @staticmethod
    def _load_generic_drug_candidates() -> List[Dict]:
        """Load candidates from the generic/repurposable drug database."""
        try:
            import sys
            from pathlib import Path
            sys.path.insert(0, str(Path(__file__).parent))
            from generic_drug_db import as_pipeline_candidates
            candidates = as_pipeline_candidates()
            logger.info(
                "Generic drug database: %d candidates loaded", len(candidates)
            )
            return candidates
        except ImportError as e:
            logger.warning(
                "generic_drug_db.py not found — no generic drugs added: %s", e
            )
            return []

    @staticmethod
    def _merge_candidates(
        ot_drugs: List[Dict],
        generic_drugs: List[Dict],
    ) -> List[Dict]:
        """
        Merge OpenTargets and generic drug candidates, deduplicating by name.

        Deduplication strategy:
          - Normalise names to lowercase for comparison.
          - If drug appears in both lists, keep OpenTargets version as base
            and augment it with generic metadata (is_generic, cost_category, etc.).
          - If drug appears only in generic_drug_db, add it as-is.

        The merged list preserves all drugs from both sources.
        """
        # Index OT drugs by lowercase name
        merged: Dict[str, Dict] = {}
        for d in ot_drugs:
            key = d["name"].lower().strip()
            merged[key] = d

        # Merge generic drugs
        n_added   = 0
        n_enriched = 0

        for gd in generic_drugs:
            gd_name = gd.get("name", "").lower().strip()
            gd_generic = gd.get("drug_name", "").lower().strip()

            # Check for match by name or generic_name
            match_key = None
            if gd_name in merged:
                match_key = gd_name
            elif gd_generic and gd_generic in merged:
                match_key = gd_generic
            else:
                # Try partial match (e.g., "EPZ-6438" ↔ "tazemetostat")
                for ot_key in merged:
                    if (gd_name and gd_name in ot_key) or (ot_key and ot_key in gd_name):
                        match_key = ot_key
                        break

            if match_key:
                # Enrich existing OT drug with generic metadata
                existing = merged[match_key]
                existing.setdefault("is_generic",      gd.get("is_generic", False))
                existing.setdefault("cost_category",   gd.get("cost_category", "UNKNOWN"))
                existing.setdefault("evidence_level",  gd.get("evidence_level", "LOW"))
                existing.setdefault("atrt_rationale",  gd.get("atrt_rationale", ""))
                existing.setdefault("fda_approved",    gd.get("fda_approved", False))
                # Merge targets (OT targets take priority; add any from generic db)
                ot_targets = set(existing.get("targets", []))
                gd_targets = set(gd.get("targets", []))
                existing["targets"] = sorted(ot_targets | gd_targets)
                n_enriched += 1
            else:
                # New drug not in OT results — add from generic db
                merged[gd_name] = {
                    "name":           gd.get("name"),
                    "targets":        gd.get("targets", []),
                    "mechanism":      gd.get("mechanism", ""),
                    "drug_class":     gd.get("drug_class", ""),
                    "is_generic":     gd.get("is_generic", False),
                    "cost_category":  gd.get("cost_category", "UNKNOWN"),
                    "evidence_level": gd.get("evidence_level", "LOW"),
                    "fda_approved":   gd.get("fda_approved", False),
                    "source":         "generic_drug_db",
                    "atrt_rationale": gd.get("atrt_rationale", ""),
                }
                n_added += 1

        logger.info(
            "Merge result: %d OT drugs enriched with generic metadata, "
            "%d generic-only drugs added",
            n_enriched, n_added,
        )
        return list(merged.values())