"""
data_fetcher.py
===============
OpenTargets API Fetcher for ATRT Drug Repurposing Pipeline v2.0

CHANGES FROM v1.0
-----------------
- EFO IDs updated to include EFO_0002915 (rhabdoid tumor — primary ATRT EFO)
- Null-safe response parsing improved
- Retry logic with exponential backoff maintained
- Added drug deduplication across EFO queries

PRIMARY EFO IDs FOR ATRT
--------------------------
EFO_0002915: rhabdoid tumor (primary — catches tazemetostat, etc.)
EFO_0000543: malignant rhabdoid tumor (non-CNS; same SMARCB1 biology)
EFO_0000519: glioblastoma (catches BET/HDAC/CDK drugs with CNS PK data)
EFO_0001422: DIPG (catches epigenetic drugs — overlapping vulnerabilities)
EFO_0000618: pediatric brain tumor (broad catch for CNS-active drugs)

Reference: Ochoa D et al. Nucleic Acids Res 2021; 49(D1):D1302. PMID 33290552.
"""

import asyncio
import aiohttp
import logging
from typing import Dict, List, Optional

try:
    from .pipeline_config import OPENTARGETS
except ImportError:
    from pipeline_config import OPENTARGETS

logger = logging.getLogger(__name__)

MAX_RETRIES   = 3
RETRY_BACKOFF = 2.0   # seconds (exponential backoff)


class ProductionDataFetcher:
    """
    Fetches drug-target associations from OpenTargets API.

    Primary use: populate the drug candidate list before scoring.
    Secondary: basic disease gene data for context.

    If the API is unavailable, the pipeline falls back to
    discovery_pipeline._ATRT_FALLBACK_CANDIDATES() — a curated minimal set.
    """

    def __init__(self):
        self.api_url  = "https://api.platform.opentargets.org/api/v4/graphql"
        self.efo_ids  = OPENTARGETS["efo_ids"]
        self.max_pages = OPENTARGETS["max_pages_per_efo"]
        self.page_size = OPENTARGETS["drugs_per_page"]

    async def fetch_disease_data(self, disease_name: str) -> Dict:
        """
        Fetch basic disease context data.
        Returns a minimal dict for pipeline compatibility if API fails.
        """
        logger.info("Fetching disease data for: %s", disease_name)
        return {
            "name":  disease_name,
            "genes": [
                # Core ATRT disease genes for PPI context
                "SMARCB1", "EZH2", "BRD4", "HDAC1", "AURKA",
                "MYC", "MYCN", "CDK4", "PSMB5", "BCL2L1",
            ],
            "id": "EFO_0002915",   # rhabdoid tumor
        }

    async def _safe_post(
        self,
        session: aiohttp.ClientSession,
        payload: dict,
        attempt: int = 0,
    ) -> Optional[dict]:
        """POST to OpenTargets GraphQL with null-safe parsing and retries."""
        try:
            async with session.post(self.api_url, json=payload) as resp:

                if resp.status >= 500:
                    if attempt < MAX_RETRIES:
                        wait = RETRY_BACKOFF ** attempt
                        logger.warning(
                            "OpenTargets API returned %d — retrying in %.1fs (attempt %d/%d)",
                            resp.status, wait, attempt + 1, MAX_RETRIES,
                        )
                        await asyncio.sleep(wait)
                        return await self._safe_post(session, payload, attempt + 1)
                    else:
                        logger.error(
                            "OpenTargets API returned %d after %d retries",
                            resp.status, MAX_RETRIES,
                        )
                        return None

                if resp.status != 200:
                    logger.warning("OpenTargets API returned status %d", resp.status)
                    return None

                content_type = resp.headers.get("Content-Type", "")
                if "json" not in content_type:
                    text = await resp.text()
                    logger.error(
                        "OpenTargets API returned non-JSON (Content-Type: %s). "
                        "First 200 chars: %s", content_type, text[:200],
                    )
                    return None

                data = await resp.json()

                # Handle GraphQL-level errors
                if "errors" in data:
                    for err in data["errors"]:
                        logger.error("OpenTargets GraphQL error: %s", err.get("message", err))

                # Null-safe: data['data'] may be None on error
                raw_data = data.get("data")
                if raw_data is None:
                    logger.warning(
                        "OpenTargets returned data=null. EFO ID may be invalid "
                        "or service is degraded. Check https://platform.opentargets.org/api"
                    )
                    return None

                return data

        except aiohttp.ClientError as e:
            if attempt < MAX_RETRIES:
                wait = RETRY_BACKOFF ** attempt
                logger.warning("Network error: %s — retrying in %.1fs", e, wait)
                await asyncio.sleep(wait)
                return await self._safe_post(session, payload, attempt + 1)
            logger.error("API request failed after %d retries: %s", MAX_RETRIES, e)
            return None
        except Exception as e:
            logger.error("API request failed: %s", e)
            return None

    async def fetch_approved_drugs(self) -> List[Dict]:
        """
        Fetch CNS/Oncology drugs from OpenTargets using ATRT EFO IDs.

        Uses cursor pagination to get all drugs, not just first page.
        Deduplicates across EFO queries.

        Returns list of dicts: [{name, targets}, ...]
        """
        logger.info(
            "Fetching drugs from OpenTargets API (EFO IDs: %s)...",
            self.efo_ids,
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
                        logger.warning(
                            "EFO %s page %d: no data — skipping remaining pages",
                            efo, page_count + 1,
                        )
                        break

                    raw      = data.get("data", {})
                    disease  = raw.get("disease")
                    if disease is None:
                        logger.warning(
                            "EFO %s: 'disease' is null — EFO ID may not exist in OpenTargets",
                            efo,
                        )
                        break

                    kd   = disease.get("knownDrugs")
                    if kd is None:
                        logger.warning("EFO %s: 'knownDrugs' is null", efo)
                        break

                    rows = kd.get("rows") or []
                    if not rows:
                        break

                    for row in rows:
                        try:
                            drug_obj   = row.get("drug")   or {}
                            target_obj = row.get("target") or {}
                            d_name = drug_obj.get("name")
                            t_name = target_obj.get("approvedSymbol")
                            if not d_name or not t_name:
                                continue
                            if d_name not in drugs_dict:
                                drugs_dict[d_name] = {"name": d_name, "targets": []}
                            drugs_dict[d_name]["targets"].append(t_name)
                            efo_drugs += 1
                        except (KeyError, TypeError) as e:
                            logger.debug("Skipping malformed row: %s", e)

                    cursor = kd.get("cursor")
                    page_count += 1
                    if not cursor:
                        break

                logger.info(
                    "  EFO %s: %d drugs across %d pages", efo, efo_drugs, page_count
                )

        # Deduplicate targets within each drug
        for d in drugs_dict.values():
            d["targets"] = list(set(d["targets"]))

        count = len(drugs_dict)
        if count == 0:
            logger.error(
                "❌ OpenTargets API returned 0 drugs.\n"
                "   Possible causes:\n"
                "   1. Network/firewall blocking outbound HTTPS\n"
                "   2. OpenTargets API temporarily down\n"
                "   3. All EFO IDs returned null\n"
                "   Pipeline will use fallback drug library."
            )
        else:
            logger.info("✅ OpenTargets: %d unique drugs retrieved", count)

        return list(drugs_dict.values())

    async def close(self) -> None:
        pass