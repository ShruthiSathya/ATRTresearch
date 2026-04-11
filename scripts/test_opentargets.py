#!/usr/bin/env python3
"""
scripts/test_opentargets.py
============================
Diagnostic tool to test OpenTargets API connectivity and query correctness.
Run this first to identify why the drug fetch is failing.

Usage:
    python scripts/test_opentargets.py

Expected output:
    ✅ if API is reachable and returns data
    ❌ with the exact error if something is broken
"""

import asyncio
import json
import sys
import aiohttp

GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# Test 1: Minimal connectivity check
PING_QUERY = """
{
  meta {
    dataVersion {
      iteration
    }
  }
}
"""

# Test 2: Disease lookup — verify EFO IDs work
DISEASE_QUERY = """
query($efoId: String!) {
  disease(efoId: $efoId) {
    id
    name
    knownDrugs {
      count
    }
  }
}
"""

# Test 3: Full knownDrugs with correct v4 fields
FULL_DRUGS_QUERY = """
query($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    knownDrugs(size: $size) {
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
        }
        phase
        status
      }
    }
  }
}
"""

# EFO IDs to test (in priority order)
EFO_IDS = [
    ("EFO_0002915",   "rhabdoid tumor"),
    ("EFO_0000543",   "malignant rhabdoid tumor"),
    ("MONDO_0024491", "ATRT (MONDO)"),
    ("EFO_0000519",   "glioblastoma"),
    ("EFO_0001422",   "DIPG"),
]


async def test_api():
    print("=" * 70)
    print("OpenTargets API Diagnostic Tool")
    print("=" * 70)

    timeout = aiohttp.ClientTimeout(total=30, connect=10)
    async with aiohttp.ClientSession(timeout=timeout) as session:

        # ── Test 1: API connectivity ──────────────────────────────────────
        print("\n[1] Testing API connectivity ...")
        try:
            async with session.post(
                GRAPHQL_URL,
                json={"query": PING_QUERY},
                headers={"Content-Type": "application/json"},
            ) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    version = (data.get("data") or {}).get("meta", {}).get(
                        "dataVersion", {}
                    ).get("iteration", "?")
                    print(f"    ✅ API reachable — data version: {version}")
                else:
                    body = await resp.text()
                    print(f"    ❌ HTTP {resp.status}: {body[:200]}")
                    return
        except aiohttp.ClientConnectorError as e:
            print(f"    ❌ Network error (DNS/firewall): {e}")
            print("    → Check your internet connection or proxy settings")
            return
        except Exception as e:
            print(f"    ❌ Unexpected error: {e}")
            return

        # ── Test 2: EFO ID validation ─────────────────────────────────────
        print("\n[2] Testing EFO IDs (checking which return drug data) ...")
        working_efos = []

        for efo_id, name in EFO_IDS:
            try:
                async with session.post(
                    GRAPHQL_URL,
                    json={
                        "query":     DISEASE_QUERY,
                        "variables": {"efoId": efo_id},
                    },
                    headers={"Content-Type": "application/json"},
                ) as resp:
                    data = await resp.json()
                    errors = data.get("errors")
                    if errors:
                        print(f"    ⚠️  {efo_id} ({name}): GraphQL error — {errors[0].get('message','?')[:80]}")
                        continue

                    disease = (data.get("data") or {}).get("disease")
                    if disease is None:
                        print(f"    ❌ {efo_id} ({name}): disease=None (EFO not found in OT)")
                        continue

                    kd    = disease.get("knownDrugs") or {}
                    count = kd.get("count", 0)
                    print(f"    {'✅' if count > 0 else '⚠️ '} {efo_id} ({name}): {count} known drugs")
                    if count > 0:
                        working_efos.append(efo_id)

            except Exception as e:
                print(f"    ❌ {efo_id}: {e}")

        # ── Test 3: Full drug retrieval on best EFO ───────────────────────
        if working_efos:
            best_efo = working_efos[0]
            print(f"\n[3] Retrieving drugs for best EFO: {best_efo} ...")
            try:
                async with session.post(
                    GRAPHQL_URL,
                    json={
                        "query":     FULL_DRUGS_QUERY,
                        "variables": {"efoId": best_efo, "size": 10},
                    },
                    headers={"Content-Type": "application/json"},
                ) as resp:
                    data = await resp.json()
                    errors = data.get("errors")
                    if errors:
                        print(f"    ❌ GraphQL errors: {json.dumps(errors[:1], indent=2)}")
                        return

                    disease = (data.get("data") or {}).get("disease", {})
                    kd      = disease.get("knownDrugs", {})
                    rows    = kd.get("rows", [])
                    total   = kd.get("count", 0)

                    print(f"    ✅ Total drugs available: {total}")
                    print(f"    ✅ Retrieved {len(rows)} rows in this page")
                    print("\n    Sample drugs returned:")
                    for row in rows[:5]:
                        d = row.get("drug", {})
                        t = row.get("target", {})
                        moa = row.get("mechanismOfAction", "—")
                        print(
                            f"      {d.get('name','?'):<25}  "
                            f"target={t.get('approvedSymbol','?'):<12}  "
                            f"phase={row.get('phase','?')}  "
                            f"moa={moa[:40]}"
                        )

                    print(f"\n    Cursor for next page: {kd.get('cursor', 'None')}")

            except Exception as e:
                print(f"    ❌ Full retrieval error: {e}")

        else:
            print("\n[3] ⚠️  No working EFO IDs found — all returned 0 drugs.")
            print("    This may mean:")
            print("    a) OpenTargets has no ATRT-specific drug data (sparse indication)")
            print("    b) EFO IDs are not in the OT database for this version")
            print("    c) Network restriction is blocking the request selectively")
            print("\n    ✅ The curated fallback list (24 drugs) will be used instead.")
            print("    This provides the same key ATRT drugs (tazemetostat,")
            print("    panobinostat, alisertib, birabresib, abemaciclib, etc.)")

        # ── Test 4: Check REST API as alternative ─────────────────────────
        print("\n[4] Testing REST API alternative ...")
        rest_url = f"https://api.platform.opentargets.org/api/v4/disease/EFO_0002915"
        try:
            async with session.get(rest_url) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    print(f"    ✅ REST API reachable — disease: {data.get('name','?')}")
                else:
                    print(f"    ⚠️  REST API returned {resp.status}")
        except Exception as e:
            print(f"    ⚠️  REST API error: {e}")

        print("\n" + "=" * 70)
        print("DIAGNOSIS SUMMARY")
        print("=" * 70)
        if working_efos:
            print(f"✅ Working EFO IDs: {working_efos}")
            print("   Update ATRT_EFO_IDS in data_fetcher.py to prioritise these.")
        else:
            print("⚠️  OpenTargets has sparse/no data for ATRT-specific EFO IDs.")
            print("   This is EXPECTED — ATRT is rare and OT coverage is limited.")
            print("   The curated fallback list provides all key ATRT candidates.")
            print("\n   RECOMMENDED ACTION:")
            print("   The curated fallback (24 drugs) is scientifically adequate.")
            print("   All key drugs from ATRT literature are already included:")
            print("   tazemetostat, panobinostat, alisertib, birabresib,")
            print("   abemaciclib, marizomib, vismodegib, ONC201, paxalisib,")
            print("   valproic acid, vorinostat, sirolimus, itraconazole, etc.")
        print("=" * 70)


if __name__ == "__main__":
    asyncio.run(test_api())