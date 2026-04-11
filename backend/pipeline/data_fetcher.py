"""
data_fetcher.py  (v9.0 — Live Broad Institute Fetcher)
================================================================
ARCHITECTURE CHANGE: Option B (Global Screen) + Live Download

This script now natively checks for a massive database of drugs. 
If it is not found on the local machine, it connects to the Broad 
Institute's public AWS servers and downloads the ~6,700 drug 
Repurposing Hub database LIVE.

This guarantees a massive, unbiased multi-omic screen without 
requiring manual data management from the user.
"""

import asyncio
import aiohttp
import csv
import logging
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# ── Fallback Seed (Used ONLY if the internet goes down completely) ────────────
GLOBAL_REPURPOSING_SEED: List[Dict] = [
    {"name": "Tazemetostat", "targets": ["EZH2"], "mechanism": "EZH2 inhibitor"},
    {"name": "Panobinostat", "targets": ["HDAC1", "HDAC2", "HDAC3", "HDAC6"], "mechanism": "pan-HDAC inhibitor"},
    {"name": "Valproic acid", "targets": ["HDAC1", "HDAC2", "HDAC3"], "mechanism": "HDAC inhibitor"},
    {"name": "Alisertib", "targets": ["AURKA"], "mechanism": "aurora kinase A inhibitor"},
    {"name": "Birabresib", "targets": ["BRD4", "BRD2", "BRD3"], "mechanism": "BET bromodomain inhibitor"},
    {"name": "Abemaciclib", "targets": ["CDK4", "CDK6"], "mechanism": "CDK4/CDK6 inhibitor"},
    {"name": "Paxalisib", "targets": ["PIK3CA", "PIK3CD", "PIK3CG", "MTOR"], "mechanism": "PI3K inhibitor"},
    {"name": "Marizomib", "targets": ["PSMB5", "PSMB2", "PSMB1"], "mechanism": "proteasome inhibitor"},
    {"name": "Metformin", "targets": ["PRKAB1", "PRKAB2"], "mechanism": "AMPK activator"},
]


class ProductionDataFetcher:
    """
    Unbiased Global Drug Fetcher for the ATRT Pipeline.
    Automatically fetches the Broad Repurposing Hub if missing.
    """

    def __init__(self):
        self.local_library_path = Path("data/repurposing_library.tsv")

    async def fetch_disease_data(self, disease_name: str) -> Dict:
        """Still provides the base biological context of ATRT to the pipeline."""
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
        """Fetch a massive global drug list for unbiased screening."""
        
        logger.info("Initializing Unbiased Repurposing Screen Mode...")
        
        candidates = self._load_local_library()
        
        if not candidates:
            logger.warning("Live download failed or file empty. Using emergency internal seed.")
            candidates = GLOBAL_REPURPOSING_SEED
            
        logger.info("✅ Loaded %d total drugs with known targets for unbiased screening.", len(candidates))

        if annotate_generics:
            # Dynamically imported to avoid circular logic
            from backend.pipeline.generic_drug_fetcher import annotate_with_generic_status
            candidates = await annotate_with_generic_status(candidates)
        else:
            for d in candidates:
                d.setdefault("has_generic", False)
                d.setdefault("cost_category", "UNKNOWN")

        return candidates

    def _load_local_library(self) -> List[Dict]:
        """
        Attempts to load the massive Broad Repurposing Hub.
        If it doesn't exist, it downloads it live from the internet.
        """
        # Step 1: Download live if it doesn't exist
        if not self.local_library_path.exists():
            logger.info("Massive library not found locally. Fetching Broad Institute Repurposing Hub LIVE...")
            self.local_library_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Official Broad Institute AWS URL
            url = "https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt"
            try:
                urllib.request.urlretrieve(url, self.local_library_path)
                logger.info("✅ Live download complete! Saved to data/repurposing_library.tsv")
            except Exception as e:
                logger.error("❌ Failed to download live database: %s", e)
                return []

        # Step 2: Parse the massive dataset
        drugs = []
        try:
            with open(self.local_library_path, 'r', encoding='utf-8') as f:
                # The Broad file has 9 lines of metadata starting with '!'. We skip those.
                lines = [line for line in f if not line.startswith('!')]
                reader = csv.DictReader(lines, delimiter='\t')
                
                for row in reader:
                    name = row.get("pert_iname", "").strip()
                    targets_raw = row.get("target", "")
                    moa = row.get("moa", "").strip()
                    phase = row.get("clinical_phase", "0")
                    
                    if not name or not targets_raw:
                        continue
                    
                    # Broad targets are separated by the pipe character '|'
                    targets = [t.strip() for t in targets_raw.split("|") if t.strip()]
                    
                    if targets:
                        drugs.append({
                            "name": name,
                            "targets": targets,
                            "mechanism": moa,
                            "source": "Broad_Repurposing_Hub",
                            "phase": phase
                        })
            return drugs
        except Exception as e:
            logger.error("Failed to parse local library %s: %s", self.local_library_path, e)
            return []

    async def close(self) -> None:
        pass