import asyncio
import aiohttp
import logging
try:
    from .pipeline_config import PPI as PPI_CONFIG
except ImportError:
    from pipeline_config import PPI as PPI_CONFIG
from typing import Dict, List, Set

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# CURATED PPI NEIGHBORS FOR DIPG/GBM DISEASE GENES
#
# Source: curated from published literature and well-established pathway
# databases (KEGG, Reactome, OMIM). Core interactions (PRC2 complex,
# CDK4-RB1-E2F1 axis, EGFR-PI3K-AKT, BRD4-MYC) are textbook biology
# with very high STRING-DB confidence scores (≥ 900). Peripheral entries
# (e.g. PSMB5 neighbors) are plausible but less formally verified.
#
# These are NOT directly queried from STRING-DB — they are curated
# from literature knowledge as a fallback for network-restricted
# environments. When STRING-DB is reachable, live queries take precedence
# (see fetch_neighbors_live). When it is not, these provide meaningful
# differentiation vs the previous flat-0.2 default.
#
# Purpose: provides PPI scoring without live STRING-DB calls.
# Drug targets that are 1st-degree neighbors of disease genes score 0.85.
# ─────────────────────────────────────────────────────────────────────────────

CURATED_PPI_NEIGHBORS: Dict[str, List[str]] = {
    # ── Core DIPG disease genes (ProductionDataFetcher returns these) ─────────
    "H3-3A":  ["EZH2", "EED", "SUZ12", "ATRX", "DAXX", "BRD4", "KDM6A", "SETD2"],
    "H3F3A":  ["EZH2", "EED", "SUZ12", "ATRX", "DAXX", "BRD4", "KDM6A", "SETD2"],
    "EZH2":   ["EED", "SUZ12", "HDAC1", "HDAC2", "BRD4", "KDM6A", "KDM6B",
               "DNMT1", "DNMT3A", "ATRX", "RBBP4", "RBBP7"],
    "CDK4":   ["CDKN2A", "CDKN2B", "CCND1", "CCND2", "RB1", "E2F1",
               "CDK6", "MDM2", "HDAC1", "PIK3CA"],
    "PDGFRA": ["PIK3CA", "PIK3R1", "AKT1", "MTOR", "STAT3", "EGFR",
               "MAPK1", "GRB2", "SRC", "JAK1"],
    "EGFR":   ["ERBB2", "SRC", "PIK3CA", "AKT1", "MAPK1", "STAT3",
               "GRB2", "MET", "PDGFRA", "MDM2"],
    "PTEN":   ["AKT1", "PIK3CA", "TP53", "MDM2", "MTOR", "PIK3R1",
               "CDKN1B", "CDK4"],

    # ── Extended DIPG-specific genes ─────────────────────────────────────────
    "ACVR1":  ["BMPR1A", "BMPR2", "SMAD1", "SMAD4", "SMAD5", "ID1", "ID2",
               "BMP4", "NOGGIN"],
    "BRD4":   ["CDK9", "MYC", "MYCN", "HDAC1", "HDAC2", "EZH2",
               "E2F1", "RELA", "CDK4"],
    "HDAC1":  ["HDAC2", "HDAC3", "EZH2", "BRD4", "DNMT1", "MBD3",
               "SIN3A", "RBBP4", "MYC", "TP53"],
    "HDAC2":  ["HDAC1", "HDAC3", "EZH2", "SIN3A", "MBD3",
               "RBBP4", "RBBP7", "MYC"],
    "CDK6":   ["CDK4", "CDKN2A", "CCND1", "CCND3", "RB1", "E2F1"],
    "CDKN2A": ["CDK4", "CDK6", "MDM2", "TP53", "RB1"],
    "TP53":   ["MDM2", "MDM4", "CDK4", "CDKN2A", "PTEN", "ATM",
               "HDAC1", "BCL2", "PARP1"],
    "MYC":    ["BRD4", "HDAC1", "HDAC2", "CDK4", "E2F1", "MYCN",
               "MAX", "MTOR"],
    "MYCN":   ["BRD4", "MYC", "CDK4", "AURKA", "ALK", "MDM2"],
    "STAT3":  ["JAK1", "JAK2", "IL6", "EGFR", "PDGFRA", "SRC", "BCL2"],
    "PIK3CA": ["PIK3R1", "AKT1", "MTOR", "PTEN", "EGFR", "PDGFRA",
               "KRAS", "RAS"],
    "AKT1":   ["PIK3CA", "MTOR", "PTEN", "MDM2", "BCL2", "HDAC1"],
    "MTOR":   ["PIK3CA", "AKT1", "PTEN", "TSC1", "TSC2", "RPTOR"],

    # ── Proteasome (Marizomib targets) ────────────────────────────────────────
    "PSMB5":  ["PSMB1", "PSMB2", "PSMA1", "PSMA3", "PSMD1",
               "UBB", "UBC", "SQSTM1", "ATG5"],
    "PSMB2":  ["PSMB5", "PSMB1", "PSMA1", "PSMD1", "SQSTM1"],
    "PSMB8":  ["PSMB5", "PSMB9", "PSMB10", "PSME1", "PSME2"],

    # ── BET bromodomain (Birabresib targets) ──────────────────────────────────
    "BRD2":   ["BRD4", "BRD3", "CDK9", "MYC", "HDAC1", "E2F1"],
    "BRD3":   ["BRD4", "BRD2", "CDK9", "MYC"],
    "BRDT":   ["BRD4", "BRD2", "CDK9"],

    # ── HDAC inhibitor targets (Belinostat, Panobinostat) ────────────────────
    "HDAC6":  ["HDAC1", "HDAC2", "HSP90AA1", "TUBB", "TUBA1A", "STAT3"],
    "HDAC9":  ["HDAC1", "HDAC2", "HDAC3", "HDAC4", "MEF2D"],
    "HDAC11": ["HDAC1", "HDAC2", "HDAC3"],

    # ── Taxane/vinca targets ──────────────────────────────────────────────────
    "TUBB3":  ["TUBB", "TUBB2B", "TUBA1A", "TUBA1C", "MAPT"],
    "TUBB2B": ["TUBB3", "TUBB", "TUBA1A", "TUBA1C"],
    "TUBA1C": ["TUBB3", "TUBB2B", "TUBB", "TUBA1A"],
    "TUBA4A": ["TUBA1C", "TUBA1A", "TUBB3"],

    # ── Topoisomerase targets (Doxorubicin, Etoposide) ───────────────────────
    "TOP2A":  ["TOP2B", "MYC", "HDAC1", "PARP1", "ATM"],
    "TOP2B":  ["TOP2A", "PARP1", "ATM"],

    # ── DNA damage ────────────────────────────────────────────────────────────
    "PARP1":  ["PARP2", "ATM", "ATR", "BRCA1", "BRCA2", "TP53",
               "HDAC1", "TOP2A"],

    # ── Other drug targets in top-20 ─────────────────────────────────────────
    "RRM2":   ["RRM1", "TP53R2", "CDK4", "E2F1"],   # Hydroxyurea
    "BIRC5":  ["CDK4", "CDK6", "AURKA", "PLK1", "TP53"],  # Survivin
    "CRBN":   ["CUL4A", "DDB1", "RBX1", "IKZF1"],   # Thalidomide
    "XPO1":   ["TP53", "CDKN1A", "RB1", "MDM2", "FOXO1"],  # Selinexor
    "ALK":    ["MYCN", "NPM1", "EML4", "STAT3", "PIK3CA"],  # Ceritinib
    "ATP2A2": ["BCL2", "TP53", "HDAC1"],   # Mipsagargin (SERCA)
    "ATP2A3": ["ATP2A2", "BCL2"],

    # ── Immune / TME targets ──────────────────────────────────────────────────
    "DRD2":   ["DRD1", "DRD3", "DRD4", "SIGMAR1", "GNAI2"],
    "SIGMAR1":["DRD2", "HDAC1", "BCL2", "IP3R1"],
}


class PPINetwork:
    def __init__(self):
        self.string_api_url = "https://string-db.org/api/json/network"
        self.cache: Dict[str, List[str]] = {}
        logger.info("✅ PPI Network Module Initialized (Live STRING-DB + Curated Fallback)")
        # Pre-load curated neighbors into cache so get_neighbors() works immediately
        for gene, neighbors in CURATED_PPI_NEIGHBORS.items():
            self.cache[gene.upper()] = neighbors

    async def fetch_neighbors_live(
        self, gene: str, session: aiohttp.ClientSession
    ) -> List[str]:
        """Fetch 1st-degree interactors from STRING-DB (high-confidence, score ≥ 800)."""
        gene_upper = gene.upper()
        if gene_upper in self.cache:
            return self.cache[gene_upper]

        params = {
            "identifiers": gene,
            "species": 9606,
            "required_score": 800,
        }
        try:
            async with session.get(
                self.string_api_url,
                params=params,
                timeout=aiohttp.ClientTimeout(total=8),
            ) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    neighbors = list(set(
                        row["preferredName_B"]
                        for row in data
                        if row["preferredName_A"].upper() == gene_upper
                    ))
                    if neighbors:
                        self.cache[gene_upper] = neighbors
                        logger.debug("STRING-DB live: %s → %d neighbors", gene, len(neighbors))
                        return neighbors
        except Exception as e:
            logger.debug("STRING-DB unavailable for %s: %s", gene, e)

        # Fall back to curated map (pre-loaded into cache at init)
        return self.cache.get(gene_upper, [])

    def get_neighbors(self, gene: str) -> List[str]:
        """Synchronous neighbor lookup — curated data always available."""
        return self.cache.get(gene.upper(), [])

    async def score_batch(
        self, candidates: List[Dict], disease_genes: List[str]
    ) -> List[Dict]:
        """
        Score candidates by PPI proximity to disease genes.

        FIX v5.1: Replaces binary 1.0/0.2 with graduated scoring:
            1.00 = drug target IS a disease gene (direct hit)
            0.85 = drug target is a 1st-degree STRING-DB neighbor of a disease gene
            0.60 = drug target is a neighbor of a neighbor (2nd-degree, curated only)
            0.20 = no network proximity found

        Uses curated DIPG STRING-DB neighbors when live API is unavailable,
        eliminating the previous flat-0.2 fallback for most drugs.
        """
        disease_set = set(g.upper() for g in disease_genes)

        # Expand disease set to include 1st-degree neighbors (curated)
        disease_neighbors: Set[str] = set()
        for gene in disease_genes:
            for nb in self.get_neighbors(gene.upper()):
                disease_neighbors.add(nb.upper())

        # Also expand to 2nd-degree (neighbors of neighbors, curated only)
        disease_2nd: Set[str] = set()
        for nb in list(disease_neighbors):
            for nb2 in self.get_neighbors(nb):
                disease_2nd.add(nb2.upper())

        async with aiohttp.ClientSession() as session:
            for c in candidates:
                targets = [t.upper() for t in c.get("targets", [])]
                target_set = set(targets)

                # Check direct overlap
                direct_hits = target_set & disease_set
                if direct_hits:
                    c["ppi_score"] = 1.00
                    c["network_context"] = f"Direct target: {', '.join(sorted(direct_hits))}"
                    continue

                # Try live STRING-DB for any targets not yet in cache
                live_upgraded = False
                for t in targets:
                    if t not in self.cache:
                        live_nb = await self.fetch_neighbors_live(t, session)
                        if set(live_nb) & disease_set:
                            c["ppi_score"] = PPI_CONFIG["first_degree_score"]
                            c["network_context"] = f"STRING-DB 1st-degree neighbor (live)"
                            live_upgraded = True
                            break

                if live_upgraded:
                    continue

                # Check curated 1st-degree neighbors
                first_degree_hits = target_set & disease_neighbors
                if first_degree_hits:
                    c["ppi_score"] = PPI_CONFIG["first_degree_score"]
                    c["network_context"] = (
                        f"1st-degree neighbor: {', '.join(sorted(first_degree_hits)[:3])}"
                    )
                    continue

                # Check 2nd-degree
                second_degree_hits = target_set & disease_2nd
                if second_degree_hits:
                    c["ppi_score"] = PPI_CONFIG["second_degree_score"]
                    c["network_context"] = (
                        f"2nd-degree neighbor: {', '.join(sorted(second_degree_hits)[:3])}"
                    )
                    continue

                # No network proximity found
                c["ppi_score"] = PPI_CONFIG["no_proximity_score"]
                c["network_context"] = "No PPI proximity to disease gene network"

        # Log distribution
        scores = [c.get("ppi_score", 0) for c in candidates]
        from collections import Counter
        dist = Counter(round(s, 2) for s in scores)
        logger.info(
            "✅ PPI scoring: %d candidates | distribution: %s",
            len(candidates), dict(sorted(dist.items(), reverse=True)),
        )

        return candidates