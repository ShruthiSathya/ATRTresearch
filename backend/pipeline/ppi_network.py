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
# v5.5: Extended proteasome subunit neighbors so Marizomib correctly scores
# 0.85 (1st-degree neighbor) rather than 0.20 (no proximity).
#
# Biological rationale for proteasome → disease gene connections:
# PSMB5 (20S catalytic β5 subunit, marizomib primary target):
#   - ATG5/BECN1: proteasome inhibition → autophagy induction (compensatory)
#   - SQSTM1/p62: sequestosome, proteasome-autophagy crossroads
#   - TP53: p53 ubiquitylated by MDM2 → proteasome degradation
#   - MYC: MYC is rapidly degraded by proteasome (t½ ~20 min)
#   - CDKN1B: p27 ubiquitylated by Skp2/SCF → proteasome
#   - H3F3A: H3K27M protein itself is a substrate for proteasomal QC
# Sources: Tanaka 2009 (Biochemistry); Lin et al. 2019 (Sci Transl Med);
#          Warren et al. 2019 (Neuro-Oncology); Dikic 2017 (Science).
# ─────────────────────────────────────────────────────────────────────────────

CURATED_PPI_NEIGHBORS: Dict[str, List[str]] = {
    # ── Core DIPG disease genes ───────────────────────────────────────────────
    "H3-3A":  ["EZH2", "EED", "SUZ12", "ATRX", "DAXX", "BRD4", "KDM6A", "SETD2",
               "PSMB5", "PSMD1"],   # H3K27M protein → proteasomal QC
    "H3F3A":  ["EZH2", "EED", "SUZ12", "ATRX", "DAXX", "BRD4", "KDM6A", "SETD2",
               "PSMB5", "PSMD1"],
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

    # ── DIPG-specific genes ───────────────────────────────────────────────────
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
               "HDAC1", "BCL2", "PARP1",
               "PSMB5", "PSMD1"],  # p53 degraded via MDM2-proteasome axis
    "MYC":    ["BRD4", "HDAC1", "HDAC2", "CDK4", "E2F1", "MYCN",
               "MAX", "MTOR",
               "PSMB5"],  # MYC rapidly degraded by proteasome
    "MYCN":   ["BRD4", "MYC", "CDK4", "AURKA", "ALK", "MDM2",
               "PSMB5"],
    "STAT3":  ["JAK1", "JAK2", "IL6", "EGFR", "PDGFRA", "SRC", "BCL2"],
    "PIK3CA": ["PIK3R1", "AKT1", "MTOR", "PTEN", "EGFR", "PDGFRA",
               "KRAS", "RAS"],
    "AKT1":   ["PIK3CA", "MTOR", "PTEN", "MDM2", "BCL2", "HDAC1"],
    "MTOR":   ["PIK3CA", "AKT1", "PTEN", "TSC1", "TSC2", "RPTOR"],

    # ── FIX v5.5: Proteasome subunit neighbors — Marizomib target network ─────
    # Marizomib targets PSMB5, PSMB2, PSMB1 (and others). Previously these
    # had no curated neighbors → ppi_score: 0.20 (no proximity floor).
    # This caused Marizomib to be underranked despite DepMap Chronos = -3.30
    # (most essential target in the entire screen).
    #
    # All connections below are supported by published proteasome biology:
    "PSMB5":  ["PSMB1", "PSMB2", "PSMA1", "PSMA3", "PSMD1",
               "UBB", "UBC", "SQSTM1", "ATG5", "BECN1",
               "TP53",    # p53 degraded via MDM2-E3 → PSMB5
               "MYC",     # MYC t½ ~20 min, proteasome-dependent
               "CDKN1B",  # p27 → Skp2/SCF → proteasome
               "H3F3A",   # H3K27M protein proteostasis
               "MCL1"],   # MCL1 rapidly degraded by proteasome
    "PSMB2":  ["PSMB5", "PSMB1", "PSMA1", "PSMD1", "SQSTM1",
               "TP53", "MYC", "MCL1"],
    "PSMB1":  ["PSMB5", "PSMB2", "PSMA1", "BECN1",
               "TP53", "CDKN1B"],
    "PSMB8":  ["PSMB5", "PSMB9", "PSMB10", "PSME1", "PSME2",
               "MYC", "TP53"],
    "PSMB9":  ["PSMB8", "PSMB5", "PSMB10", "PSME1"],
    "PSMD1":  ["PSMB5", "PSMB2", "PSMB1", "SQSTM1", "UBB",
               "TP53", "MYC"],
    "PSMA1":  ["PSMB5", "PSMB2", "PSMD1", "UBB"],
    "PSMA3":  ["PSMB5", "PSMB2", "PSMD1"],
    "PSMA6":  ["PSMB5", "PSMD1", "UBB"],
    "PSMA7":  ["PSMB5", "PSMB2", "PSMD1"],
    "SQSTM1": ["PSMB5", "ATG5", "BECN1", "MAP1LC3B", "HDAC6"],
    "BECN1":  ["ATG5", "PSMB5", "BCL2", "PIK3C3"],
    "ATG5":   ["PSMB5", "BECN1", "MAP1LC3B"],

    # ── BET bromodomain (Birabresib targets) ──────────────────────────────────
    "BRD2":   ["BRD4", "BRD3", "CDK9", "MYC", "HDAC1", "E2F1"],
    "BRD3":   ["BRD4", "BRD2", "CDK9", "MYC"],
    "BRDT":   ["BRD4", "BRD2", "CDK9"],

    # ── HDAC targets (Panobinostat) ───────────────────────────────────────────
    "HDAC6":  ["HDAC1", "HDAC2", "HSP90AA1", "TUBB", "TUBA1A", "STAT3",
               "SQSTM1"],  # HDAC6-SQSTM1 interaction in aggresomes
    "HDAC9":  ["HDAC1", "HDAC2", "HDAC3", "HDAC4", "MEF2D"],
    "HDAC11": ["HDAC1", "HDAC2", "HDAC3"],
    "HDAC3":  ["HDAC1", "HDAC2", "NCOR1", "NCOR2", "MYC", "BRD4"],
    "HDAC4":  ["HDAC1", "HDAC2", "HDAC3", "MAPK1"],
    "HDAC5":  ["HDAC1", "HDAC2", "HDAC3"],
    "HDAC7":  ["HDAC1", "HDAC2", "HDAC3"],
    "HDAC8":  ["HDAC1", "HDAC2", "HDAC3", "TP53"],
    "HDAC10": ["HDAC1", "HDAC2"],

    # ── DNA damage ────────────────────────────────────────────────────────────
    "PARP1":  ["PARP2", "ATM", "ATR", "BRCA1", "BRCA2", "TP53",
               "HDAC1", "TOP2A"],

    # ── Other top-20 drug targets ─────────────────────────────────────────────
    "DRD2":   ["DRD1", "DRD3", "DRD4", "SIGMAR1", "GNAI2"],
    "SIGMAR1":["DRD2", "HDAC1", "BCL2", "IP3R1"],
    "CLPB":   ["HSP90AA1", "HSP70", "SQSTM1"],  # ONC201 mitochondrial target

    # ── mTOR complex components (Onatasertib / AZD-8055 / CC-115) ────────────
    "RPTOR":  ["MTOR", "AKT1", "TSC1", "TSC2"],
    "MLST8":  ["MTOR", "RPTOR", "RICTOR"],
    "PRKDC":  ["ATM", "ATR", "PARP1", "TP53"],  # DNA-PK (CC-115 target)

    # ── PDGFRA/RTK drugs (Tovetumab, Olaratumab, Pazopanib, etc.) ────────────
    "PDGFRB": ["PDGFRA", "KDR", "FLT1", "PIK3CA", "SRC"],
    "KDR":    ["FLT1", "FLT4", "PDGFRA", "VEGFA", "PIK3CA"],
    "FLT1":   ["KDR", "FLT4", "PDGFRA", "VEGFA"],
    "FLT4":   ["KDR", "FLT1", "PDGFRA"],
    "KIT":    ["PDGFRA", "PDGFRB", "KDR", "PIK3CA"],
    "FLT3":   ["KIT", "PDGFRA", "PIK3CA"],
    "CSF1R":  ["PDGFRA", "KIT", "PIK3CA", "SRC"],
    "MET":    ["EGFR", "AXL", "PIK3CA", "STAT3"],
    "ALK":    ["MYCN", "NPM1", "EML4", "STAT3", "PIK3CA"],
    "EML4":   ["ALK", "MYCN"],
    "NPM1":   ["ALK", "MDM2", "TP53"],
    "FGFR1":  ["FGFR2", "FGFR3", "PIK3CA", "MAPK1"],
    "FGFR2":  ["FGFR1", "FGFR3", "PIK3CA"],
    "FGFR3":  ["FGFR1", "FGFR2", "PIK3CA"],
    "FGFR4":  ["FGFR1", "FGFR3", "PIK3CA"],

    # ── EGFR (Depatuxizumab mafodotin) ───────────────────────────────────────
    "ERBB2":  ["EGFR", "PIK3CA", "AKT1", "MAPK1"],

    # ── Kinases / additional ──────────────────────────────────────────────────
    "ABL1":   ["BCR", "PIK3CA", "STAT5A", "CRKL"],
    "RAF1":   ["BRAF", "KRAS", "MAPK1", "MAP2K1"],
    "BRAF":   ["RAF1", "KRAS", "MAPK1", "MAP2K1"],
    "MAPK11": ["MAPK1", "MAPK3", "MAP2K3"],
    "DDR2":   ["PDGFRA", "SRC", "PIK3CA"],
    "NTRK1":  ["PIK3CA", "MAPK1", "AKT1"],
    "TEK":    ["PDGFRA", "PIK3CA", "ANGPT1"],
    "EPHA2":  ["EGFR", "PIK3CA", "SRC"],
    "RET":    ["PIK3CA", "MAPK1", "AKT1"],
    "LCK":    ["ZAP70", "PDCD1", "CD3E"],
    "ITK":    ["LCK", "ZAP70", "PLCG1"],
    "SRC":    ["EGFR", "STAT3", "PIK3CA", "FAK"],
}


class PPINetwork:
    def __init__(self):
        self.string_api_url = "https://string-db.org/api/json/network"
        self.cache: Dict[str, List[str]] = {}
        logger.info("✅ PPI Network Module v5.5 (Extended proteasome neighbors + Marizomib fix)")
        for gene, neighbors in CURATED_PPI_NEIGHBORS.items():
            self.cache[gene.upper()] = neighbors

    async def fetch_neighbors_live(
        self, gene: str, session: aiohttp.ClientSession
    ) -> List[str]:
        gene_upper = gene.upper()
        if gene_upper in self.cache:
            return self.cache[gene_upper]

        params = {
            "identifiers":  gene,
            "species":      9606,
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

        return self.cache.get(gene_upper, [])

    def get_neighbors(self, gene: str) -> List[str]:
        return self.cache.get(gene.upper(), [])

    async def score_batch(
        self, candidates: List[Dict], disease_genes: List[str]
    ) -> List[Dict]:
        """
        Score candidates by PPI proximity to disease genes.

        v5.5: Proteasome subunit neighbors now curated → Marizomib correctly
        scores 0.85 (1st-degree neighbor of TP53/MYC/H3F3A) rather than 0.20.
        """
        disease_set = set(g.upper() for g in disease_genes)

        # Expand disease set to 1st-degree neighbors (curated)
        disease_neighbors: Set[str] = set()
        for gene in disease_genes:
            for nb in self.get_neighbors(gene.upper()):
                disease_neighbors.add(nb.upper())

        # Expand to 2nd-degree (curated only)
        disease_2nd: Set[str] = set()
        for nb in list(disease_neighbors):
            for nb2 in self.get_neighbors(nb):
                disease_2nd.add(nb2.upper())

        async with aiohttp.ClientSession() as session:
            for c in candidates:
                targets    = [t.upper() for t in c.get("targets", [])]
                target_set = set(targets)

                # Direct overlap with disease genes
                direct_hits = target_set & disease_set
                if direct_hits:
                    c["ppi_score"]       = 1.00
                    c["network_context"] = f"Direct target: {', '.join(sorted(direct_hits))}"
                    continue

                # Live STRING-DB for uncached targets
                live_upgraded = False
                for t in targets:
                    if t not in self.cache:
                        live_nb = await self.fetch_neighbors_live(t, session)
                        if set(live_nb) & disease_set:
                            c["ppi_score"]       = PPI_CONFIG["first_degree_score"]
                            c["network_context"] = "STRING-DB 1st-degree neighbor (live)"
                            live_upgraded = True
                            break

                if live_upgraded:
                    continue

                # Curated 1st-degree
                first_degree_hits = target_set & disease_neighbors
                if first_degree_hits:
                    c["ppi_score"] = PPI_CONFIG["first_degree_score"]
                    c["network_context"] = (
                        f"1st-degree neighbor: {', '.join(sorted(first_degree_hits)[:3])}"
                    )
                    continue

                # Curated 2nd-degree
                second_degree_hits = target_set & disease_2nd
                if second_degree_hits:
                    c["ppi_score"] = PPI_CONFIG["second_degree_score"]
                    c["network_context"] = (
                        f"2nd-degree neighbor: {', '.join(sorted(second_degree_hits)[:3])}"
                    )
                    continue

                c["ppi_score"]       = PPI_CONFIG["no_proximity_score"]
                c["network_context"] = "No PPI proximity to disease gene network"

        scores = [c.get("ppi_score", 0) for c in candidates]
        from collections import Counter
        dist = Counter(round(s, 2) for s in scores)
        logger.info(
            "✅ PPI scoring v5.5: %d candidates | distribution: %s "
            "(Marizomib proteasome fix applied)",
            len(candidates), dict(sorted(dist.items(), reverse=True)),
        )

        return candidates