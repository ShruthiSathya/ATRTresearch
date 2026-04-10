"""
depmap_essentiality.py
======================
ATRT-specific DepMap CRISPR essentiality scoring.

Filters DepMap for ATRT/rhabdoid cell lines — all SMARCB1-null, directly
relevant to SMARCB1-loss dependency scoring.

KEY ATRT CELL LINES (verify against your Model.csv):
  BT16    : CNS ATRT, SMARCB1-null, ATRT-MYC subgroup — most studied
  BT37    : CNS ATRT, SMARCB1-null, ATRT-TYR subgroup
  BT12    : CNS ATRT, SMARCB1-null
  G401    : Renal rhabdoid, SMARCB1-null — same biology, commonly used
  A204    : Rhabdomyosarcoma/rhabdoid, SMARCB1-null
  MON     : Rhabdoid, SMARCB1-null
  CHLA02/04/06 : COG ATRT lines, SMARCB1-null
  KP-MRT-NS/RY : Mouse rhabdoid — may not be in human DepMap

WHY G401 AND A204 ARE INCLUDED:
  G401 and A204 are renal rhabdoid tumors (non-CNS) but share the same
  defining molecular feature: biallelic SMARCB1 loss. The SMARCB1-EZH2
  synthetic lethality was originally validated in G401/A204 (Knutson 2013
  PNAS). Drug sensitivities in these lines predict CNS ATRT responses.

CHRONOS INTERPRETATION:
  Chronos score ≤ -1.0 : strongly essential (best signal for targeting)
  Chronos -0.5 to -1.0 : moderately essential
  Chronos 0 to -0.5    : weakly essential
  Chronos ≥ 0          : non-essential / anti-essential

REFERENCES
-----------
Knutson 2013 : PNAS 110(19):7922. PMID 23620515. [EZH2 in G401/A204]
Tsherniak 2017: Cell 170(3):564. PMID 28753430. [DepMap methods]
Behan 2019   : Nature 568:511. PMID 30971826. [CRISPR essentiality]
DepMap Public 23Q4: https://depmap.org/portal/download/all/
"""

import asyncio
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List

try:
    from .pipeline_config import (
        PATHS, ATRT_SUBTYPE_TERMS, ATRT_LINEAGE_TERMS,
        ATRT_CELL_LINE_NAMES, MIN_LINES_SUBTYPE,
    )
except ImportError:
    from pipeline_config import (
        PATHS, ATRT_SUBTYPE_TERMS, ATRT_LINEAGE_TERMS,
        ATRT_CELL_LINE_NAMES, MIN_LINES_SUBTYPE,
    )

logger = logging.getLogger(__name__)


class DepMapEssentiality:
    """
    ATRT DepMap CRISPR essentiality scorer.

    Filters CRISPRGeneEffect.csv for SMARCB1-null rhabdoid/ATRT cell lines
    and computes per-gene Chronos essentiality scores.

    Three-tier matching strategy:
      1. OncotreeSubtype filter (MRT, ATRT, rhabdoid) — most precise
      2. Known cell line name filter — catches lines with generic subtypes
      3. Lineage fallback (brain/CNS/pediatric) — if 1+2 yield too few lines

    All three tiers are merged (union) before CRISPR data loading.
    """

    def __init__(self, data_dir: str = None):
        self.effect_file = Path(PATHS["depmap_effect"])
        self.model_file  = Path(PATHS["depmap_model"])
        self.is_ready    = False
        self.gene_scores: Dict[str, float] = {}
        self._loaded_cell_lines: List[str] = []
        logger.info("DepMap ATRT Essentiality Module initialized")

    async def _load_data_if_needed(self, disease: str = "atrt"):
        if self.is_ready:
            return
        if not self.effect_file.exists() or not self.model_file.exists():
            logger.warning(
                "DepMap CSVs not found.\n"
                "Download from: https://depmap.org/portal/download/all/\n"
                "  CRISPRGeneEffect.csv\n"
                "  Model.csv\n"
                "Place in: data/depmap/"
            )
            return

        logger.info("Loading DepMap CRISPR dataset for ATRT/rhabdoid lines...")

        models = pd.read_csv(self.model_file)
        models.columns = models.columns.str.strip()

        # Find column names in Model.csv
        subtype_col  = next((c for c in models.columns
                             if "oncotreesubtype" in c.lower() or c.lower() == "subtype"), None)
        lineage_col  = next((c for c in models.columns
                             if "lineage" in c.lower()), None)
        name_col     = next((c for c in models.columns
                             if c.lower() in ["celllinename", "cell_line_name",
                                              "stripped_cell_line_name", "displayname"]), None)
        model_id_col = next((c for c in models.columns
                             if c in ["ModelID", "DepMap_ID", "model_id"]), None)

        if model_id_col is None:
            logger.error(
                "Cannot find ModelID column in Model.csv.\n"
                "Available columns: %s", models.columns.tolist()
            )
            return

        atrt_cell_lines: List[str] = []

        # ── Tier 1: OncotreeSubtype filter ─────────────────────────────────────
        if subtype_col:
            subtype_mask = models[subtype_col].astype(str).str.lower().apply(
                lambda s: any(term in s for term in ATRT_SUBTYPE_TERMS)
            )
            tier1_lines = models[subtype_mask][model_id_col].tolist()
            logger.info(
                "OncotreeSubtype filter ('%s'): matched %d cell lines",
                subtype_col, len(tier1_lines)
            )
            atrt_cell_lines = list(set(atrt_cell_lines) | set(tier1_lines))

        # ── Tier 2: Known ATRT/rhabdoid cell line names ─────────────────────────
        if name_col:
            known_upper = {n.upper() for n in ATRT_CELL_LINE_NAMES}
            name_mask = models[name_col].astype(str).str.upper().str.strip().isin(known_upper)
            tier2_lines = models[name_mask][model_id_col].tolist()
            if tier2_lines:
                logger.info(
                    "Known cell line name filter ('%s'): matched %d additional lines",
                    name_col, len(set(tier2_lines) - set(atrt_cell_lines))
                )
                atrt_cell_lines = list(set(atrt_cell_lines) | set(tier2_lines))
        else:
            # Try OncotreeCode or other identifier columns
            for col in models.columns:
                if "code" in col.lower() or "id" in col.lower():
                    known_upper = {n.upper() for n in ATRT_CELL_LINE_NAMES}
                    mask = models[col].astype(str).str.upper().str.strip().isin(known_upper)
                    if mask.sum() > 0:
                        extra = models[mask][model_id_col].tolist()
                        atrt_cell_lines = list(set(atrt_cell_lines) | set(extra))
                        logger.info(
                            "Known line name via '%s': matched %d lines", col, mask.sum()
                        )
                        break

        # ── Tier 3: Lineage fallback if too few lines ───────────────────────────
        if len(atrt_cell_lines) < MIN_LINES_SUBTYPE and lineage_col:
            lineage_mask = models[lineage_col].astype(str).str.lower().apply(
                lambda s: any(term in s for term in ATRT_LINEAGE_TERMS)
            )
            lineage_lines = models[lineage_mask][model_id_col].tolist()
            added = set(lineage_lines) - set(atrt_cell_lines)
            if added:
                logger.info(
                    "Lineage fallback filter ('%s'): adding %d lines (brain/CNS/pediatric)",
                    lineage_col, len(added)
                )
                atrt_cell_lines = list(set(atrt_cell_lines) | set(lineage_lines))

        if not atrt_cell_lines:
            logger.warning(
                "No ATRT/rhabdoid cell lines found.\n"
                "Expected OncotreeSubtype values: MRT, ATRT, rhabdoid\n"
                "Known cell line names searched: %s\n"
                "Sample subtype values in Model.csv: %s",
                ATRT_CELL_LINE_NAMES[:5],
                models[subtype_col].dropna().unique()[:10].tolist() if subtype_col else "N/A"
            )
            return

        logger.info(
            "Total ATRT/rhabdoid cell lines identified: %d", len(atrt_cell_lines)
        )

        # ── Load CRISPR Gene Effect scores ─────────────────────────────────────
        logger.info(
            "Loading CRISPRGeneEffect.csv for %d lines (may take 30–60s)...",
            len(atrt_cell_lines)
        )
        effect_df = pd.read_csv(self.effect_file, index_col=0)

        # Intersect with available cell lines
        available = set(effect_df.index)
        matched   = [cl for cl in atrt_cell_lines if cl in available]

        logger.info(
            "Cell line matching: %d identified → %d present in CRISPRGeneEffect.csv",
            len(atrt_cell_lines), len(matched)
        )

        if not matched:
            logger.warning(
                "No overlap between Model.csv and CRISPRGeneEffect.csv IDs.\n"
                "Model.csv sample IDs: %s\n"
                "CRISPRGeneEffect.csv sample IDs: %s",
                atrt_cell_lines[:5], list(effect_df.index[:5])
            )
            return

        relevant_effects = effect_df.loc[matched]

        # Parse column format: "EGFR (1956)" → "EGFR"
        for col in relevant_effects.columns:
            gene_symbol = col.split(" ")[0].upper().strip()
            self.gene_scores[gene_symbol] = float(
                relevant_effects[col].median()
            )

        self._loaded_cell_lines = matched
        self.is_ready           = True

        # Log essential genes found (Chronos ≤ -1.0)
        essential = {g: s for g, s in self.gene_scores.items() if s <= -1.0}
        top_essential = sorted(essential.items(), key=lambda x: x[1])[:15]

        logger.info(
            "DepMap loaded: %d ATRT/rhabdoid cell lines, %d gene scores\n"
            "  Cell lines: %s\n"
            "  Top essential genes (Chronos ≤ -1.0): %d total\n"
            "  Examples: %s",
            len(matched), len(self.gene_scores),
            matched,
            len(essential),
            [(g, round(s, 2)) for g, s in top_essential[:8]]
        )

        # Validate key ATRT targets are present
        key_atrt_targets = {
            "EZH2", "BRD4", "HDAC1", "AURKA", "CDK4",
            "PSMB5", "MYC", "MYCN",
        }
        found_key = key_atrt_targets & set(self.gene_scores.keys())
        missing_key = key_atrt_targets - found_key
        if missing_key:
            logger.warning(
                "Key ATRT targets missing from DepMap gene scores: %s\n"
                "These will receive neutral prior (0.50).",
                sorted(missing_key)
            )
        else:
            logger.info(
                "All key ATRT targets found in DepMap: %s",
                {g: round(self.gene_scores[g], 2) for g in sorted(found_key)}
            )

    async def score_batch(
        self, candidates: List[Dict], disease: str = "atrt"
    ) -> List[Dict]:
        await self._load_data_if_needed(disease)

        for c in candidates:
            if not self.is_ready:
                c["depmap_score"] = 0.50
                c["depmap_note"]  = "DepMap data not loaded — neutral prior 0.50"
                continue

            targets       = c.get("targets", [])
            target_scores = [
                self.gene_scores.get(t.upper(), 0.0) for t in targets
            ]

            if target_scores:
                best_chronos = min(target_scores)   # lower = more lethal

                if best_chronos <= -1.0:
                    c["depmap_score"] = 1.00
                    c["depmap_note"]  = (
                        f"Essential in ATRT/rhabdoid lines "
                        f"(Chronos {best_chronos:.2f} ≤ -1.0)"
                    )
                elif best_chronos <= -0.5:
                    c["depmap_score"] = 0.80
                    c["depmap_note"]  = (
                        f"Strongly essential (Chronos {best_chronos:.2f})"
                    )
                elif best_chronos < 0:
                    c["depmap_score"] = 0.50
                    c["depmap_note"]  = (
                        f"Modestly essential (Chronos {best_chronos:.2f})"
                    )
                else:
                    c["depmap_score"] = 0.10
                    c["depmap_note"]  = (
                        f"Non-essential in ATRT lines "
                        f"(Chronos {best_chronos:.2f} ≥ 0)"
                    )
            else:
                c["depmap_score"] = 0.10
                c["depmap_note"]  = "No targets matched DepMap gene list"

        return candidates

    def get_coverage_report(self) -> str:
        if not self.is_ready:
            return (
                "DepMap: not loaded (CRISPRGeneEffect.csv / Model.csv missing "
                "or no ATRT/rhabdoid cell lines matched)"
            )
        essential_count = sum(1 for s in self.gene_scores.values() if s <= -1.0)
        return (
            f"DepMap ATRT: {len(self._loaded_cell_lines)} rhabdoid/ATRT cell lines "
            f"| {len(self.gene_scores)} gene Chronos scores "
            f"| {essential_count} essential (Chronos ≤ -1.0)\n"
            f"Lines: {self._loaded_cell_lines}"
        )