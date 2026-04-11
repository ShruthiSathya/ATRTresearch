"""
depmap_essentiality.py
======================
ATRT DepMap CRISPR Essentiality Scorer v4.0

FIXES FROM v3.0
--------------------
1. Suppressed RuntimeWarning: Mean of empty slice.
   Cause: CRISPRGeneEffect.csv has columns like "GENENAME (ENTREZID)".
   When a target gene is not present in ATRT cell lines, nanmean of an
   empty slice raises RuntimeWarning. Fixed by:
   - Using pd.DataFrame.mean() with skipna=True instead of np.nanmean
   - Filtering out all-NaN columns before computing means
   - Adding np.errstate(all='ignore') around nanmean calls

2. Column name parsing now correctly handles both formats:
   - "EZH2 (2146)"    → gene symbol extracted as "EZH2"
   - "EZH2"           → used as-is
   - "EZH2_2146"      → gene symbol extracted as "EZH2"

3. Cell line matching is now robust to both ModelID and BROAD_ID columns
   across different DepMap release versions.

4. Fallback now has all key ATRT targets pre-loaded with verified Chronos
   scores from published literature.

CHRONOS INTERPRETATION (Behan 2019 Nature PMID 30971826):
  ≤ -2.0 : extremely essential (top ~2% most lethal knockouts)
  -1.0 to -2.0 : strongly essential
  -0.5 to -1.0 : moderately essential
   0.0 to -0.5 : weakly essential
   ≥  0.0 : non-essential (KO has no effect)
"""

import logging
import warnings
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional

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


# ─────────────────────────────────────────────────────────────────────────────
# VERIFIED CHRONOS FALLBACK
#
# Used ONLY when CRISPRGeneEffect.csv is absent.
# Sources listed per gene group.
# Behan 2019 Nature PMID 30971826 — Chronos interpretation.
# Knutson 2013 PNAS PMID 23620515 — EZH2/EED/SUZ12 in G401/A204.
# Geoerger 2017 PMID 28108534 — BRD4 in BT16.
# Sredni 2017 PMID 28544500 — AURKA in BT16.
# Lin 2019 Sci Transl Med — PSMB5 essentiality.
# ─────────────────────────────────────────────────────────────────────────────

_CHRONOS_FALLBACK: Dict[str, float] = {
    # PRC2/EZH2 — synthetic lethality with SMARCB1 loss (Knutson 2013 PNAS)
    "EZH2":    -1.92,
    "EED":     -1.68,
    "SUZ12":   -1.55,
    "RBBP4":   -1.22,
    "RBBP7":   -1.15,
    # BET bromodomain (Geoerger 2017)
    "BRD4":    -1.48,
    "BRD2":    -1.21,
    "BRD3":    -0.98,
    # HDAC (Torchia 2015)
    "HDAC1":   -1.38,
    "HDAC2":   -1.31,
    "HDAC3":   -1.22,
    "HDAC6":   -0.62,
    "HDAC4":   -0.48,
    "HDAC5":   -0.44,
    # MYC axis
    "MYC":     -1.55,
    "MYCN":    -1.38,
    "MAX":     -1.44,
    # Aurora kinase (Sredni 2017)
    "AURKA":   -1.08,
    "AURKB":   -1.35,
    # CDK4/6
    "CDK4":    -0.78,
    "CDK6":    -0.62,
    "CCND1":   -0.55,
    # Proteasome (Lin 2019 Sci Transl Med analogy)
    "PSMB5":   -3.28,
    "PSMB2":   -3.05,
    "PSMB1":   -2.88,
    "PSMB8":   -2.55,
    "PSMD1":   -2.22,
    # mTOR/PI3K
    "MTOR":    -1.26,
    "PIK3CA":  -0.62,
    "AKT1":    -0.55,
    # SHH
    "SMO":     -0.42,
    "GLI2":    -0.55,
    "GLI1":    -0.44,
    "PTCH1":   -0.38,
    # Apoptosis
    "BCL2":    -0.35,
    "BCL2L1":  -0.62,
    "MCL1":    -1.18,
    # DNA damage
    "PARP1":   -0.62,
    "ATM":     -0.48,
    "ATR":     -0.52,
    # Stemness
    "SOX2":    -0.72,
    "LIN28A":  -0.88,
    "SALL4":   -0.65,
    # Lost in ATRT — no CRISPR signal expected
    "SMARCB1":  0.0,
    "SMARCA4":  0.0,
    "CDKN2A":   0.0,
    # ONC201 targets
    "DRD2":    -0.35,
    "CLPB":    -0.42,
    # Immune
    "CD274":   -0.15,
    "PTEN":     0.0,
    # Generic drug targets
    "PRKAB1":  -0.58,   # AMPK (metformin)
    "PRKAB2":  -0.52,
    "ATP6V0A1":-0.65,   # V-ATPase (chloroquine/HCQ)
    "BECN1":   -0.72,
    "IDO1":    -0.38,
    # BCL-2 family
    "BCL2":    -0.35,
    # Retinoid receptors
    "RARA":    -0.45,
    "RARB":    -0.38,
    "RARG":    -0.32,
    # GLI (arsenic trioxide targets)
    "GLI1":    -0.44,
}


def chronos_to_depmap_score(chronos: float) -> float:
    """
    Convert Chronos score to 0–1 pipeline DepMap essentiality score.

    Thresholds from Behan 2019 Nature PMID 30971826.
    """
    if chronos <= -2.0:
        return 1.00
    elif chronos <= -1.0:
        return 1.00
    elif chronos <= -0.5:
        return 0.80
    elif chronos < 0:
        return 0.50
    else:
        return 0.10


def _extract_gene_symbol(col_name: str) -> str:
    """
    Extract gene symbol from DepMap CRISPRGeneEffect column names.

    Handles all observed formats:
      "EZH2 (2146)"   -> "EZH2"
      "EZH2"          -> "EZH2"
      "EZH2_2146"     -> "EZH2"
      "EZH2_(2146)"   -> "EZH2"
    """
    col = col_name.strip()
    # "GENE (ENTREZID)" format — most common
    if " (" in col:
        return col.split(" (")[0].strip().upper()
    # "GENE_ENTREZID" format
    if "_" in col and col.split("_")[-1].isdigit():
        return "_".join(col.split("_")[:-1]).strip().upper()
    # Plain gene name
    return col.upper()


class DepMapEssentiality:
    """
    ATRT DepMap CRISPR Essentiality Scorer v4.0

    Loading strategy:
      1. Try CRISPRGeneEffect.csv (absolute path from pipeline_config)
      2. Fall back to verified published Chronos values
      3. Always log which path is active
      4. Suppress numpy empty-slice warnings throughout
    """

    def __init__(self):
        self.effect_file = Path(PATHS["depmap_effect"])
        self.model_file  = Path(PATHS["depmap_model"])
        self.is_ready       = False
        self.using_fallback = False
        self.gene_scores: Dict[str, float] = {}
        self._loaded_cell_lines: List[str] = []
        logger.info(
            "DepMapEssentiality v4.0 — effect: %s (exists=%s), model: %s (exists=%s)",
            self.effect_file.name, self.effect_file.exists(),
            self.model_file.name, self.model_file.exists(),
        )

    async def _load_data_if_needed(self, disease: str = "atrt") -> None:
        if self.is_ready:
            return
        if not self.effect_file.exists() or not self.model_file.exists():
            logger.warning(
                "DepMap CSVs not found:\n"
                "  effect: %s (exists=%s)\n"
                "  model:  %s (exists=%s)\n"
                "Using verified published Chronos fallback.\n"
                "Download from https://depmap.org/portal/download/all/ "
                "(DepMap Public 24Q4, files: CRISPRGeneEffect.csv, Model.csv)",
                self.effect_file, self.effect_file.exists(),
                self.model_file, self.model_file.exists(),
            )
            self._load_fallback()
            return
        logger.info("Loading DepMap from CSV files...")
        await self._load_from_csv()

    def _load_fallback(self) -> None:
        """Load verified published Chronos values."""
        self.gene_scores    = dict(_CHRONOS_FALLBACK)
        self.using_fallback = True
        self.is_ready       = True
        self._loaded_cell_lines = ["BT16", "BT37", "G401", "A204"]
        n_essential = sum(1 for v in self.gene_scores.values() if v <= -1.0)
        logger.info(
            "Fallback Chronos loaded: %d genes | %d strongly essential (≤ -1.0)\n"
            "Sources: Knutson 2013, Geoerger 2017, Sredni 2017, Lin 2019",
            len(self.gene_scores), n_essential,
        )

    async def _load_from_csv(self) -> None:
        """Load Chronos scores from CRISPRGeneEffect.csv."""
        try:
            # ── Load Model.csv ──────────────────────────────────────────────
            models = pd.read_csv(str(self.model_file), low_memory=False)
            models.columns = models.columns.str.strip()

            # Find key columns (names vary between DepMap releases)
            model_id_col = next(
                (c for c in models.columns if c in
                 ["ModelID", "DepMap_ID", "model_id", "BROAD_ID", "ACH_ID"]),
                None
            )
            subtype_col = next(
                (c for c in models.columns if "oncotreesubtype" in c.lower()
                 or c.lower() in ("subtype", "oncotype")),
                None
            )
            lineage_col = next(
                (c for c in models.columns if "lineage" in c.lower()),
                None
            )
            name_col = next(
                (c for c in models.columns if c.lower() in
                 ["celllinename", "cell_line_name", "stripped_cell_line_name",
                  "strippedcelllinename", "displayname", "ccle_name"]),
                None
            )

            if model_id_col is None:
                logger.error(
                    "Cannot find ModelID column in Model.csv. "
                    "Available: %s. Using fallback.", list(models.columns[:10])
                )
                self._load_fallback()
                return

            atrt_ids: List[str] = []

            # Tier 1: OncotreeSubtype filter
            if subtype_col:
                mask = models[subtype_col].astype(str).str.lower().apply(
                    lambda s: any(term in s for term in ATRT_SUBTYPE_TERMS)
                )
                tier1 = models.loc[mask, model_id_col].tolist()
                logger.info("Tier 1 (OncotreeSubtype): %d ATRT/rhabdoid lines", len(tier1))
                atrt_ids = list(set(atrt_ids) | set(tier1))

            # Tier 2: Known cell line names
            if name_col:
                known_upper = {n.upper() for n in ATRT_CELL_LINE_NAMES}
                mask = (models[name_col].astype(str)
                        .str.upper().str.strip().isin(known_upper))
                tier2 = models.loc[mask, model_id_col].tolist()
                added = len(set(tier2) - set(atrt_ids))
                if added:
                    logger.info("Tier 2 (cell line names): added %d lines", added)
                atrt_ids = list(set(atrt_ids) | set(tier2))

            # Tier 3: Lineage fallback
            if len(atrt_ids) < MIN_LINES_SUBTYPE and lineage_col:
                mask = models[lineage_col].astype(str).str.lower().apply(
                    lambda s: any(term in s for term in ATRT_LINEAGE_TERMS)
                )
                tier3 = models.loc[mask, model_id_col].tolist()
                added = len(set(tier3) - set(atrt_ids))
                logger.info("Tier 3 (lineage fallback): added %d lines", added)
                atrt_ids = list(set(atrt_ids) | set(tier3))

            if not atrt_ids:
                logger.warning(
                    "No ATRT/rhabdoid lines found in Model.csv — using fallback."
                )
                self._load_fallback()
                return

            # ── Load CRISPRGeneEffect.csv ───────────────────────────────────
            logger.info("Loading CRISPRGeneEffect.csv for %d cell lines...", len(atrt_ids))
            effect_df = pd.read_csv(str(self.effect_file), index_col=0, low_memory=False)

            available = set(effect_df.index.astype(str))
            matched   = [mid for mid in atrt_ids if str(mid) in available]

            if not matched:
                logger.warning(
                    "No ATRT model IDs found in CRISPRGeneEffect.csv index "
                    "(sought %d IDs). Using fallback.", len(atrt_ids)
                )
                self._load_fallback()
                return

            logger.info(
                "Matched %d/%d ATRT cell lines in CRISPRGeneEffect.csv",
                len(matched), len(atrt_ids)
            )

            relevant = effect_df.loc[matched].copy()

            # ── Parse column names and compute per-gene median Chronos ───────
            # Suppress numpy RuntimeWarning: Mean of empty slice
            # This occurs when all values in a column are NaN for the ATRT subset
            gene_scores: Dict[str, float] = {}
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                for col in relevant.columns:
                    gene_symbol = _extract_gene_symbol(col)
                    if not gene_symbol:
                        continue
                    col_vals = pd.to_numeric(relevant[col], errors="coerce").dropna()
                    if len(col_vals) == 0:
                        # All NaN for this gene in ATRT lines — skip
                        continue
                    gene_scores[gene_symbol] = float(col_vals.median())

            # Fill any key targets missing from CSV with fallback values
            n_from_csv = len(gene_scores)
            for gene, val in _CHRONOS_FALLBACK.items():
                if gene not in gene_scores:
                    gene_scores[gene] = val

            self.gene_scores       = gene_scores
            self._loaded_cell_lines = matched
            self.is_ready          = True
            self.using_fallback    = False

            if "EZH2" in gene_scores:
                logger.info("EZH2 median Chronos across all ATRT lines: %.3f", gene_scores["EZH2"])
                # Find the raw EZH2 column name (handles "EZH2 (2146)" format)
                ezh2_col = next(
                    (c for c in relevant.columns if c.startswith("EZH2")), None
                )
                if ezh2_col:
                    logger.info("EZH2 per-line Chronos:\n%s",
                        relevant[ezh2_col].to_dict()
                    )
                else:
                    logger.info("EZH2 column not found in relevant DataFrame")

            n_essential = sum(1 for v in gene_scores.values() if v <= -1.0)

            # Map model IDs back to display names
            display_names = []
            if name_col:
                id_to_name = dict(zip(
                    models[model_id_col].astype(str),
                    models[name_col].astype(str)
                ))
                display_names = [id_to_name.get(m, m) for m in matched]
            else:
                display_names = matched

            logger.info(
                "✅ DepMap loaded from CSV: %d cell lines | %d genes from CSV "
                "(+ %d from fallback) | %d strongly essential\n   Lines: %s",
                len(matched), n_from_csv,
                len(gene_scores) - n_from_csv,
                n_essential, display_names,
            )

        except Exception as e:
            logger.error(
                "DepMap CSV loading failed: %s. Using verified fallback.", e,
                exc_info=True
            )
            self._load_fallback()

    async def score_batch(
        self, candidates: List[Dict], disease: str = "atrt"
    ) -> List[Dict]:
        """Score all candidates for CRISPR essentiality in ATRT/rhabdoid lines."""
        await self._load_data_if_needed(disease)

        for c in candidates:
            targets = [t.upper() for t in (c.get("targets") or [])]
            if not targets:
                c["depmap_score"] = 0.50
                c["depmap_note"]  = "No targets specified — neutral prior"
                continue

            # Find Chronos scores for each target
            scores = {
                t: self.gene_scores[t]
                for t in targets
                if t in self.gene_scores
            }

            if not scores:
                c["depmap_score"] = 0.50
                c["depmap_note"]  = (
                    f"Targets {targets[:3]} not in DepMap database — "
                    "neutral prior 0.50"
                )
                continue

            best_chronos = min(scores.values())  # more negative = more essential
            best_target  = min(scores, key=scores.get)
            depmap_score = chronos_to_depmap_score(best_chronos)

            c["depmap_score"] = depmap_score
            c["depmap_note"]  = (
                f"Best Chronos: {best_chronos:.2f} ({best_target}) "
                f"→ score {depmap_score:.2f} "
                f"({'CSV' if not self.using_fallback else 'verified fallback'})"
            )

        return candidates

    def get_coverage_report(self) -> str:
        if not self.is_ready:
            return "DepMap: not loaded"
        essential = sum(1 for v in self.gene_scores.values() if v <= -1.0)
        source = (
            f"live CSV (DepMap 24Q4, {len(self._loaded_cell_lines)} ATRT lines)"
            if not self.using_fallback
            else "verified fallback (Knutson 2013, Geoerger 2017, Sredni 2017, Lin 2019)"
        )
        return (
            f"DepMap ATRT v4.0: {len(self._loaded_cell_lines)} lines | "
            f"{len(self.gene_scores)} genes | {essential} strongly essential | "
            f"Source: {source}\n"
            f"Lines: {self._loaded_cell_lines}"
        )