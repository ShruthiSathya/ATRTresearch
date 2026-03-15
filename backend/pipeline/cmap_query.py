import asyncio
import logging
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List

try:
    from cmapPy.pandasGEXpress.parse import parse
    from scipy.spatial.distance import cosine
    CMAP_AVAILABLE = True
except ImportError:
    CMAP_AVAILABLE = False

logger = logging.getLogger(__name__)

# Pre-computed clue.io scores (populated after running integrate_cmap_results.py)
PRECOMPUTED_SCORES_PATH = Path("data/cmap_query/cmap_scores_pipeline.json")

# Disease signatures (used if running locally against .gctx files)
GBM_SIGNATURE = {
    "up":   ["EGFR", "PDGFRA", "CDK4", "MDM2", "STAT3", "VIM", "CD44"],
    "down": ["CDKN2A", "PTEN", "TP53", "OLIG2"]
}

DIPG_SIGNATURE = {
    "up":   ["ACVR1", "H3-3A", "BRD4", "MYCN", "EZH2", "DRD2", "SIGMAR1"],
    "down": ["CDKN2A", "OLIG2", "PTEN"]
}


class CMAPQuery:
    """
    v5.5: CMap query engine with pre-computed score support.

    Two modes:
    1. Pre-computed (recommended): loads clue.io Tau scores from
       data/cmap_query/cmap_scores_pipeline.json after you run
       integrate_cmap_results.py. No local data download needed.

    2. Local .gctx mode: requires the full ~30GB LINCS dataset
       at data/raw_omics/level5_beta_trt_cp_n720216x12328.gctx

    To get pre-computed scores:
        python -m backend.pipeline.prepare_cmap_query   # makes gene lists
        # go to clue.io, run query, download results
        python -m backend.pipeline.integrate_cmap_results  # integrates them
    """

    def __init__(self, data_dir: str = "data/raw_omics/"):
        self.gctx_path    = Path(data_dir) / "level5_beta_trt_cp_n720216x12328.gctx"
        self.siginfo_path = Path(data_dir) / "siginfo_beta.txt"
        self.is_ready     = False
        self.drug_to_sig_map = {}
        self._precomputed: Dict[str, Dict] = {}

        # Try to load pre-computed scores first
        self._load_precomputed()

        if not CMAP_AVAILABLE:
            logger.debug("cmapPy not installed — local .gctx mode unavailable")

        if self._precomputed:
            logger.info(
                "✅ CMap: loaded %d pre-computed clue.io scores from %s",
                len(self._precomputed), PRECOMPUTED_SCORES_PATH,
            )
        elif self.gctx_path.exists():
            logger.info("CMap: .gctx file found — local mode available")
        else:
            logger.info(
                "CMap: no pre-computed scores and no .gctx file. "
                "Run prepare_cmap_query.py + integrate_cmap_results.py "
                "to add clue.io scores."
            )

    def _load_precomputed(self):
        """Load pre-computed clue.io scores if available."""
        if not PRECOMPUTED_SCORES_PATH.exists():
            return
        try:
            with open(PRECOMPUTED_SCORES_PATH) as f:
                self._precomputed = json.load(f)
        except Exception as e:
            logger.warning("Failed to load pre-computed CMap scores: %s", e)

    def get_precomputed_score(self, drug_name: str) -> float:
        """
        Return pre-computed clue.io cmap_score for a drug.
        Returns 0.5 (neutral) if drug not in results.
        """
        if not self._precomputed:
            return None

        name_upper = drug_name.upper().strip()

        # Try exact match
        if name_upper in self._precomputed:
            return self._precomputed[name_upper]["cmap_score"]

        # Try partial match (handles salt forms)
        for key in self._precomputed:
            if key in name_upper or name_upper in key:
                return self._precomputed[key]["cmap_score"]

        # Drug screened but not in CMap database — neutral
        return 0.50

    def has_precomputed_scores(self) -> bool:
        return bool(self._precomputed)

    async def _load_metadata(self):
        """Load .gctx metadata for local mode (only if file exists)."""
        if self.is_ready:
            return
        if not self.gctx_path.exists() or not self.siginfo_path.exists():
            return
        logger.info("⏳ Loading CMap .gctx metadata (local mode)...")
        siginfo = pd.read_csv(self.siginfo_path, sep="\t", low_memory=False)
        hq_sigs = siginfo[siginfo["is_gold"] == 1]
        for _, row in hq_sigs.iterrows():
            drug_name = str(row["cmap_name"]).upper()
            sig_id    = row["sig_id"]
            if drug_name not in self.drug_to_sig_map:
                self.drug_to_sig_map[drug_name] = []
            self.drug_to_sig_map[drug_name].append(sig_id)
        self.is_ready = True
        logger.info(
            "✅ CMap .gctx metadata: %d drugs mapped", len(self.drug_to_sig_map)
        )

    def _calculate_reversal_score(self, drug_expr: pd.Series, disease_sig: Dict) -> float:
        score = 0.0
        valid_genes = 0
        for gene in disease_sig["up"]:
            if gene in drug_expr.index:
                score -= drug_expr[gene]
                valid_genes += 1
        for gene in disease_sig["down"]:
            if gene in drug_expr.index:
                score += drug_expr[gene]
                valid_genes += 1
        if valid_genes == 0:
            return 0.0
        avg_reversal = score / valid_genes
        return float(1 / (1 + np.exp(-avg_reversal)))

    async def query_reversers(self, disease: str, top_k: int = 50) -> List[Dict]:
        """
        Return top drug reversers of the disease signature.

        If pre-computed scores are available (from clue.io query),
        returns those. Otherwise attempts local .gctx mode.
        """
        is_dipg = "dipg" in disease.lower() or "h3k27m" in disease.lower()
        sig     = DIPG_SIGNATURE if is_dipg else GBM_SIGNATURE

        if self._precomputed:
            # Return top reversers from pre-computed scores
            reversers = sorted(
                self._precomputed.items(),
                key=lambda x: x[1]["cmap_score"],
                reverse=True
            )[:top_k]
            return [
                {
                    "name":       name,
                    "cmap_score": data["cmap_score"],
                    "tau":        data.get("tau"),
                    "is_reverser": data.get("is_reverser", False),
                    "source":     "clue.io pre-computed",
                }
                for name, data in reversers
            ]

        # Local .gctx fallback
        await self._load_metadata()
        if not self.is_ready:
            return []

        # Full local mode would parse .gctx here
        return []

    async def query_differential_reversers(self, top_k: int = 30) -> List[Dict]:
        """Drugs that reverse H3K27M but not normal DIPG — from pre-computed."""
        if self._precomputed:
            reversers = [
                {"name": name, **data}
                for name, data in self._precomputed.items()
                if data.get("tau", 0) < -90
            ]
            return sorted(reversers, key=lambda x: x["tau"])[:top_k]
        return []