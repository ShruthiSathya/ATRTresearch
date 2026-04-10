"""
normal_brain_fetcher.py
========================
Loads the GTEx v8 normal brain reference for differential expression.

WHY GTEx FOR NORMAL BRAIN?
----------------------------
GSE70678 contains only 49 ATRT tumour samples — no matched normal brain controls
are included. Computing differential expression requires a normal brain baseline.

GTEx v8 is the current gold-standard reference for human tissue expression:
  n = 209 cerebellum donors, n = 255 cortex donors (post-mortem, non-diseased)
  Gene-level RNA-seq, median TPM per tissue (no individual sample noise)
  Open access: https://gtexportal.org

After log2(TPM + 1) transformation, GTEx values are roughly comparable to the
Affymetrix log2 scale from GSE70678, with the important caveat that they are
from RNA-seq vs microarray — this is corrected for by using z-score
differential rather than raw fold-change where possible.

BRAIN TISSUES USED
-------------------
Selected on anatomical relevance to ATRT tumour locations:
  Brain - Cerebellum             (infratentorial ATRT, ~50%)
  Brain - Cerebellar Hemisphere  (infratentorial ATRT)
  Brain - Cortex                 (supratentorial ATRT, ~35%)
  Brain - Frontal Cortex (BA9)   (supratentorial ATRT)

Mean across selected tissues = composite "normal brain" reference.

References
-----------
GTEx Consortium (2020). The GTEx Consortium atlas of genetic regulatory effects
  across human tissues. Science 369(6509):1318-1330. PMID 32913098.
Frühwald MC et al. (2020). CNS Oncology 9(2):CNS56. PMID 32432484. [ATRT location]
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Paths and constants
# ─────────────────────────────────────────────────────────────────────────────

GTEX_PROCESSED_PATH = Path("data/raw_omics/GTEx_brain_normal_reference.tsv")
GTEX_RAW_PATH       = Path("data/raw_omics/GTEx_v8_gene_median_tpm.gct.gz")

# Brain tissues to include in the normal reference composite.
# Names must match the column headers in the GTEx GCT file exactly.
BRAIN_TISSUES = [
    "Brain - Cerebellum",
    "Brain - Cerebellar Hemisphere",
    "Brain - Cortex",
    "Brain - Frontal Cortex (BA9)",
]

# How to combine multiple tissues into a single "normal brain" vector
# Options: "mean" | "median" | "cerebellum_only"
COMPOSITE_STRATEGY = "mean"


class NormalBrainFetcher:
    """
    Provides a gene-level normal brain expression reference from GTEx v8.

    Usage
    -----
    >>> fetcher = NormalBrainFetcher()
    >>> ref = fetcher.get_reference()   # dict: gene_symbol → normal_log2_expr
    >>> diff = fetcher.compute_differential(atrt_gene_df)
    """

    def __init__(
        self,
        processed_path: Path = GTEX_PROCESSED_PATH,
        tissues: List[str] = BRAIN_TISSUES,
        strategy: str = COMPOSITE_STRATEGY,
    ):
        self.processed_path = Path(processed_path)
        self.tissues        = tissues
        self.strategy       = strategy
        self._ref: Optional[pd.Series] = None     # gene → mean log2 expr
        self._tissue_df: Optional[pd.DataFrame] = None   # gene × tissue
        self._is_loaded = False

    # ── Public API ────────────────────────────────────────────────────────────

    def load(self) -> None:
        """Load the GTEx brain reference into memory."""
        if self._is_loaded:
            return

        if self.processed_path.exists():
            self._load_processed()
        elif GTEX_RAW_PATH.exists():
            logger.info(
                "Processed GTEx reference not found — processing raw file (%s)...",
                GTEX_RAW_PATH,
            )
            from data_downloader import process_gtex
            if process_gtex():
                self._load_processed()
            else:
                raise RuntimeError(
                    "Failed to process GTEx raw file. "
                    "See data_downloader.process_gtex() log for details."
                )
        else:
            raise FileNotFoundError(
                f"GTEx brain reference not found at {self.processed_path} and "
                f"raw file not found at {GTEX_RAW_PATH}.\n"
                "Download with:\n"
                "  python -m backend.pipeline.data_downloader --dataset gtex_brain\n"
                "  python -m backend.pipeline.data_downloader --process gtex"
            )

        self._is_loaded = True
        logger.info(
            "Normal brain reference loaded: %d genes | strategy=%s | tissues=%s",
            len(self._ref), self.strategy, list(self._tissue_df.columns),
        )

    def get_reference(self) -> pd.Series:
        """
        Return a Series mapping gene symbol → composite normal brain log2 expression.

        Values are log2(median_TPM + 1) averaged across selected brain tissues.
        Index = uppercase HGNC gene symbols.
        """
        if not self._is_loaded:
            self.load()
        return self._ref.copy()

    def get_expression(self, gene: str) -> Optional[float]:
        """Return the normal brain log2 expression for a single gene; None if absent."""
        if not self._is_loaded:
            self.load()
        return self._ref.get(gene.upper())

    def compute_differential(
        self,
        atrt_gene_df: pd.DataFrame,
        method: str = "mean_diff",
    ) -> pd.Series:
        """
        Compute differential expression: ATRT mean − normal brain reference.

        Parameters
        ----------
        atrt_gene_df : DataFrame
            Rows = gene symbols (uppercase), columns = sample IDs.
            Values = log2 expression (same scale as Affymetrix data).
        method : str
            "mean_diff"  — simple mean(ATRT) − mean(normal)  [default]
            "zscore"     — (mean(ATRT) − mean(normal)) / std(normal)

        Returns
        -------
        pd.Series
            Index = gene symbols, values = differential expression scores.
            Positive = upregulated in ATRT vs normal brain.
            Negative = downregulated in ATRT vs normal brain.
        """
        if not self._is_loaded:
            self.load()

        # Ensure gene index is uppercase
        atrt_mean = atrt_gene_df.apply(pd.to_numeric, errors="coerce").mean(axis=1)
        atrt_mean.index = atrt_mean.index.str.upper()

        # Align to genes present in both datasets
        common = atrt_mean.index.intersection(self._ref.index)
        if len(common) == 0:
            raise ValueError(
                "No genes in common between ATRT data and GTEx reference. "
                "Check that both use uppercase HGNC gene symbols."
            )

        atrt_aligned   = atrt_mean.loc[common]
        normal_aligned = self._ref.loc[common]

        if method == "mean_diff":
            diff = atrt_aligned - normal_aligned
        elif method == "zscore":
            # Use per-tissue std across tissue columns for the normal variation
            if self._tissue_df is not None:
                tissue_aligned = self._tissue_df.loc[
                    self._tissue_df.index.intersection(common)
                ]
                normal_std = tissue_aligned.std(axis=1).clip(lower=0.01)
                diff = (atrt_aligned - normal_aligned) / normal_std
            else:
                diff = atrt_aligned - normal_aligned
                logger.warning(
                    "method='zscore' requested but tissue-level SD unavailable — "
                    "falling back to mean_diff"
                )
        else:
            raise ValueError(f"method must be 'mean_diff' or 'zscore', got '{method}'")

        n_up   = int((diff > 0).sum())
        n_down = int((diff < 0).sum())
        logger.info(
            "Differential expression (ATRT vs normal brain, %s): "
            "%d genes total | %d upregulated | %d downregulated",
            method, len(diff), n_up, n_down,
        )
        return diff

    def get_top_upregulated(self, diff: pd.Series, n: int = 150) -> List[str]:
        """Return top-n upregulated gene symbols sorted by differential score."""
        return diff.sort_values(ascending=False).head(n).index.tolist()

    def get_top_downregulated(self, diff: pd.Series, n: int = 50) -> List[str]:
        """Return top-n downregulated gene symbols sorted by differential score."""
        return diff.sort_values(ascending=True).head(n).index.tolist()

    def coverage_report(self, gene_list: List[str]) -> Dict:
        """Report what fraction of query genes have normal brain reference values."""
        if not self._is_loaded:
            self.load()
        query   = {g.upper() for g in gene_list}
        covered = query & set(self._ref.index)
        return {
            "n_query":   len(query),
            "n_covered": len(covered),
            "pct":       round(100 * len(covered) / max(len(query), 1), 1),
            "missing":   sorted(query - covered),
        }

    # ── Private helpers ───────────────────────────────────────────────────────

    def _load_processed(self) -> None:
        """Load the compact processed GTEx brain TSV."""
        df = pd.read_csv(self.processed_path, sep="\t", index_col=0)
        df.index = df.index.astype(str).str.upper().str.strip()
        df.index.name = "gene_symbol"

        # Select columns matching our requested tissues
        avail_tissues = list(df.columns)
        selected = []
        for tissue in self.tissues:
            matches = [c for c in avail_tissues if tissue.lower() in c.lower()]
            selected.extend(matches)
        selected = list(dict.fromkeys(selected))   # deduplicate, preserve order

        if not selected:
            logger.warning(
                "None of the requested tissues found in GTEx reference. "
                "Available: %s. Using all available columns.",
                avail_tissues,
            )
            selected = avail_tissues

        logger.info("GTEx brain tissues selected: %s", selected)

        tissue_df = df[selected].apply(pd.to_numeric, errors="coerce")
        self._tissue_df = tissue_df

        # Build composite reference vector
        if self.strategy == "mean":
            self._ref = tissue_df.mean(axis=1)
        elif self.strategy == "median":
            self._ref = tissue_df.median(axis=1)
        elif self.strategy == "cerebellum_only":
            cereb_cols = [c for c in selected if "cerebellum" in c.lower()]
            if cereb_cols:
                self._ref = tissue_df[cereb_cols].mean(axis=1)
            else:
                logger.warning("No cerebellum column found; using all tissues mean")
                self._ref = tissue_df.mean(axis=1)
        else:
            raise ValueError(f"Unknown strategy: {self.strategy}")

        self._ref = self._ref.dropna()
        self._ref.name = "normal_brain_log2"


# ─────────────────────────────────────────────────────────────────────────────
# Convenience function
# ─────────────────────────────────────────────────────────────────────────────

_DEFAULT_FETCHER: Optional[NormalBrainFetcher] = None


def get_normal_brain_reference() -> pd.Series:
    """
    Module-level convenience: return the default normal brain reference Series.
    Lazy-loads and caches a single NormalBrainFetcher instance.
    """
    global _DEFAULT_FETCHER
    if _DEFAULT_FETCHER is None:
        _DEFAULT_FETCHER = NormalBrainFetcher()
    if not _DEFAULT_FETCHER._is_loaded:
        _DEFAULT_FETCHER.load()
    return _DEFAULT_FETCHER.get_reference()