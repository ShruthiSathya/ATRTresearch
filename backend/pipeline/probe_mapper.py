"""
probe_mapper.py
===============
Maps Affymetrix HuGene 1.0 ST (GPL6244) probe IDs → HGNC gene symbols.

WHY THIS IS NEEDED
-------------------
GSE70678 (Torchia 2015 Cancer Cell) was generated on the Affymetrix HuGene 1.0
ST Array (GPL6244). The GEO series_matrix.txt file uses Affymetrix probe set IDs
(e.g., "200099_s_at") as the row index — NOT gene symbols (e.g., "EZH2").

This module downloads the GPL6244 platform annotation file from NCBI GEO,
parses the probe → gene symbol mapping, and transforms a probe-indexed
expression DataFrame into a gene-symbol-indexed one.

ANNOTATION FILE
----------------
Source : NCBI GEO FTP
URL    : https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz
Size   : ~5 MB compressed
License: Public domain (NCBI GEO)

The annotation file contains columns including:
  ID           — Affymetrix probe set ID   (e.g., "200099_s_at")
  Gene symbol  — HGNC gene symbol(s)       (e.g., "UBC" or "BRCA1 /// BRCA2")
  Gene ID      — NCBI Entrez Gene ID

Multi-gene probes (e.g., "BRCA1 /// BRCA2") are handled by taking the
first annotated symbol. Control probes (AFFX-*, no gene annotation) are
filtered out.

AGGREGATION
------------
Multiple probe sets often target the same gene. After mapping, probes are
aggregated per gene using the median across all probes for that gene.
Median is preferred over mean because it is more robust to outlier probes
caused by cross-hybridisation or probe design artefacts.

References
-----------
Carvalho BS & Irizarry RA (2010). A framework for oligonucleotide microarray
  preprocessing. Bioinformatics 26(19):2363-2367. PMID 20688976.
Torchia J et al. (2015). Cancer Cell 30(6):891-908. PMID 26609405.
"""

import gzip
import io
import logging
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

GPL6244_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
    "GPL6nnn/GPL6244/annot/GPL6244.annot.gz"
)
DEFAULT_CACHE = Path("data/raw_omics/GPL6244_probe_map.tsv")

# Column names that appear in the GPL6244 annotation file
_ID_CANDIDATES     = ["ID", "id", "probe_id", "ProbeID", "Probe Set ID"]
_SYMBOL_CANDIDATES = ["Gene symbol", "gene_symbol", "Symbol", "GENE_SYMBOL",
                       "GeneSymbol", "gene symbol"]
_ENTREZ_CANDIDATES = ["Gene ID", "gene_id", "ENTREZ_GENE_ID", "Entrez Gene",
                       "entrez_gene_id"]

# Affymetrix control probe prefixes — always excluded
_CONTROL_PREFIXES  = ("AFFX-", "affx-")


# ─────────────────────────────────────────────────────────────────────────────
# Helper
# ─────────────────────────────────────────────────────────────────────────────

def _find_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """Return the first column name matching any candidate (case-insensitive)."""
    col_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        if cand.lower() in col_lower:
            return col_lower[cand.lower()]
    return None


def _first_symbol(raw: str) -> str:
    """
    Extract the first gene symbol from a possibly multi-valued annotation.

    GPL6244 uses ' /// ' as a delimiter for multi-gene probes, e.g.:
      "BRCA1 /// BRCA2"  →  "BRCA1"
      "TP53"             →  "TP53"
      "---"              →  ""   (unmapped probe)
    """
    if not raw or str(raw).strip() in ("---", "nan", "NaN", "", "None"):
        return ""
    sym = str(raw).split(" /// ")[0].strip()
    # Some annotations use '/' as a sub-separator (isoforms) — take first
    sym = sym.split("/")[0].strip()
    return sym.upper()


# ─────────────────────────────────────────────────────────────────────────────
# ProbeMapper
# ─────────────────────────────────────────────────────────────────────────────

class ProbeMapper:
    """
    Affymetrix HuGene 1.0 ST (GPL6244) probe-to-gene mapper.

    Usage
    -----
    >>> mapper = ProbeMapper()
    >>> mapper.load()                        # downloads annotation if needed
    >>> gene_df = mapper.map_probes_to_genes(probe_df)  # probe → gene matrix

    The mapper is stateless after `load()` — call once per process.
    """

    def __init__(self, cache_path: Path = DEFAULT_CACHE):
        self.cache_path            = Path(cache_path)
        self._probe_to_symbol: Dict[str, str] = {}
        self._probe_to_entrez: Dict[str, str] = {}
        self._is_loaded            = False

    # ── Public API ────────────────────────────────────────────────────────────

    def load(self, force_download: bool = False) -> None:
        """
        Load the GPL6244 probe→symbol mapping.

        First checks for a local cache (data/raw_omics/GPL6244_probe_map.tsv).
        Downloads from NCBI GEO FTP if not cached or force_download=True.

        Parameters
        ----------
        force_download : bool
            Re-download even if cache exists.
        """
        if self._is_loaded and not force_download:
            return

        if self.cache_path.exists() and not force_download:
            logger.info("Loading GPL6244 probe map from cache: %s", self.cache_path)
            self._load_from_cache()
        else:
            logger.info(
                "Downloading GPL6244 annotation from NCBI GEO (%s)...", GPL6244_URL
            )
            self._download_and_parse()

        n_total  = len(self._probe_to_symbol)
        n_mapped = sum(1 for s in self._probe_to_symbol.values() if s)
        logger.info(
            "GPL6244 annotation loaded: %d total probes, %d with gene symbols (%.1f%%)",
            n_total, n_mapped, 100 * n_mapped / max(n_total, 1),
        )
        self._is_loaded = True

    def get_symbol(self, probe_id: str) -> str:
        """Return gene symbol for a single probe ID; empty string if unmapped."""
        if not self._is_loaded:
            self.load()
        return self._probe_to_symbol.get(str(probe_id).strip(), "")

    def map_probes_to_genes(
        self,
        probe_df: pd.DataFrame,
        aggregation: str = "median",
    ) -> pd.DataFrame:
        """
        Convert a probe-indexed expression DataFrame to gene-symbol-indexed.

        Parameters
        ----------
        probe_df : pd.DataFrame
            Rows  = probe IDs (e.g., "200099_s_at").
            Cols  = sample identifiers (e.g., "GSM1722234").
            Values assumed to be log2-transformed expression.
        aggregation : str
            How to combine multiple probes mapping to the same gene.
            "median" (default, recommended) | "mean" | "max"

        Returns
        -------
        pd.DataFrame
            Rows  = HGNC gene symbols (uppercase).
            Cols  = same sample identifiers as input.
            Index.name = "gene_symbol"

        Notes
        -----
        * Control probes (AFFX-*) are excluded.
        * Probes with no gene annotation are excluded.
        * Genes covered by only 1 probe still appear in the output.
        """
        if not self._is_loaded:
            self.load()

        probe_df = probe_df.copy()
        probe_df.index = probe_df.index.astype(str).str.strip()

        # Annotate each probe with its gene symbol
        probe_df.insert(
            0, "_sym",
            probe_df.index.map(lambda p: self._probe_to_symbol.get(p, ""))
        )

        # Filter unmapped and control probes
        n_before = len(probe_df)
        probe_df = probe_df[
            (probe_df["_sym"] != "")
            & ~probe_df.index.str.startswith(_CONTROL_PREFIXES)
        ]
        n_after  = len(probe_df)
        logger.info(
            "Probe mapping: %d probes → kept %d (removed %d unmapped/control)",
            n_before, n_after, n_before - n_after,
        )

        # Ensure numeric values only in sample columns
        sample_cols = [c for c in probe_df.columns if c != "_sym"]
        for col in sample_cols:
            probe_df[col] = pd.to_numeric(probe_df[col], errors="coerce")

        # Aggregate multiple probes per gene
        agg_fn = {"median": "median", "mean": "mean", "max": "max"}.get(aggregation)
        if agg_fn is None:
            raise ValueError(f"aggregation must be 'median', 'mean', or 'max'; got '{aggregation}'")

        gene_df = probe_df.groupby("_sym")[sample_cols].agg(agg_fn)
        gene_df.index.name = "gene_symbol"

        logger.info(
            "After probe aggregation (%s): %d unique genes × %d samples",
            aggregation, len(gene_df), len(sample_cols),
        )
        return gene_df

    def get_coverage(self, probe_ids: List[str]) -> Dict:
        """Summarise what fraction of the supplied probe IDs were successfully mapped."""
        if not self._is_loaded:
            self.load()
        total  = len(probe_ids)
        mapped = sum(1 for p in probe_ids if self._probe_to_symbol.get(str(p)))
        return {
            "total_probes":  total,
            "mapped_probes": mapped,
            "unmapped":      total - mapped,
            "coverage_pct":  round(100 * mapped / max(total, 1), 1),
        }

    # ── Private helpers ───────────────────────────────────────────────────────

    def _download_and_parse(self) -> None:
        """Download GPL6244.annot.gz, parse it, build mappings, cache result."""
        try:
            with urllib.request.urlopen(GPL6244_URL, timeout=120) as resp:
                raw = resp.read()
            logger.info("Downloaded GPL6244 annotation (%.1f MB)", len(raw) / 1e6)
        except Exception as e:
            raise RuntimeError(
                f"Failed to download GPL6244 annotation from {GPL6244_URL}: {e}\n"
                "Check your internet connection or download manually and place at:\n"
                f"  {self.cache_path}"
            ) from e

        # Decompress
        try:
            with gzip.open(io.BytesIO(raw), "rt", encoding="utf-8", errors="replace") as fh:
                lines = fh.readlines()
        except Exception as e:
            raise RuntimeError(f"Failed to decompress GPL6244 annotation: {e}") from e

        # Separate comment header from data
        # The GPL6244 annotation format:
        #   Lines beginning with # are metadata/comments
        #   First non-comment line is the column header
        #   Subsequent lines are data rows
        data_lines = [ln for ln in lines if not ln.startswith("#") and ln.strip()]
        if not data_lines:
            raise ValueError("GPL6244 annotation: no data lines found after comment header")

        logger.info("Parsing GPL6244 annotation (%d data lines)...", len(data_lines))

        try:
            df = pd.read_csv(
                io.StringIO("".join(data_lines)),
                sep="\t",
                dtype=str,
                low_memory=False,
                na_filter=False,
            )
        except Exception as e:
            raise RuntimeError(f"Failed to parse GPL6244 annotation as TSV: {e}") from e

        df.columns = df.columns.str.strip()
        logger.debug("GPL6244 columns: %s", list(df.columns[:10]))

        id_col     = _find_col(df, _ID_CANDIDATES)
        symbol_col = _find_col(df, _SYMBOL_CANDIDATES)
        entrez_col = _find_col(df, _ENTREZ_CANDIDATES)

        if id_col is None:
            raise ValueError(
                f"Could not find probe ID column. Available columns: {list(df.columns)}"
            )
        if symbol_col is None:
            raise ValueError(
                f"Could not find gene symbol column. Available columns: {list(df.columns)}"
            )

        logger.info("Probe ID column: '%s' | Symbol column: '%s'", id_col, symbol_col)

        # Build mappings
        for _, row in df.iterrows():
            probe = str(row[id_col]).strip()
            sym   = _first_symbol(str(row[symbol_col]))
            self._probe_to_symbol[probe] = sym
            if entrez_col:
                self._probe_to_entrez[probe] = str(row[entrez_col]).strip()

        # Cache to disk
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_df = pd.DataFrame([
            {
                "probe_id":    k,
                "gene_symbol": v,
                "entrez_id":   self._probe_to_entrez.get(k, ""),
            }
            for k, v in self._probe_to_symbol.items()
        ])
        cache_df.to_csv(self.cache_path, sep="\t", index=False)
        logger.info("GPL6244 probe map cached → %s (%d rows)", self.cache_path, len(cache_df))

    def _load_from_cache(self) -> None:
        """Load previously parsed probe→symbol mapping from TSV cache."""
        df = pd.read_csv(self.cache_path, sep="\t", dtype=str, na_filter=False)
        self._probe_to_symbol = dict(zip(df["probe_id"], df["gene_symbol"]))
        if "entrez_id" in df.columns:
            self._probe_to_entrez = dict(zip(df["probe_id"], df["entrez_id"]))


# ─────────────────────────────────────────────────────────────────────────────
# Series matrix parser (needed before ProbeMapper.map_probes_to_genes)
# ─────────────────────────────────────────────────────────────────────────────

def parse_geo_series_matrix(path: str | Path) -> tuple[pd.DataFrame, dict]:
    """
    Parse a GEO series_matrix.txt file into expression data + sample metadata.

    GEO series matrix format:
        Lines starting with '!' are key-value metadata.
        Data is between '!series_matrix_table_begin' and '!series_matrix_table_end'.
        First data row is the column header ('ID_REF' + GSM accessions).
        Subsequent rows are probe ID + float expression values.

    Parameters
    ----------
    path : str or Path
        Path to the series_matrix.txt file (plain text or .gz).

    Returns
    -------
    (expr_df, meta)
        expr_df : DataFrame with probe IDs as index, GSM IDs as columns.
        meta    : dict of sample metadata keyed by field name.
    """
    path = Path(path)
    logger.info("Parsing GEO series matrix: %s", path)

    # Handle gzipped files
    if str(path).endswith(".gz"):
        import gzip as _gz
        open_fn = lambda p: _gz.open(p, "rt", encoding="utf-8", errors="replace")
    else:
        open_fn = lambda p: open(p, "r", encoding="utf-8", errors="replace")

    metadata: dict = {}
    in_table = False
    table_lines: list[str] = []

    with open_fn(path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            # Metadata lines
            if line.startswith("!"):
                key, _, val = line.partition("\t")
                key = key.lstrip("!")
                if key not in metadata:
                    metadata[key] = []
                # GEO metadata values are quoted — strip them
                metadata[key].append(val.strip('"'))
                continue

            # Data table markers
            if line.startswith("!series_matrix_table_begin") or line.startswith("\"!series_matrix_table_begin\""):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end") or line.startswith("\"!series_matrix_table_end\""):
                break

            if in_table and line.strip():
                table_lines.append(line)

    if not table_lines:
        raise ValueError(
            f"No data table found in {path}. "
            "Expected '!series_matrix_table_begin' / '!series_matrix_table_end' markers."
        )

    # Parse the expression table
    expr_df = pd.read_csv(
        io.StringIO("\n".join(table_lines)),
        sep="\t",
        index_col=0,
        dtype=str,
        na_filter=False,
    )
    expr_df.index.name = "probe_id"
    # Strip quotes from index and column names (GEO sometimes quotes them)
    expr_df.index   = expr_df.index.str.strip('"')
    expr_df.columns = expr_df.columns.str.strip('"')

    # Convert to numeric
    for col in expr_df.columns:
        expr_df[col] = pd.to_numeric(expr_df[col], errors="coerce")

    logger.info(
        "Series matrix parsed: %d probes × %d samples",
        len(expr_df), len(expr_df.columns),
    )
    return expr_df, metadata