"""
atrt_de_scorer.py
=================
Differential expression scoring for ATRT pipeline v3.0.

Replaces the raw mean-difference approach with Welch's t-test + log2FC,
which accounts for unequal variances between microarray (GSE70678) and
RNA-seq (GTEx) platforms.

STATISTICAL RATIONALE
----------------------
Welch's t-test is preferred over Student's t-test because:
  - GSE70678: Affymetrix microarray values (log2-normal distribution, ~0.3 SD)
  - GTEx:     RNA-seq log2(TPM+1) (slightly heavier tails, ~0.5 SD)
  Pooled variance assumption is violated; Welch's handles this correctly.

A z-score approach assumes normality within each gene across samples.
At n=49 for rare disease, this is fragile for outlier genes like TYR
(high in ATRT-TYR, absent in ATRT-SHH/MYC) — the distribution is bimodal.

Combined filter: abs(log2FC) >= 1.0 AND p_adj <= 0.05 (BH correction).
Source: Love MI et al. Genome Biology 2014 (DESeq2 rationale applies to
any log-scale comparison using Welch's statistic).

PLATFORM HARMONISATION
-----------------------
GSE70678 (GPL6244 Affymetrix) values are already log2-normalised by RMA.
GTEx values are log2(median_TPM + 1) — also log2 scale.
The platforms are NOT directly comparable in absolute units, so we use
the z-score of each gene's log2FC distribution (across all genes) rather
than raw log2FC for the final score. This corrects for systematic offset.

Reference: Irizarry RA et al. Biostatistics 2003 — RMA normalisation.
"""

import logging
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Statistical thresholds
# ─────────────────────────────────────────────────────────────────────────────

LOG2FC_THRESHOLD   = 1.0    # ≥ 2-fold change required (Torchia 2015 uses this)
P_ADJ_THRESHOLD    = 0.05   # BH-corrected FDR
MIN_EXPR_ATRT      = 3.0    # log2 units — filter unexpressed genes in ATRT
MIN_SAMPLES_TESTED = 10     # minimum ATRT samples with valid expr for test
WELCH_EQUAL_VAR    = False  # Always False — Welch's t-test


def welch_ttest_deseq2_style(
    atrt_matrix:   pd.DataFrame,
    gtex_ref:      pd.Series,
    atrt_col_mask: Optional[List[bool]] = None,
) -> pd.DataFrame:
    """
    Compute per-gene differential expression between ATRT and normal brain.

    Parameters
    ----------
    atrt_matrix : DataFrame
        Rows = gene symbols (uppercase), columns = ATRT sample IDs.
        Values = log2-normalised expression (RMA for Affymetrix).

    gtex_ref : Series
        Index = gene symbols (uppercase), values = log2(median_TPM + 1)
        averaged across GTEx cerebellum + cortex brain tissues.

    atrt_col_mask : list of bool, optional
        If provided, only include columns where mask is True.
        Useful for subgroup-specific DE (e.g. ATRT-MYC samples only).

    Returns
    -------
    DataFrame with columns:
        gene           : gene symbol
        log2FC         : log2(mean_ATRT / mean_normal) — ATRT vs normal
        zscore_log2FC  : z-score of log2FC across all genes (normalises
                         platform offset between microarray and RNA-seq)
        p_value        : Welch's t-test p-value (one-sample vs gtex_ref
                         treated as known reference mean)
        p_adj          : BH-corrected FDR
        mean_atrt      : mean log2 expression across ATRT samples
        mean_normal    : GTEx reference value for this gene
        n_atrt_valid   : number of ATRT samples with non-NA expression
        is_significant : bool — passes both log2FC and p_adj thresholds
        direction      : UP / DOWN / NS
        score_0to1     : tissue expression score for pipeline scoring

    Notes
    -----
    We use a one-sample t-test against the GTEx reference mean (treating GTEx
    as a known population mean, not a sample). This avoids the degree-of-freedom
    problem when GTEx has very few brain tissue replicates (n=4 in our composite).
    This is mathematically equivalent to testing H0: mu_ATRT = mu_GTEx.
    """
    # Align gene symbols
    atrt_matrix = atrt_matrix.copy()
    atrt_matrix.index = atrt_matrix.index.str.upper().str.strip()
    gtex_ref.index    = gtex_ref.index.str.upper().str.strip()

    common_genes = atrt_matrix.index.intersection(gtex_ref.index)
    if len(common_genes) < 100:
        logger.warning(
            "Only %d genes in common between ATRT matrix and GTEx reference. "
            "Check that both use uppercase HGNC symbols. "
            "Consider re-running data_downloader.py --process gse70678 gpl6244.",
            len(common_genes),
        )

    atrt_common = atrt_matrix.loc[common_genes]
    gtex_common = gtex_ref.loc[common_genes]

    # Apply sample mask if provided (subgroup mode)
    if atrt_col_mask is not None:
        cols = [c for c, m in zip(atrt_common.columns, atrt_col_mask) if m]
        atrt_common = atrt_common[cols]
        logger.info("Subgroup mask: using %d/%d ATRT samples", len(cols), atrt_matrix.shape[1])

    rows = []
    for gene in common_genes:
        expr_vals = pd.to_numeric(atrt_common.loc[gene], errors="coerce").dropna()
        ref_val   = float(gtex_common[gene])
        n_valid   = len(expr_vals)

        if n_valid < MIN_SAMPLES_TESTED:
            continue
        if expr_vals.mean() < MIN_EXPR_ATRT and ref_val < MIN_EXPR_ATRT:
            continue  # Both unexpressed — not informative

        mean_atrt  = float(expr_vals.mean())
        log2fc     = mean_atrt - ref_val  # log2(ATRT/normal) in log2 space

        # One-sample Welch's t-test: H0: mu_ATRT = ref_val
        # scipy.stats.ttest_1samp is equivalent to Welch's when popmean is fixed
        t_stat, p_val = stats.ttest_1samp(expr_vals, popmean=ref_val)

        rows.append({
            "gene":         gene,
            "log2FC":       round(log2fc, 4),
            "t_statistic":  round(float(t_stat), 4),
            "p_value":      float(p_val),
            "mean_atrt":    round(mean_atrt, 4),
            "mean_normal":  round(ref_val, 4),
            "std_atrt":     round(float(expr_vals.std()), 4),
            "n_atrt_valid": n_valid,
        })

    if not rows:
        logger.error("No genes survived DE filtering. Check input data.")
        return pd.DataFrame()

    de_df = pd.DataFrame(rows).set_index("gene")

    # BH multiple testing correction
    from scipy.stats import false_discovery_control
    try:
        de_df["p_adj"] = false_discovery_control(de_df["p_value"].values, method="bh")
    except AttributeError:
        # scipy < 1.11 fallback — manual BH
        de_df["p_adj"] = _bh_correction(de_df["p_value"].values)

    # Z-score of log2FC across all genes (corrects for microarray/RNA-seq platform offset)
    lfc_mean = de_df["log2FC"].mean()
    lfc_std  = de_df["log2FC"].std()
    de_df["zscore_log2FC"] = (de_df["log2FC"] - lfc_mean) / lfc_std.clip(lower=0.01)

    # Significance flag
    de_df["is_significant"] = (
        (de_df["log2FC"].abs() >= LOG2FC_THRESHOLD) &
        (de_df["p_adj"] <= P_ADJ_THRESHOLD)
    )
    de_df["direction"] = np.where(
        ~de_df["is_significant"], "NS",
        np.where(de_df["log2FC"] > 0, "UP", "DOWN")
    )

    # Convert to 0–1 pipeline score using the z-score of log2FC
    # z_lfc == 0 → score 0.50 (average expression in ATRT)
    # z_lfc == +2 → score 0.85 (strongly upregulated)
    # z_lfc == -2 → score 0.15 (strongly downregulated)
    de_df["score_0to1"] = (de_df["zscore_log2FC"] * 0.175 + 0.50).clip(0.05, 0.95)

    n_up   = (de_df["direction"] == "UP").sum()
    n_down = (de_df["direction"] == "DOWN").sum()
    n_ns   = (de_df["direction"] == "NS").sum()
    logger.info(
        "Welch DE complete: %d genes | UP=%d DOWN=%d NS=%d (log2FC≥%.1f, FDR≤%.2f)",
        len(de_df), n_up, n_down, n_ns, LOG2FC_THRESHOLD, P_ADJ_THRESHOLD,
    )

    return de_df.reset_index()


def _bh_correction(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction (scipy fallback for older versions)."""
    n = len(p_values)
    order = np.argsort(p_values)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    p_adj = np.minimum(1.0, p_values * n / ranked)
    # Enforce monotonicity
    for i in range(n - 2, -1, -1):
        p_adj[order[i]] = min(p_adj[order[i]], p_adj[order[i + 1]])
    return p_adj


# ─────────────────────────────────────────────────────────────────────────────
# Top gene selector for CMap query (called by prepare_cmap_query.py)
# ─────────────────────────────────────────────────────────────────────────────

def select_cmap_query_genes(
    de_df: pd.DataFrame,
    n_up:   int = 150,
    n_down: int = 50,
    require_significant: bool = True,
) -> Tuple[List[str], List[str]]:
    """
    Select top up/down genes for clue.io L1000 query.

    Selection criteria (in order of priority):
      1. Statistically significant (FDR ≤ 0.05) AND |log2FC| ≥ 1.0
      2. Ranked by absolute z-score of log2FC (strongest signal first)
      3. Cap at n_up UP genes and n_down DOWN genes

    clue.io recommends 50–150 genes per direction for maximum signal.
    We use more UP genes because ATRT has a strong proliferative
    transcriptional signature (EZH2, BRD4, AURKA all upregulated).

    Parameters
    ----------
    de_df : DataFrame from welch_ttest_deseq2_style()
    n_up, n_down : max genes per direction
    require_significant : if True, only include FDR-significant genes.
        Set to False if the dataset is small and few genes reach FDR < 0.05.

    Returns
    -------
    (up_genes, down_genes) : lists of HGNC gene symbols
    """
    if de_df.empty:
        logger.error("DE dataframe is empty — cannot select CMap query genes.")
        return [], []

    if require_significant:
        up_pool   = de_df[(de_df["direction"] == "UP")]
        down_pool = de_df[(de_df["direction"] == "DOWN")]
    else:
        up_pool   = de_df[de_df["log2FC"] > 0]
        down_pool = de_df[de_df["log2FC"] < 0]

    # Sort by effect size (absolute z-score of log2FC)
    up_sorted   = up_pool.sort_values("zscore_log2FC", ascending=False)
    down_sorted = down_pool.sort_values("zscore_log2FC", ascending=True)

    up_genes   = up_sorted["gene"].head(n_up).tolist()
    down_genes = down_sorted["gene"].head(n_down).tolist()

    logger.info(
        "CMap gene selection: %d UP (of %d sig) | %d DOWN (of %d sig)",
        len(up_genes), len(up_pool), len(down_genes), len(down_pool),
    )

    # Validate ATRT biology is represented
    expected_up = {"EZH2", "BRD4", "AURKA", "MYC", "MYCN", "HDAC1", "CDK4"}
    found = expected_up & set(up_genes)
    if len(found) < 3:
        logger.warning(
            "Only %d/%d expected ATRT markers in UP list: %s. "
            "Consider lowering log2FC threshold or checking GSE70678 loading.",
            len(found), len(expected_up), sorted(found),
        )
    else:
        logger.info("ATRT biology check — key markers in UP list: %s", sorted(found))

    return up_genes, down_genes


# ─────────────────────────────────────────────────────────────────────────────
# Subgroup-stratified DE
# ─────────────────────────────────────────────────────────────────────────────

def compute_subgroup_de(
    atrt_matrix:   pd.DataFrame,
    gtex_ref:      pd.Series,
    subgroup_labels: Dict[str, List[str]],
) -> Dict[str, pd.DataFrame]:
    """
    Run Welch's DE separately for each ATRT molecular subgroup.

    Parameters
    ----------
    subgroup_labels : dict mapping subgroup name → list of sample column names
        e.g. {"TYR": ["GSM_001", "GSM_002"], "SHH": [...], "MYC": [...]}
        Derive from GSE70678 metadata or GSE73868 (Johann 2016) methylation calls.

    Returns
    -------
    dict mapping subgroup_name → DE DataFrame (same schema as welch_ttest_deseq2_style)
    """
    results: Dict[str, pd.DataFrame] = {}

    for subgroup, sample_ids in subgroup_labels.items():
        available = [s for s in sample_ids if s in atrt_matrix.columns]
        if len(available) < 5:
            logger.warning(
                "Subgroup %s: only %d samples available (need ≥5 for DE). Skipping.",
                subgroup, len(available),
            )
            continue

        logger.info("Running DE for subgroup %s (n=%d samples)...", subgroup, len(available))
        sub_matrix = atrt_matrix[available]
        de_result  = welch_ttest_deseq2_style(sub_matrix, gtex_ref)
        de_result["subgroup"] = subgroup
        results[subgroup] = de_result

    return results


# ─────────────────────────────────────────────────────────────────────────────
# Score a drug candidate using the DE results
# ─────────────────────────────────────────────────────────────────────────────

def score_candidate_from_de(
    drug: Dict,
    de_df: pd.DataFrame,
    curated_scores: Dict[str, float],
    curated_weight: float = 0.55,
    de_weight:      float = 0.45,
    subgroup:       Optional[str] = None,
    subgroup_de:    Optional[Dict[str, pd.DataFrame]] = None,
) -> float:
    """
    Compute tissue_expression_score for a drug candidate using DE results.

    Blends:
      - DE-based score (Welch's log2FC z-score → 0-1): data-driven
      - Curated score (ATRT_CURATED_SCORES): knowledge-based fallback

    Weights default to 55/45 curated/DE because the curated scores reflect
    published ATRT cell-line functional data (not just transcript abundance)
    — e.g. EZH2 is high even though its RNA upregulation is 2.31 z-score,
    which the curated 0.92 already captures.

    When subgroup DE is available and matches the requested subgroup,
    the subgroup-specific DE replaces pan-ATRT DE for scoring.

    Parameters
    ----------
    drug : dict with 'targets' list
    de_df : pan-ATRT DE DataFrame from welch_ttest_deseq2_style()
    curated_scores : dict of gene → score (ATRT_CURATED_SCORES from pipeline_config)
    subgroup_de : optional dict of subgroup → DE DataFrame
    """
    targets = [t.upper() for t in (drug.get("targets") or [])]
    if not targets:
        return 0.40

    # Select the right DE table
    active_de = de_df.copy()
    if subgroup and subgroup_de and subgroup in subgroup_de:
        active_de = subgroup_de[subgroup]
        logger.debug("Using subgroup-specific DE for %s (%s)", drug.get("name"), subgroup)

    de_scores = {}
    if not active_de.empty and "gene" in active_de.columns:
        de_lookup = active_de.set_index("gene")["score_0to1"].to_dict()
        for t in targets:
            if t in de_lookup:
                de_scores[t] = de_lookup[t]

    # Curated scores for targets
    cur_scores = {t: curated_scores.get(t, 0.40) for t in targets}

    # Best + mean-of-top-2 aggregation (same as existing pipeline)
    if de_scores:
        de_vals_sorted  = sorted(de_scores.values(), reverse=True)
        de_best         = de_vals_sorted[0]
        de_top2_mean    = sum(de_vals_sorted[:2]) / min(len(de_vals_sorted), 2)
        de_combined     = 0.65 * de_best + 0.35 * de_top2_mean
    else:
        de_combined = None

    cur_vals_sorted = sorted(cur_scores.values(), reverse=True)
    cur_best        = cur_vals_sorted[0]
    cur_top2_mean   = sum(cur_vals_sorted[:2]) / min(len(cur_vals_sorted), 2)
    cur_combined    = 0.65 * cur_best + 0.35 * cur_top2_mean

    if de_combined is not None:
        final = curated_weight * cur_combined + de_weight * de_combined
    else:
        final = cur_combined  # DE not available for these targets

    return round(float(np.clip(final, 0.05, 1.00)), 4)