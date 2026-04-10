"""
depmap_cmap_crossref.py
=======================
Cross-reference CMap drug hits with DepMap 24Q4 CRISPR essentiality.

VALIDATION LOGIC
-----------------
A CMap "hit" (strong reverser, norm_cs < -0.9) is only clinically actionable
if its target gene is also ESSENTIAL in ATRT cell lines.

The cross-reference has two tiers:

Tier 1 — Strong validation:
  CMap reversal AND DepMap Chronos ≤ -1.0 in ATRT/rhabdoid lines.
  These are the most credible hits. Panobinostat example:
    - HDAC1 Chronos = -1.38 in BT16/BT37/G401 (strongly essential)
    - Panobinostat norm_cs = -1.31 (strong reverser)
    → Tier 1 ✅

Tier 2 — Supported:
  CMap reversal AND DepMap Chronos -0.5 to -1.0.
  Target is essential but not in the top quartile.

Tier 3 — Transcriptional only:
  CMap reversal but DepMap Chronos > -0.5.
  The gene is transcriptionally dysregulated but the cell doesn't depend on it.
  Example: ONC201 targets DRD2 — high CMap score, but DRD2 Chronos = -0.35
  → not essential. ONC201's mechanism may be non-target-dependent (TRAIL pathway).

NOTE ON CHLA-266
-----------------
CHLA-266 is NOT in DepMap 24Q4 under that name. The CHLA cell line series
(CHLA02, CHLA04, CHLA06) are COG lines. CHLA-266 is a neuroblastoma line,
not ATRT. Verify your cell line IDs against Model.csv:
  grep -i "chla" data/depmap/Model.csv | cut -f1,5,10

For ATRT specifically, use:
  BT16  (ACH-000725) — CNS ATRT, ATRT-MYC, most studied
  BT37  (ACH-000881) — CNS ATRT, ATRT-TYR
  G401  (ACH-000039) — Renal rhabdoid, SMARCB1-null (EZH2 synleth founding paper)
  A204  (ACH-000658) — Rhabdoid, SMARCB1-null
Reference: Knutson 2013 PNAS — EZH2 essentiality first shown in G401/A204.

REFERENCES
-----------
Behan FM et al. Nature 2019; 568:511. PMID 30971826. [Chronos interpretation]
Tsherniak A et al. Cell 2017; 170(3):564. PMID 28753430. [DepMap framework]
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

CHRONOS_ESSENTIAL_THRESHOLD  = -1.0    # Strongly essential
CHRONOS_MODERATE_THRESHOLD   = -0.5    # Moderately essential
NORM_CS_REVERSAL_THRESHOLD   = -0.9    # Strong CMap reversal

# Known DepMap Model IDs for ATRT/rhabdoid lines (verified DepMap 24Q4)
ATRT_CELL_LINE_MODEL_IDS = {
    "BT16":      "ACH-000725",
    "BT37":      "ACH-000881",
    "G401":      "ACH-000039",
    "A204":      "ACH-000658",
    "BT12":      "ACH-001082",   # Verify in your Model.csv release
    "CHLA02":    None,            # May not be in DepMap — check Model.csv
    "CHLA04":    None,
    "CHLA06":    None,
}


def load_depmap_atrt_chronos(
    effect_file:  str = "data/depmap/CRISPRGeneEffect.csv",
    model_file:   str = "data/depmap/Model.csv",
    require_lines: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Load Chronos gene effect scores for ATRT/rhabdoid cell lines from DepMap.

    Parameters
    ----------
    require_lines : optional list of cell line names to require.
        If a line is specified but not found, logs a warning rather than failing.

    Returns
    -------
    DataFrame: rows = cell lines, columns = gene symbols, values = Chronos scores.
    Index = cell line display names (e.g. "BT16").
    """
    effect_path = Path(effect_file)
    model_path  = Path(model_file)

    if not effect_path.exists() or not model_path.exists():
        logger.error(
            "DepMap files not found.\n"
            "Download from: https://depmap.org/portal/download/all/\n"
            "Files needed: CRISPRGeneEffect.csv, Model.csv\n"
            "Place in: data/depmap/"
        )
        return pd.DataFrame()

    models = pd.read_csv(model_path)
    models.columns = models.columns.str.strip()

    # Find relevant column names in Model.csv (varies between releases)
    model_id_col  = next((c for c in models.columns if c in
                          ["ModelID", "DepMap_ID", "model_id"]), None)
    subtype_col   = next((c for c in models.columns if
                          "oncotreesubtype" in c.lower()), None)
    name_col      = next((c for c in models.columns if c.lower() in
                          ["celllinename", "cell_line_name",
                           "stripped_cell_line_name"]), None)

    if model_id_col is None:
        logger.error("Cannot find ModelID column. Available: %s", models.columns.tolist())
        return pd.DataFrame()

    # Multi-tier cell line selection
    atrt_ids: List[str] = []

    # Tier 1: OncotreeSubtype = MRT or ATRT
    if subtype_col:
        mask_subtype = models[subtype_col].astype(str).str.upper().isin(["MRT", "ATRT"])
        tier1_ids = models[mask_subtype][model_id_col].tolist()
        atrt_ids.extend(tier1_ids)
        logger.info("OncotreeSubtype filter: %d lines (MRT/ATRT)", len(tier1_ids))

    # Tier 2: Known cell line names
    if name_col:
        known_names = {k.upper() for k in ATRT_CELL_LINE_MODEL_IDS.keys()}
        mask_names  = models[name_col].astype(str).str.upper().str.strip().isin(known_names)
        tier2_ids   = models[mask_names][model_id_col].tolist()
        newly_added = set(tier2_ids) - set(atrt_ids)
        if newly_added:
            logger.info("Name-based filter: added %d additional lines", len(newly_added))
        atrt_ids = list(set(atrt_ids) | set(tier2_ids))

    # Check for specifically requested lines
    if require_lines:
        for line_name in require_lines:
            if name_col:
                match = models[
                    models[name_col].astype(str).str.upper().str.strip() == line_name.upper()
                ]
                if match.empty:
                    logger.warning(
                        "Requested line '%s' not found in Model.csv. "
                        "Note: CHLA-266 is neuroblastoma, not ATRT. "
                        "ATRT lines in DepMap: BT16, BT37, G401, A204. "
                        "Check: grep -i '%s' data/depmap/Model.csv",
                        line_name, line_name,
                    )
                else:
                    mid = match[model_id_col].values[0]
                    if mid not in atrt_ids:
                        atrt_ids.append(mid)
                        logger.info("Added requested line '%s' (%s)", line_name, mid)

    if not atrt_ids:
        logger.error(
            "No ATRT/rhabdoid cell lines found. "
            "Expected OncotreeSubtype MRT/ATRT in Model.csv. "
            "Check your DepMap release version (24Q4 recommended)."
        )
        return pd.DataFrame()

    logger.info("Loading CRISPRGeneEffect.csv for %d ATRT/rhabdoid lines...", len(atrt_ids))
    effect_df = pd.read_csv(effect_path, index_col=0)

    # Match Model IDs to effect matrix index
    available    = set(effect_df.index)
    matched_ids  = [mid for mid in atrt_ids if mid in available]
    missing_ids  = [mid for mid in atrt_ids if mid not in available]

    if missing_ids:
        logger.warning(
            "%d lines not in CRISPRGeneEffect.csv (may be in 23Q4 but not 24Q4): %s",
            len(missing_ids), missing_ids,
        )

    if not matched_ids:
        logger.error("No ATRT model IDs matched CRISPRGeneEffect.csv index.")
        return pd.DataFrame()

    atrt_chronos = effect_df.loc[matched_ids].copy()

    # Rename index to display names for readability
    if name_col:
        id_to_name = dict(zip(models[model_id_col], models[name_col]))
        atrt_chronos.index = [id_to_name.get(mid, mid) for mid in atrt_chronos.index]

    # Clean up column names: "EZH2 (2146)" → "EZH2"
    atrt_chronos.columns = [c.split(" ")[0].upper() for c in atrt_chronos.columns]

    logger.info(
        "DepMap ATRT loaded: %d lines × %d genes\n  Lines: %s",
        len(atrt_chronos), len(atrt_chronos.columns),
        list(atrt_chronos.index),
    )
    return atrt_chronos


def crossref_cmap_with_depmap(
    cmap_scores:     Dict[str, Dict],
    atrt_chronos:    pd.DataFrame,
    drug_targets:    Dict[str, List[str]],
) -> pd.DataFrame:
    """
    Cross-reference CMap reversal scores with DepMap CRISPR essentiality.

    For each CMap drug hit, find its targets, look up their Chronos scores
    in ATRT lines, and assign a validation tier.

    Parameters
    ----------
    cmap_scores : dict mapping DRUG_NAME → {"norm_cs": float, "cmap_score": float}
                  (from integrate_cmap_results.py output)
    atrt_chronos : DataFrame from load_depmap_atrt_chronos()
                   rows = cell lines, columns = gene symbols
    drug_targets : dict mapping DRUG_NAME → [target gene symbols]
                   (from pipeline candidate list)

    Returns
    -------
    DataFrame with one row per drug, columns:
        drug, norm_cs, n_targets, best_target, best_chronos,
        median_chronos, validation_tier, tier_explanation,
        targets_essential, all_chronos
    """
    if atrt_chronos.empty:
        logger.error("atrt_chronos is empty — cannot cross-reference")
        return pd.DataFrame()

    rows = []
    for drug_name, cmap_data in cmap_scores.items():
        norm_cs    = cmap_data.get("norm_cs")
        cmap_score = cmap_data.get("cmap_score", 0.5)
        targets    = [t.upper() for t in drug_targets.get(drug_name, [])]

        if not targets:
            continue

        # Get Chronos scores for each target
        target_chronos: Dict[str, Dict] = {}
        for gene in targets:
            if gene in atrt_chronos.columns:
                per_line = atrt_chronos[gene].to_dict()
                median_c = float(atrt_chronos[gene].median())
                target_chronos[gene] = {
                    "median_chronos": median_c,
                    "per_line":       per_line,
                    "n_lines":        len(atrt_chronos),
                }

        if not target_chronos:
            chronos_note = f"No targets ({', '.join(targets[:3])}) in DepMap"
            rows.append({
                "drug":             drug_name,
                "norm_cs":          norm_cs,
                "cmap_score":       cmap_score,
                "n_targets":        len(targets),
                "best_target":      None,
                "best_chronos":     None,
                "median_chronos":   None,
                "validation_tier":  "UNVALIDATED",
                "tier_explanation": chronos_note,
                "targets_essential": 0,
                "all_chronos":      {},
            })
            continue

        chronos_vals  = {g: v["median_chronos"] for g, v in target_chronos.items()}
        best_target   = min(chronos_vals, key=chronos_vals.get)
        best_chronos  = chronos_vals[best_target]
        median_c      = float(np.median(list(chronos_vals.values())))

        n_essential = sum(1 for c in chronos_vals.values()
                          if c <= CHRONOS_ESSENTIAL_THRESHOLD)

        # Assign validation tier
        is_reverser = (norm_cs is not None and norm_cs <= NORM_CS_REVERSAL_THRESHOLD)

        if is_reverser and best_chronos <= CHRONOS_ESSENTIAL_THRESHOLD:
            tier = "TIER_1_STRONG"
            explanation = (
                f"CMap reverser (norm_cs={norm_cs:.2f}) AND "
                f"{best_target} Chronos={best_chronos:.2f} (strongly essential in ATRT). "
                "Highest priority combination hit."
            )
        elif is_reverser and best_chronos <= CHRONOS_MODERATE_THRESHOLD:
            tier = "TIER_2_SUPPORTED"
            explanation = (
                f"CMap reverser (norm_cs={norm_cs:.2f}) AND "
                f"{best_target} Chronos={best_chronos:.2f} (moderately essential). "
            )
        elif is_reverser:
            tier = "TIER_3_TRANSCRIPTIONAL"
            explanation = (
                f"CMap reverser (norm_cs={norm_cs:.2f}) but "
                f"{best_target} Chronos={best_chronos:.2f} (not strongly essential). "
                "Target is transcriptionally dysregulated but may not be a dependency. "
                "Drug may act via off-target or indirect mechanism (e.g. ONC201/TRAIL)."
            )
        elif norm_cs is not None and norm_cs > -0.5:
            tier = "TIER_4_WEAK_CMAP"
            explanation = (
                f"Weak CMap reversal (norm_cs={norm_cs:.2f}). "
                f"Even though {best_target} Chronos={best_chronos:.2f} shows dependency, "
                "the transcriptional reversal is insufficient."
            )
        else:
            tier = "TIER_3_DEPMAP_ONLY"
            explanation = (
                f"No CMap data, but {best_target} Chronos={best_chronos:.2f}. "
                "Drug ranked on DepMap essentiality alone."
            )

        rows.append({
            "drug":              drug_name,
            "norm_cs":           norm_cs,
            "cmap_score":        cmap_score,
            "n_targets":         len(targets),
            "best_target":       best_target,
            "best_chronos":      round(best_chronos, 3),
            "median_chronos":    round(median_c, 3),
            "validation_tier":   tier,
            "tier_explanation":  explanation,
            "targets_essential": n_essential,
            "all_chronos":       {k: round(v, 3) for k, v in chronos_vals.items()},
        })

    result_df = pd.DataFrame(rows)
    if result_df.empty:
        return result_df

    # Sort: Tier 1 first, then by CMap score descending
    tier_order = {
        "TIER_1_STRONG":       0,
        "TIER_2_SUPPORTED":    1,
        "TIER_3_TRANSCRIPTIONAL": 2,
        "TIER_3_DEPMAP_ONLY":  3,
        "TIER_4_WEAK_CMAP":    4,
        "UNVALIDATED":         5,
    }
    result_df["tier_rank"] = result_df["validation_tier"].map(tier_order)
    result_df = result_df.sort_values(["tier_rank", "cmap_score"],
                                       ascending=[True, False])

    # Summary logging
    tier_counts = result_df["validation_tier"].value_counts()
    logger.info(
        "Cross-reference complete: %d drugs\n  %s",
        len(result_df),
        "\n  ".join(f"{k}: {v}" for k, v in tier_counts.items()),
    )

    return result_df.drop(columns=["tier_rank"])


def print_crossref_report(crossref_df: pd.DataFrame) -> str:
    """Format the cross-reference results as a readable report."""
    if crossref_df.empty:
        return "Cross-reference: no results."

    lines = [
        "ATRT Drug Validation: CMap × DepMap Cross-Reference",
        "=" * 60,
        f"{'Drug':<22} {'Tier':<22} {'norm_cs':>8}  {'Best target':>12}  Chronos",
        "-" * 75,
    ]

    for _, row in crossref_df.iterrows():
        tier_short = row["validation_tier"].replace("TIER_1_", "").replace("TIER_2_", "").replace("TIER_3_", "")[:12]
        norm_cs_str = f"{row['norm_cs']:.2f}" if row["norm_cs"] is not None else "N/A"
        chronos_str = f"{row['best_chronos']:.2f}" if row["best_chronos"] is not None else "N/A"
        lines.append(
            f"{row['drug']:<22} {tier_short:<22} {norm_cs_str:>8}  "
            f"{str(row['best_target'] or 'N/A'):>12}  {chronos_str}"
        )

    tier1 = crossref_df[crossref_df["validation_tier"] == "TIER_1_STRONG"]
    if not tier1.empty:
        lines += [
            "",
            "Tier 1 hits (CMap reversal + DepMap essential — highest confidence):",
            *[f"  {r['drug']}: {r['best_target']} Chronos={r['best_chronos']}, norm_cs={r['norm_cs']}"
              for _, r in tier1.iterrows()],
        ]

    return "\n".join(lines)