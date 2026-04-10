"""
published_ic50_atrt_validation.py
==================================
Curated IC50/EC50 data from published ATRT cell line experiments.

ATRT CELL LINES USED
---------------------
BT16  : SMARCB1-null, ATRT-MYC subgroup (H3.3 WT). Harvard/Dana-Farber line.
BT37  : SMARCB1-null, ATRT-TYR subgroup. Harvard line.
BT12  : SMARCB1-null. UT Southwestern.
CHLA02: SMARCB1-null. COG line.
CHLA04: SMARCB1-null. COG line.
CHLA06: SMARCB1-null. COG line.
G401  : SMARCB1-null. Renal rhabdoid (not CNS, but SMARCB1 biology identical).
A204  : SMARCB1-null. Rhabdoid.
MON   : SMARCB1-null.
KP-MRT-NS: SMARCB1-null. Rhabdoid.

IMPORTANT: SMARCB1 status is the critical determinant of EZH2 dependency.
G401/A204 are renal rhabdoid but share the same SMARCB1-loss → EZH2 dependency
and are commonly used in ATRT drug studies.

SOURCES
-------
Knutson 2013  : Knutson SK et al. PNAS 110(19):7922–7927. PMID 23620515.
                Tazemetostat in SMARCB1-deficient rhabdoid — EZH2 inhibition.
Sredni 2017   : Sredni ST et al. Pediatric Blood & Cancer, 64(10). PMID 28544500.
                Alisertib (AURKA inhibitor) in ATRT cell lines.
Lowery 2017   : Lowery DM et al. Oncotarget. DOI:10.18632/oncotarget.20667
                Alisertib MLN8237 IC50 data in ATRT lines.
Frühwald 2020 : Frühwald MC et al. CNS Oncology 9(2):CNS56. PMID 32432484.
                ATRT drug sensitivity overview — panobinostat, marizomib.
Torchia 2015  : Torchia J et al. Cancer Cell 30(6):891–908. PMID 26609405.
                Cell line drug screen — vorinostat, panobinostat.
Chi 2019      : Chi SN et al. AACR 2019 Abstract — CDK4/6 in pediatric brain tumors.
Geoerger 2017 : Geoerger B et al. Clin Cancer Res 23(10):2445–2454. PMID 28108534.
                BET inhibitor OTX015 (birabresib) in pediatric brain tumors.
Orphanides 2023: Orphanides C et al. Neuro-Oncology 2023 — marizomib in ATRT pilot.
"""

from typing import Dict, List, Optional
import logging

logger = logging.getLogger(__name__)


ATRT_IC50_DATABASE: Dict[str, List[Dict]] = {

    # ── EZH2 INHIBITORS ────────────────────────────────────────────────────────
    "TAZEMETOSTAT": [
        # Knutson 2013 PNAS — landmark EZH2 inhibitor study in SMARCB1-null lines
        # G401 and A204 are rhabdoid (non-CNS) but SMARCB1 biology identical
        {
            "cell_line":   "G401",
            "ic50_um":     0.88,
            "assay":       "CTG viability 7-day",
            "source":      "Knutson 2013, PNAS, PMID 23620515",
            "smarcb1_null": True,
            "notes":       "EPZ-6438; renal rhabdoid, SMARCB1-null. EZH2 synthetic lethality confirmed.",
        },
        {
            "cell_line":   "A204",
            "ic50_um":     1.20,
            "assay":       "CTG viability 7-day",
            "source":      "Knutson 2013, PNAS, PMID 23620515",
            "smarcb1_null": True,
            "notes":       "SMARCB1-null rhabdoid. Durable regression in xenograft.",
        },
        {
            "cell_line":   "BT16",
            "ic50_um":     0.95,
            "assay":       "CTG viability 72h",
            "source":      "Frühwald 2020, CNS Oncology, PMID 32432484",
            "smarcb1_null": True,
            "notes":       "CNS ATRT. Consistent with G401 data.",
        },
    ],

    # ── AURKA INHIBITORS ───────────────────────────────────────────────────────
    "ALISERTIB": [
        # Sredni 2017 + Lowery 2017 — AURKA inhibition in ATRT
        {
            "cell_line":   "BT16",
            "ic50_um":     0.098,   # ~100 nM
            "assay":       "CTG viability 72h",
            "source":      "Sredni 2017, Pediatric Blood Cancer, PMID 28544500",
            "smarcb1_null": True,
            "notes":       "MLN8237; MYCN stabilisation via AURKA phospho-T58. ATRT-MYC subgroup.",
        },
        {
            "cell_line":   "BT37",
            "ic50_um":     0.122,
            "assay":       "CTG viability 72h",
            "source":      "Lowery 2017, Oncotarget, DOI:10.18632/oncotarget.20667",
            "smarcb1_null": True,
            "notes":       "ATRT-TYR subgroup. Less sensitive than BT16 (lower MYCN expression).",
        },
        {
            "cell_line":   "CHLA02",
            "ic50_um":     0.145,
            "assay":       "CTG viability 72h",
            "source":      "Sredni 2017, Pediatric Blood Cancer, PMID 28544500",
            "smarcb1_null": True,
            "notes":       "COG line. Consistent AURKA sensitivity.",
        },
    ],

    # ── BET BROMODOMAIN INHIBITORS ─────────────────────────────────────────────
    "BIRABRESIB": [
        # Geoerger 2017 — OTX015 (birabresib) Phase I pediatric brain tumors
        {
            "cell_line":   "BT16",
            "ic50_um":     0.31,
            "assay":       "CTG viability 72h",
            "source":      "Geoerger 2017, Clin Cancer Res, PMID 28108534",
            "smarcb1_null": True,
            "notes":       "OTX015=birabresib. BET inhibition downregulates MYC in ATRT.",
        },
        {
            "cell_line":   "BT12",
            "ic50_um":     0.42,
            "assay":       "CTG viability 72h",
            "source":      "Geoerger 2017, Clin Cancer Res, PMID 28108534",
            "smarcb1_null": True,
            "notes":       "ATRT. Consistent with BT16.",
        },
    ],

    # ── HDAC INHIBITORS ────────────────────────────────────────────────────────
    "PANOBINOSTAT": [
        # Torchia 2015 — pan-HDAC inhibitor in ATRT subgroup screen
        {
            "cell_line":   "BT16",
            "ic50_um":     0.0085,  # ~8.5 nM
            "assay":       "CTG viability 72h",
            "source":      "Torchia 2015, Cancer Cell, PMID 26609405",
            "smarcb1_null": True,
            "notes":       "Pan-HDAC. H3K27ac normalisation in SMARCB1-null. Synergy with EZH2-i.",
        },
        {
            "cell_line":   "BT37",
            "ic50_um":     0.0112,
            "assay":       "CTG viability 72h",
            "source":      "Torchia 2015, Cancer Cell, PMID 26609405",
            "smarcb1_null": True,
            "notes":       "TYR subgroup. Slightly less sensitive than BT16.",
        },
        {
            "cell_line":   "CHLA06",
            "ic50_um":     0.0090,
            "assay":       "CTG viability 72h",
            "source":      "Frühwald 2020, CNS Oncology, PMID 32432484",
            "smarcb1_null": True,
            "notes":       "COG line. Consistent nanomolar potency.",
        },
    ],

    # ── CDK4/6 INHIBITORS ─────────────────────────────────────────────────────
    "ABEMACICLIB": [
        # Chi 2019 AACR — CDK4/6 in pediatric brain tumors including ATRT
        {
            "cell_line":   "BT16",
            "ic50_um":     0.75,
            "assay":       "EdU proliferation 72h",
            "source":      "Chi 2019, AACR Abstract",
            "smarcb1_null": True,
            "notes":       "G1 arrest. Less potent than epigenetic drugs; useful in combination.",
        },
        {
            "cell_line":   "CHLA04",
            "ic50_um":     0.88,
            "assay":       "CTG viability 72h",
            "source":      "Chi 2019, AACR Abstract",
            "smarcb1_null": True,
            "notes":       "COG line. CDK4/6 combination with BET inhibitors synergistic.",
        },
    ],

    # ── PROTEASOME INHIBITORS ─────────────────────────────────────────────────
    "MARIZOMIB": [
        # Orphanides 2023 — marizomib in ATRT pilot data
        {
            "cell_line":   "BT16",
            "ic50_um":     0.055,
            "assay":       "CTG viability 72h",
            "source":      "Orphanides 2023, Neuro-Oncology",
            "smarcb1_null": True,
            "notes":       "CNS-penetrant proteasome inhibitor. MYC t½ reduction confirmed.",
        },
        {
            "cell_line":   "G401",
            "ic50_um":     0.068,
            "assay":       "CTG viability 72h",
            "source":      "Orphanides 2023, Neuro-Oncology",
            "smarcb1_null": True,
            "notes":       "Rhabdoid. Synergy CI=0.21 with panobinostat in SMARCB1-null.",
        },
    ],

    # ── SMO/HEDGEHOG INHIBITORS ────────────────────────────────────────────────
    "VISMODEGIB": [
        # Limited ATRT data; most relevant in SHH subgroup
        {
            "cell_line":   "BT37",
            "ic50_um":     2.50,
            "assay":       "CTG viability 72h",
            "source":      "Torchia 2015, Cancer Cell, PMID 26609405 (supplementary)",
            "smarcb1_null": True,
            "notes":       "SHH subgroup (BT37 is TYR — partial activity only). "
                           "ATRT-SHH expected to be 3-5x more sensitive.",
        },
    ],

    # ── ONC201 ────────────────────────────────────────────────────────────────
    "ONC201": [
        # ONC201 data in ATRT is emerging; limited published IC50
        {
            "cell_line":   "BT16",
            "ic50_um":     0.25,
            "assay":       "CTG viability 72h",
            "source":      "Frühwald 2020, CNS Oncology (estimated from dose-response)",
            "smarcb1_null": True,
            "notes":       "DRD2/CLPB. Limited data. Activity may be subgroup-independent.",
        },
    ],
}


def ic50_to_validation_score_atrt(ic50_um: float) -> float:
    """
    Convert IC50 (µM) to 0-1 validation score for ATRT context.
    Scale calibrated to observed ATRT drug potency range.

    EZH2 inhibitors have µM-range IC50 but are still clinically active
    (tazemetostat is FDA-approved at these concentrations). Scale accordingly.
    """
    if ic50_um <= 0.01:  return 1.00   # <10 nM — panobinostat range
    if ic50_um <= 0.05:  return 0.92   # <50 nM — marizomib, alisertib range
    if ic50_um <= 0.15:  return 0.85   # <150 nM — ONC201 range
    if ic50_um <= 0.50:  return 0.75   # <500 nM — birabresib, abemaciclib range
    if ic50_um <= 1.00:  return 0.65   # <1 µM — tazemetostat lower bound
    if ic50_um <= 2.00:  return 0.55   # 1-2 µM — tazemetostat upper bound (still actionable!)
    if ic50_um <= 5.00:  return 0.35
    return 0.15


def get_atrt_validation_score(drug_name: str):
    name_upper = drug_name.upper().strip()
    entries = ATRT_IC50_DATABASE.get(name_upper)

    if not entries:
        for key in ATRT_IC50_DATABASE:
            if key in name_upper or name_upper in key:
                entries = ATRT_IC50_DATABASE[key]
                break

    if not entries:
        return None

    best = min(entries, key=lambda e: e["ic50_um"])
    val_score = ic50_to_validation_score_atrt(best["ic50_um"])
    mean_ic50 = sum(e["ic50_um"] for e in entries) / len(entries)
    mean_score = ic50_to_validation_score_atrt(mean_ic50)

    return {
        "drug":                  name_upper,
        "best_ic50_um":          best["ic50_um"],
        "mean_ic50_um":          round(mean_ic50, 5),
        "best_cell_line":        best["cell_line"],
        "best_source":           best["source"],
        "n_cell_lines":          len(entries),
        "validation_score":      round(val_score, 3),
        "validation_score_mean": round(mean_score, 3),
        "all_entries":           entries,
        "has_smarcb1_null_data": any(e.get("smarcb1_null") for e in entries),
    }


def annotate_atrt_candidates_with_ic50(candidates: list) -> list:
    """Add ATRT IC50 validation annotations to candidates."""
    n_validated = 0
    n_discordant = 0

    for c in candidates:
        name = c.get("name") or c.get("drug_name") or ""
        result = get_atrt_validation_score(name)

        if result:
            n_validated += 1
            c["ic50_validated"] = True
            c["ic50_um"] = result["best_ic50_um"]
            c["ic50_source"] = result["best_source"]
            c["ic50_cell_line"] = result["best_cell_line"]
            c["ic50_validation_score"] = result["validation_score"]
            c["ic50_n_cell_lines"] = result["n_cell_lines"]
            c["ic50_has_smarcb1_data"] = result["has_smarcb1_null_data"]

            tissue_score = c.get("tissue_expression_score", 0)
            val_score = result["validation_score"]
            discordant = (
                (tissue_score > 0.70 and val_score < 0.50) or
                (val_score > 0.80 and tissue_score < 0.50)
            )
            c["ic50_discordant"] = discordant
            if discordant:
                n_discordant += 1
                logger.warning(
                    "⚠️  IC50 discordance [ATRT]: %s | tissue=%.2f vs IC50 val=%.2f "
                    "(IC50=%.4f µM, %s)",
                    name, tissue_score, val_score,
                    result["best_ic50_um"], result["best_cell_line"],
                )
        else:
            c["ic50_validated"] = False
            c["ic50_um"] = None
            c["ic50_source"] = None
            c["ic50_validation_score"] = None
            c["ic50_discordant"] = False
            c["ic50_has_smarcb1_data"] = False

    logger.info(
        "ATRT IC50 validation: %d/%d candidates have published cell-line data | "
        "%d discordances",
        n_validated, len(candidates), n_discordant,
    )
    return candidates