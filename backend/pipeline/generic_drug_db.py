"""
generic_drug_db.py
==================
Curated database of FDA-approved generic and repurposable drugs for ATRT.

WHAT "GENERIC DRUGS" MEANS HERE
---------------------------------
Three categories are included:

1. GENERIC EQUIVALENTS of key brand-name drugs
   e.g., vorinostat (Zolinza) is a generic HDAC inhibitor available as vorinostat.
   Panobinostat is still mostly brand-name (Farydak) but is the compound name.

2. REPURPOSABLE APPROVED GENERICS — drugs approved for other indications
   that have mechanistic rationale in ATRT. These are typically cheap,
   available, and have known safety profiles in children.
   e.g., valproic acid (HDAC inhibitor), metformin (AMPK), ATRA (differentiation).

3. TOOL COMPOUNDS — extensively validated in ATRT cell lines.
   e.g., EPZ-6438 = tazemetostat (same compound, development name).

WHY INCLUDE GENERICS?
-----------------------
ATRT is a fatal pediatric disease affecting 30 new US patients/year.
Branded drugs are often inaccessible. Generic/repurposable drugs offer:
  - Established paediatric safety data
  - Immediate availability for clinical trials
  - Lower cost for families and trial sponsors
  - Known pharmacokinetics (no Phase I needed for dose-finding)

EVIDENCE GRADING
-----------------
Each drug is tagged with an evidence level:
  HIGH     — published ATRT/rhabdoid cell line IC50 or in vivo data
  MODERATE — published activity in related SWI/SNF-deficient tumours
  LOW      — mechanistic rationale only; no primary ATRT data

References
-----------
Frühwald MC et al. (2020). ATRT — current biology. CNS Oncology 9(2):CNS56.
Torchia J et al. (2015). Cancer Cell 30(6):891-908. PMID 26609405.
Knutson SK et al. (2013). PNAS 110(19):7922. PMID 23620515.
Göttlicher M et al. (2001). EMBO J 20(24):6969 — Valproic acid HDAC.
Shim JS & Liu JO (2014). Front Oncol 4:75 — Drug repurposing overview.
"""

from typing import Dict, List, Optional

# ─────────────────────────────────────────────────────────────────────────────
# Generic / repurposable drug catalogue
# ─────────────────────────────────────────────────────────────────────────────

GENERIC_DRUG_DATABASE: List[Dict] = [

    # ════════════════════════════════════════════════════════════════════════
    # HDAC INHIBITORS (generic equivalents + repurposable)
    # ════════════════════════════════════════════════════════════════════════

    {
        "name":              "Vorinostat",
        "generic_name":      "vorinostat",
        "brand_names":       ["Zolinza"],
        "targets":           ["HDAC1", "HDAC2", "HDAC3", "HDAC6"],
        "mechanism":         "pan-hdac inhibitor",
        "drug_class":        "HDAC inhibitor",
        "fda_approved":      True,
        "approved_indication": "Cutaneous T-cell lymphoma",
        "is_generic_available": False,   # Still brand-name
        "bbb_penetrance":    "MODERATE",
        "evidence_level":    "MODERATE",
        "atrt_rationale": (
            "Pan-HDAC inhibition reverses aberrant chromatin compaction caused by "
            "SMARCB1 loss. Less potent than panobinostat in BT16/BT37 but approved "
            "and has known paediatric PK (PBTC trial). Galanis 2009."
        ),
        "smarcb1_relevant":  True,
        "cost_category":     "BRAND",
        "source":            "FDA Orange Book + PBTC trial experience",
    },

    {
        "name":              "Valproic acid",
        "generic_name":      "valproic acid",
        "brand_names":       ["Depakote", "Depakene", "Convulex"],
        "targets":           ["HDAC1", "HDAC2", "HDAC3"],
        "mechanism":         "pan-hdac inhibitor class I/II",
        "drug_class":        "HDAC inhibitor (repurposed anticonvulsant)",
        "fda_approved":      True,
        "approved_indication": "Epilepsy, bipolar disorder, migraine",
        "is_generic_available": True,
        "bbb_penetrance":    "HIGH",
        "evidence_level":    "MODERATE",
        "atrt_rationale": (
            "Established paediatric HDAC inhibitor. Crosses BBB. Well-tolerated "
            "at anticonvulsant doses (50-100 µg/mL serum). Göttlicher 2001 EMBO J "
            "(PMID 11742990). Relevant for ATRT due to SMARCB1-loss chromatin context. "
            "Readily available as generic; used in children with epilepsy."
        ),
        "smarcb1_relevant":  True,
        "cost_category":     "GENERIC",
        "ic50_um":           None,
        "ic50_note":         "IC50 ~0.4 mM (anticonvulsant class; µM vs nM for other HDAC-i)",
        "source":            "Göttlicher 2001 EMBO J PMID 11742990",
    },

    {
        "name":              "Entinostat",
        "generic_name":      "entinostat",
        "brand_names":       ["Syndax"],
        "targets":           ["HDAC1", "HDAC2", "HDAC3"],
        "mechanism":         "class I hdac inhibitor",
        "drug_class":        "Class I HDAC inhibitor",
        "fda_approved":      False,   # NDA filed
        "approved_indication": "Breast cancer (NDA under review)",
        "is_generic_available": False,
        "bbb_penetrance":    "MODERATE",
        "evidence_level":    "LOW",
        "atrt_rationale": (
            "Class I selective HDAC inhibition specifically targets HDAC1/2/3 "
            "which are the ATRT-relevant HDAC subunits (Torchia 2015). "
            "Selective class I profile may reduce haematologic toxicity vs pan-HDAC."
        ),
        "smarcb1_relevant":  True,
        "cost_category":     "BRAND",
        "source":            "Syndax NDA; mechanistic rationale",
    },

    {
        "name":              "Romidepsin",
        "generic_name":      "romidepsin",
        "brand_names":       ["Istodax"],
        "targets":           ["HDAC1", "HDAC2"],
        "mechanism":         "class I hdac inhibitor (bicyclic depsipeptide)",
        "drug_class":        "Class I HDAC inhibitor",
        "fda_approved":      True,
        "approved_indication": "Cutaneous/peripheral T-cell lymphoma",
        "is_generic_available": False,
        "bbb_penetrance":    "LOW",
        "evidence_level":    "LOW",
        "atrt_rationale":    "Class I HDAC inhibitor; poor BBB limits ATRT utility.",
        "smarcb1_relevant":  True,
        "cost_category":     "BRAND",
        "source":            "Mechanistic rationale only",
    },

    # ════════════════════════════════════════════════════════════════════════
    # EZH2 INHIBITORS (SMARCB1 synthetic lethality — core ATRT target)
    # ════════════════════════════════════════════════════════════════════════

    {
        "name":              "GSK126",
        "generic_name":      "gsk126",
        "brand_names":       [],
        "targets":           ["EZH2"],
        "mechanism":         "ezh2 inhibitor",
        "drug_class":        "EZH2 inhibitor (tool compound)",
        "fda_approved":      False,
        "approved_indication": "Preclinical only",
        "is_generic_available": False,
        "bbb_penetrance":    "UNKNOWN",
        "evidence_level":    "HIGH",
        "atrt_rationale": (
            "First validated EZH2 inhibitor used to establish synthetic lethality "
            "in SMARCB1-null G401/A204 cells. IC50 ~9.9 nM vs EZH2. "
            "Knutson 2013 PNAS PMID 23620515 — the founding study. "
            "Not clinically approved but validates the target; tazemetostat "
            "is the clinical-stage analogue."
        ),
        "smarcb1_relevant":  True,
        "cost_category":     "TOOL_COMPOUND",
        "ic50_um":           0.038,
        "ic50_cell_line":    "G401",
        "source":            "Knutson 2013 PNAS PMID 23620515",
    },

    {
        "name":              "EPZ-6438",
        "generic_name":      "tazemetostat",   # same compound, development name
        "brand_names":       ["Tazverik"],
        "targets":           ["EZH2"],
        "mechanism":         "ezh2 inhibitor",
        "drug_class":        "EZH2 inhibitor",
        "fda_approved":      True,
        "approved_indication": "Epithelioid sarcoma (EZH2 wt/mut), follicular lymphoma",
        "is_generic_available": False,
        "bbb_penetrance":    "MODERATE",
        "evidence_level":    "HIGH",
        "atrt_rationale": (
            "EPZ-6438 = tazemetostat (same molecule, development code). "
            "FDA Breakthrough Therapy Designation for SMARCB1-deficient tumours. "
            "Knutson 2013 PNAS; Gounder 2020 JCO."
        ),
        "smarcb1_relevant":  True,
        "cost_category":     "BRAND",
        "source":            "Knutson 2013 PNAS; FDA BT Designation",
    },

    # ════════════════════════════════════════════════════════════════════════
    # BET BROMODOMAIN INHIBITORS
    # ════════════════════════════════════════════════════════════════════════

    {
        "name":              "JQ1",
        "generic_name":      "jq1",
        "brand_names":       [],
        "targets":           ["BRD4", "BRD2", "BRD3", "BRDT"],
        "mechanism":         "bet bromodomain inhibitor",
        "drug_class":        "BET bromodomain inhibitor (tool compound)",
        "fda_approved":      False,
        "approved_indication": "Preclinical only",
        "is_generic_available": False,
        "bbb_penetrance":    "MODERATE",
        "evidence_level":    "HIGH",
        "atrt_rationale": (
            "Pan-BET inhibitor; extensively validated in ATRT cell lines. "
            "IC50 ~50-100 nM in BT16. Clinical successor is birabresib (OTX015). "
            "Hennika 2017 PLOS ONE. Important as a research tool and class validator."
        ),
        "smarcb1_relevant":  True,
        "cost_category":     "TOOL_COMPOUND",
        "ic50_um":           0.08,
        "ic50_cell_line":    "BT16",
        "source":            "Hennika 2017 PLOS ONE PMID 28056098",
    },

    {
        "name":              "Molibresib",
        "generic_name":      "molibresib",
        "brand_names":       ["GSK525762"],
        "targets":           ["BRD4", "BRD2", "BRD3"],
        "mechanism":         "bet bromodomain inhibitor",
        "drug_class":        "BET bromodomain inhibitor",
        "fda_approved":      False,
        "approved_indication": "Phase I/II solid tumours",
        "is_generic_available": False,
        "bbb_penetrance":    "MODERATE",
        "evidence_level":    "LOW",
        "atrt_rationale":    "Pan-BET; same class as birabresib. Limited ATRT data.",
        "smarcb1_relevant":  True,
        "cost_category":     "INVESTIGATIONAL",
        "source":            "GSK development; class rationale",
    },

    # ════════════════════════════════════════════════════════════════════════
    # AURORA KINASE INHIBITORS
    # ════════════════════════════════════════════════════════════════════════

    {
        "name":              "Barasertib",
        "generic_name":      "barasertib",
        "brand_names":       ["AZD1152"],
        "targets":           ["AURKB"],
        "mechanism":         "aurora kinase b inhibitor",
        "drug_class":        "Aurora kinase B inhibitor",
        "fda_approved":      False,
        "approved_indication": "Phase II AML",
        "is_generic_available": False,
        "bbb_penetrance":    "MODERATE",
        "evidence_level":    "LOW",
        "atrt_rationale": (
            "AURKB inhibition. AURKB is upregulated in ATRT (GSE70678 z-score 1.35). "
            "Less ATRT-specific than alisertib (AURKA + MYCN stabilisation mechanism). "
            "Some CNS tumour data. Sredni 2017 includes AURKB context."
        ),
        "smarcb1_relevant":  True,
        "cost_category":     "INVESTIGATIONAL",
        "source":            "Sredni 2017 PMID 28544500 — context",
    },

    # ════════════════════════════════════════════════════════════════════════
    # CDK INHIBITORS
    # ════════════════════════════════════════════════════════════════════════

    {
        "name":              "Palbociclib",
        "generic_name":      "palbociclib",
        "brand_names":       ["Ibrance"],
        "targets":           ["CDK4", "CDK6"],
        "mechanism":         "cdk4/6 inhibitor",
        "drug_class":        "CDK4/6 inhibitor",
        "fda_approved":      True,
        "approved_indication": "HR+ HER2- breast cancer",
        "is_generic_available": False,
        "bbb_penetrance":    "LOW",
        "evidence_level":    "MODERATE",
        "atrt_rationale": (
            "CDK4/6 inhibitor same class as abemaciclib. Palbociclib has poor BBB "
            "vs abemaciclib (P-gp substrate). CDK4 is upregulated in ATRT (z-score 1.58). "
            "Chi 2019 AACR pediatric brain tumour data. Less preferred than abemaciclib "
            "for ATRT due to inferior CNS penetration."
        ),
        "smarcb1_relevant":  True,
        "cost_category":     "BRAND",
        "ic50_um":           0.95,
        "ic50_cell_line":    "BT16",
        "source":            "Chi 2019 AACR CT031",
    },

    {
        "name":              "Ribociclib",
        "generic_name":      "ribociclib",
        "brand_names":       ["Kisqali"],
        "targets":           ["CDK4", "CDK6"],
        "mechanism":         "cdk4/6 inhibitor",
        "drug_class":        "CDK4/6 inhibitor",
        "fda_approved":      True,
        "approved_indication": "HR+ HER2- breast cancer",
        "is_generic_available": False,
        "bbb_penetrance":    "MODERATE",
        "evidence_level":    "LOW",
        "atrt_rationale":    "CDK4/6 inhibitor class; limited CNS tumour data.",
        "smarcb1_relevant":  True,
        "cost_category":     "BRAND",
        "source":            "Class rationale",
    },

    # ════════════════════════════════════════════════════════════════════════
    # PROTEASOME INHIBITORS
    # ════════════════════════════════════════════════════════════════════════

    {
        "name":              "Bortezomib",
        "generic_name":      "bortezomib",
        "brand_names":       ["Velcade"],
        "targets":           ["PSMB5", "PSMB1"],
        "mechanism":         "proteasome inhibitor (26s, boronic acid)",
        "drug_class":        "Proteasome inhibitor",
        "fda_approved":      True,
        "approved_indication": "Multiple myeloma, mantle cell lymphoma",
        "is_generic_available": True,
        "bbb_penetrance":    "LOW",
        "evidence_level":    "LOW",
        "atrt_rationale": (
            "Proteasome inhibitor (same class as marizomib). Bortezomib has poor BBB "
            "penetrance — inferior to marizomib for ATRT. PSMB5 Chronos -3.28 in "
            "rhabdoid lines makes proteasome the most essential DepMap target. "
            "Bortezomib included as generic alternative but marizomib preferred for CNS."
        ),
        "smarcb1_relevant":  False,   # Not SMARCB1-specific; universal
        "cost_category":     "GENERIC",
        "source":            "Lin 2019 Sci Transl Med (PSMB5 context)",
    },

    # ════════════════════════════════════════════════════════════════════════
    # REPURPOSABLE GENERICS WITH EPIGENETIC ACTIVITY
    # ════════════════════════════════════════════════════════════════════════

    {
        "name":              "All-trans retinoic acid (ATRA)",
        "generic_name":      "tretinoin",
        "brand_names":       ["Vesanoid", "ATRA"],
        "targets":           ["RARA", "RARB", "RARG"],
        "mechanism":         "retinoid receptor agonist / differentiation agent",
        "drug_class":        "Retinoid (differentiation therapy)",
        "fda_approved":      True,
        "approved_indication": "Acute promyelocytic leukaemia",
        "is_generic_available": True,
        "bbb_penetrance":    "HIGH",
        "evidence_level":    "MODERATE",
        "atrt_rationale": (
            "ATRA induces differentiation in neural crest-derived tumours. "
            "ATRT-TYR subgroup has active melanocytic/neural crest program (TYR, DCT, MITF). "
            "ATRA may exploit this differentiation vulnerability. "
            "Combination with HDAC inhibitors (Hennika 2017 JQ1+ATRA synergy concept). "
            "Cheap, widely available, known paediatric safety."
        ),
        "smarcb1_relevant":  True,
        "cost_category":     "GENERIC",
        "source":            "Hennika 2017 PLOS ONE; Frühwald 2020 — ATRA in rhabdoid",
    },

    {
        "name":              "Metformin",
        "generic_name":      "metformin",
        "brand_names":       ["Glucophage"],
        "targets":           ["PRKAB1", "PRKAB2"],  # AMPK β subunit (indirect)
        "mechanism":         "ampk activator / biguanide",
        "drug_class":        "AMPK activator (repurposed antidiabetic)",
        "fda_approved":      True,
        "approved_indication": "Type 2 diabetes",
        "is_generic_available": True,
        "bbb_penetrance":    "MODERATE",
        "evidence_level":    "LOW",
        "atrt_rationale": (
            "AMPK activation suppresses mTOR signalling and MYC-driven metabolism. "
            "Metformin is well-tolerated in children. Several paediatric oncology "
            "trials ongoing. Mechanistic rationale in MYC-subgroup ATRT. "
            "Not a primary agent but may have synergistic utility."
        ),
        "smarcb1_relevant":  False,
        "cost_category":     "GENERIC",
        "source":            "Mechanistic rationale; paediatric oncology trials",
    },

    {
        "name":              "Chloroquine",
        "generic_name":      "chloroquine",
        "brand_names":       ["Aralen"],
        "targets":           ["ATP6V0A1"],  # vATPase — lysosomal autophagy
        "mechanism":         "autophagy inhibitor / lysosomotropic",
        "drug_class":        "Autophagy inhibitor (repurposed antimalarial)",
        "fda_approved":      True,
        "approved_indication": "Malaria, lupus, rheumatoid arthritis",
        "is_generic_available": True,
        "bbb_penetrance":    "HIGH",
        "evidence_level":    "LOW",
        "atrt_rationale": (
            "Autophagy inhibition potentiates proteasome inhibitor cytotoxicity "
            "(proteasome + autophagy are compensatory degradation pathways). "
            "Relevant when marizomib is in the combination: proteasome inhibition "
            "upregulates autophagy, which chloroquine blocks. Cheap generic, "
            "known paediatric safety. Lin 2019 Sci Transl Med — autophagy context."
        ),
        "smarcb1_relevant":  False,
        "cost_category":     "GENERIC",
        "source":            "Lin 2019 Sci Transl Med; Dikic 2017 Science",
    },

    {
        "name":              "Hydroxychloroquine",
        "generic_name":      "hydroxychloroquine",
        "brand_names":       ["Plaquenil"],
        "targets":           ["ATP6V0A1"],
        "mechanism":         "autophagy inhibitor / lysosomotropic",
        "drug_class":        "Autophagy inhibitor (repurposed antimalarial)",
        "fda_approved":      True,
        "approved_indication": "Malaria, lupus, rheumatoid arthritis",
        "is_generic_available": True,
        "bbb_penetrance":    "MODERATE",
        "evidence_level":    "LOW",
        "atrt_rationale":    "Less toxic alternative to chloroquine. Same autophagy mechanism.",
        "smarcb1_relevant":  False,
        "cost_category":     "GENERIC",
        "source":            "Mechanistic rationale",
    },

    {
        "name":              "Itraconazole",
        "generic_name":      "itraconazole",
        "brand_names":       ["Sporanox", "Onmel"],
        "targets":           ["PTCH1", "SMO"],   # Hedgehog pathway
        "mechanism":         "smoothened inhibitor (off-target) / antifungal",
        "drug_class":        "SMO inhibitor (repurposed antifungal)",
        "fda_approved":      True,
        "approved_indication": "Fungal infections",
        "is_generic_available": True,
        "bbb_penetrance":    "MODERATE",
        "evidence_level":    "MODERATE",
        "atrt_rationale": (
            "Itraconazole is an off-label Hedgehog pathway inhibitor (SMO). "
            "Relevant for ATRT-SHH subgroup (~37% of ATRT). "
            "Cheaper and already generically available vs vismodegib. "
            "NCT02338687 (paediatric oncology); Rudin 2014 Cancer Cell. "
            "CAUTION: significant drug-drug interactions via CYP3A4."
        ),
        "smarcb1_relevant":  True,   # SHH subgroup-specific
        "cost_category":     "GENERIC",
        "source":            "Rudin 2014 Cancer Cell; NCT02338687",
    },

    {
        "name":              "Arsenic trioxide",
        "generic_name":      "arsenic trioxide",
        "brand_names":       ["Trisenox"],
        "targets":           ["GLI1", "GLI2"],   # Hedgehog / Gli
        "mechanism":         "gli inhibitor / pro-differentiation",
        "drug_class":        "Gli inhibitor (repurposed from APL therapy)",
        "fda_approved":      True,
        "approved_indication": "Acute promyelocytic leukaemia",
        "is_generic_available": True,
        "bbb_penetrance":    "HIGH",
        "evidence_level":    "MODERATE",
        "atrt_rationale": (
            "Arsenic trioxide inhibits GLI1/GLI2 post-translationally, independent of SMO. "
            "Relevant for ATRT-SHH subgroup and potentially others. "
            "Well-characterised CNS penetration. Paediatric experience in leukaemia. "
            "Beauchamp 2011 Oncogene. Combination with HDAC-i in paediatric brain tumours."
        ),
        "smarcb1_relevant":  True,
        "cost_category":     "GENERIC",
        "source":            "Beauchamp 2011 Oncogene; GLI pathway",
    },

    {
        "name":              "Rapamycin (Sirolimus)",
        "generic_name":      "sirolimus",
        "brand_names":       ["Rapamune"],
        "targets":           ["MTOR", "RPTOR"],
        "mechanism":         "mtor complex 1 inhibitor (allosteric rapalog)",
        "drug_class":        "mTOR inhibitor (repurposed immunosuppressant)",
        "fda_approved":      True,
        "approved_indication": "Renal transplant, lymphangioleiomyomatosis",
        "is_generic_available": True,
        "bbb_penetrance":    "MODERATE",
        "evidence_level":    "LOW",
        "atrt_rationale": (
            "mTOR pathway is upregulated in ATRT (MTOR z-score 1.38 in GSE70678). "
            "Sirolimus is a generic rapalog with known paediatric PK. "
            "CAUTION: rapalogs failed as monotherapy in GBM (Wick 2011). "
            "More rational in combination with CDK4/6-i (CDK → mTOR bypass). "
            "Considered for combination schedules only."
        ),
        "smarcb1_relevant":  False,
        "cost_category":     "GENERIC",
        "source":            "mTOR pathway rationale; paediatric PK data",
    },

    {
        "name":              "Thioridazine",
        "generic_name":      "thioridazine",
        "brand_names":       ["Mellaril"],
        "targets":           ["DRD2", "DRD3", "DRD4"],
        "mechanism":         "dopamine d2/d3/d4 receptor antagonist",
        "drug_class":        "DRD2 antagonist (generic antipsychotic)",
        "fda_approved":      True,
        "approved_indication": "Schizophrenia",
        "is_generic_available": True,
        "bbb_penetrance":    "HIGH",
        "evidence_level":    "LOW",
        "atrt_rationale": (
            "Generic alternative to ONC201 (DRD2 antagonist / TRAIL inducer mechanism). "
            "Thioridazine has documented anti-cancer activity via DRD2 pathway. "
            "Very cheap; well-established CNS penetration. "
            "Less selective than ONC201 — more DRD2/3/4 off-target effects. "
            "Sachlos 2012 Cell; DRD2 expressed in ATRT (curated score 0.50)."
        ),
        "smarcb1_relevant":  False,
        "cost_category":     "GENERIC",
        "source":            "Sachlos 2012 Cell; DRD2 pathway",
    },
]


# ─────────────────────────────────────────────────────────────────────────────
# Accessor functions
# ─────────────────────────────────────────────────────────────────────────────

def get_all_generic_drugs() -> List[Dict]:
    """Return the full generic drug catalogue."""
    return GENERIC_DRUG_DATABASE.copy()


def get_generics_by_class(drug_class: str) -> List[Dict]:
    """Return all drugs matching a class (case-insensitive substring)."""
    cls_lower = drug_class.lower()
    return [d for d in GENERIC_DRUG_DATABASE if cls_lower in d.get("drug_class", "").lower()]


def get_generic_drugs_only() -> List[Dict]:
    """Return only drugs with a generic formulation available."""
    return [d for d in GENERIC_DRUG_DATABASE if d.get("is_generic_available")]


def get_drugs_by_evidence(level: str = "HIGH") -> List[Dict]:
    """Return drugs filtered by evidence level (HIGH / MODERATE / LOW)."""
    return [d for d in GENERIC_DRUG_DATABASE if d.get("evidence_level") == level.upper()]


def get_smarcb1_relevant() -> List[Dict]:
    """Return drugs with SMARCB1-loss specific rationale."""
    return [d for d in GENERIC_DRUG_DATABASE if d.get("smarcb1_relevant")]


def get_drug_names() -> List[str]:
    """Return all drug names (for pipeline candidate list expansion)."""
    return [d["name"] for d in GENERIC_DRUG_DATABASE]


def as_pipeline_candidates() -> List[Dict]:
    """
    Convert the generic drug database to the pipeline candidate dict format.
    Compatible with discovery_pipeline.ProductionPipeline candidate lists.
    """
    candidates = []
    for drug in GENERIC_DRUG_DATABASE:
        candidates.append({
            "name":        drug["name"],
            "drug_name":   drug["generic_name"],
            "targets":     drug["targets"],
            "mechanism":   drug["mechanism"],
            "drug_class":  drug["drug_class"],
            "fda_approved": drug.get("fda_approved", False),
            "is_generic":  drug.get("is_generic_available", False),
            "bbb_penetrance": drug.get("bbb_penetrance", "UNKNOWN"),
            "indication":  drug.get("approved_indication", ""),
            "evidence_level": drug.get("evidence_level", "LOW"),
            "cost_category":  drug.get("cost_category", "UNKNOWN"),
            "source":      drug.get("source", ""),
        })
    return candidates