"""
scorer.py  (ATRT v2.0)
=======================
Drug scoring utility — ATRT SMARCB1-loss context.

FIXES FROM v1.0
---------------
1. PATHWAY_WEIGHTS corrected for ATRT:
   REMOVED weights of 1.00 for DIPG-specific pathways:
     - "ACVR1 signaling": 1.00  ← DIPG-specific (ACVR1 gain-of-function mutations)
     - "BMP signaling pathway": 1.00  ← DIPG/ACVR1 mutant specific
     - "BMP-SMAD signaling": 0.95  ← DIPG-specific
     - "ALK2 signaling": 0.95  ← DIPG-specific (ALK2 = ACVR1)

   ADDED/BOOSTED ATRT-specific weights:
     - "SWI/SNF complex": 1.00  ← primary ATRT driver (SMARCB1/SMARCA4 loss)
     - "SMARCB1 chromatin remodelling": 1.00  ← defines the disease
     - "EZH2 histone methyltransferase": 1.00  ← synthetic lethality target
     - "Aurora kinase signaling": 0.88  ← AURKA-MYCN axis (Sredni 2017)
     - "SHH signaling": 0.80  ← ATRT-SHH subgroup (~37%)

2. get_pathway_weight() and DrugScorer unchanged — logic is correct.
   Only the PATHWAY_WEIGHTS dictionary values changed.

RATIONALE
----------
ACVR1 (ALK2) gain-of-function mutations are found in ~20-25% of DIPG
(Wu 2014 Nat Genet; Buczkowicz 2014 Nat Genet), NOT in ATRT. Using
ACVR1/BMP pathway weights of 1.00 in ATRT scoring is incorrect and
would incorrectly boost drugs targeting these DIPG-specific pathways.

In ATRT, the defining event is SMARCB1 biallelic loss (~95%) causing:
  1. PRC2/EZH2 hyperactivity → EZH2 inhibition is synthetic lethal
  2. BET/super-enhancer dependency → BRD4 inhibition rational
  3. AURKA/MYCN axis upregulation → AURKA inhibition rational
  4. SHH pathway activation in ATRT-SHH subgroup (~37%)

References:
  Knutson 2013 PNAS PMID 23620515 — EZH2 synthetic lethality
  Torchia 2015 Cancer Cell PMID 26609405 — HDAC/BET vulnerability
  Johann 2016 Cancer Cell PMID 26923874 — three subgroup vulnerability map
  Sredni 2017 Pediatric Blood Cancer PMID 28544500 — AURKA-MYCN
  Wilson & Roberts 2011 Nat Rev Cancer PMID 21654818 — SWI/SNF biology
  Frühwald 2020 CNS Oncology PMID 32432484 — ATRT biology review
"""

import logging
import math
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# ATRT PATHWAY WEIGHTS
#
# Scale: 1.00 = defining/essential pathway in ATRT
#        0.80-0.99 = strongly relevant
#        0.60-0.79 = moderately relevant
#        0.40-0.59 = contextually relevant
#        <0.40 = low relevance in ATRT
#
# Source for each weight group cited below.
# ─────────────────────────────────────────────────────────────────────────────

PATHWAY_WEIGHTS: Dict[str, float] = {

    # ── SWI/SNF complex — ATRT PRIMARY DRIVER ─────────────────────────────────
    # SMARCB1/SMARCA4 loss defines ATRT in ~100% of cases.
    # Source: Wilson & Roberts 2011 Nat Rev Cancer PMID 21654818;
    #         Hasselblatt 2011 Acta Neuropathol PMID 20625942
    "SWI/SNF complex":                   1.00,
    "SMARCB1 chromatin remodelling":     1.00,
    "BAF complex":                       0.98,
    "Chromatin remodeling":              0.92,

    # ── PRC2/EZH2 — synthetic lethality with SMARCB1 loss ─────────────────────
    # EZH2 is hyperactive when SMARCB1 is lost; EZH2 inhibition is the
    # primary synthetic lethality target in ATRT.
    # Source: Knutson 2013 PNAS PMID 23620515 (FDA BTD for tazemetostat)
    "EZH2 histone methyltransferase":    1.00,
    "EZH2 histone methyltransferase activity": 1.00,
    "EZH2 inhibition":                   1.00,
    "PRC2 complex":                      1.00,
    "H3K27 methylation":                 0.98,
    "Epigenetic regulation":             0.92,
    "Epigenetic regulation of gene expression": 0.90,
    "Histone modification":              0.88,

    # ── BET bromodomain — super-enhancer dependency ────────────────────────────
    # BRD4 maintains MYC super-enhancers in SMARCB1-null cells.
    # Source: Geoerger 2017 Clin Cancer Res PMID 28108534;
    #         Johann 2016 Cancer Cell PMID 26923874
    "BET bromodomain":                   0.95,
    "BRD4 signaling":                    0.95,
    "Super-enhancer regulation":         0.90,
    "Transcription factor activity":     0.78,

    # ── HDAC — epigenetic compensation target ──────────────────────────────────
    # Pan-HDAC inhibition reverses compensatory chromatin compaction in
    # SMARCB1-null cells; synergy with EZH2-i confirmed (CI 0.38–0.41 in BT16/BT37).
    # Source: Torchia 2015 Cancer Cell PMID 26609405 (Fig 4)
    "HDAC deacetylase activity":         0.95,
    "HDAC inhibition":                   0.95,
    "Histone deacetylation":             0.92,
    "Histone acetylation":               0.88,

    # ── Aurora kinase — MYCN stabilisation ────────────────────────────────────
    # AURKA phosphorylates MYCN at T58, protecting from proteasomal degradation.
    # Alisertib IC50 = 0.098 µM BT16 (Sredni 2017); synergy with BET-i CI=0.49.
    # Source: Sredni 2017 Pediatric Blood Cancer PMID 28544500;
    #         Lowery 2017 Oncotarget
    "Aurora kinase signaling":           0.88,
    "AURKA signaling":                   0.90,
    "Aurora kinase A inhibitor":         0.88,

    # ── MYC/MYCN axis ─────────────────────────────────────────────────────────
    # MYC and MYCN are transcriptional drivers across ATRT subgroups;
    # ATRT-MYC subgroup has amplification/high expression.
    # Source: Johann 2016; Torchia 2015 (GSE70678 log2FC data)
    "MYC signaling":                     0.88,
    "MYCN signaling":                    0.88,
    "Cancer stem cell signaling":        0.72,

    # ── CDK4/6 — cell cycle ───────────────────────────────────────────────────
    # CDK4 upregulated in ATRT; CDK4/6-i + tazemetostat synergy CI ~0.62 (BT16).
    # Source: Chi 2019 AACR Abstract CT031;
    #         Johann 2016 Cancer Cell PMID 26923874
    "CDK4/6 signaling":                  0.84,
    "Cell cycle regulation":             0.82,
    "Rb phosphorylation":                0.78,

    # ── Proteasome ────────────────────────────────────────────────────────────
    # PSMB5 Chronos = -3.28 in rhabdoid lines (most essential in screen).
    # MYC/MCL1 rapidly degraded by proteasome — proteasome inhibition is lethal.
    # Source: Lin 2019 Sci Transl Med; Dikic 2017 Annu Rev Biochem
    "Proteasome inhibition":             0.82,
    "Ubiquitin-proteasome pathway":      0.78,
    "Protein quality control":           0.68,

    # ── SHH pathway — ATRT-SHH subgroup (~37%) ────────────────────────────────
    # GLI2 amplification drives SHH pathway in ATRT-SHH subgroup.
    # Vismodegib IC50 = 2.50 µM BT37 (Torchia 2015); more potent in SHH lines.
    # Source: Johann 2016 Cancer Cell PMID 26923874;
    #         Torchia 2015 Cancer Cell PMID 26609405
    "SHH signaling":                     0.80,
    "Sonic Hedgehog pathway":            0.80,
    "Hedgehog pathway inhibitor":        0.80,
    "Smoothened inhibitor":              0.78,
    "GLI inhibitor":                     0.78,

    # ── Apoptosis/BCL-2 ───────────────────────────────────────────────────────
    # BCL2/BCL2L1 upregulated in ATRT; HDAC-i reduces BCL2 expression.
    # MCL1 proteasome-sensitive (t½ ~30-60 min).
    # Source: Torchia 2015; Dikic 2017
    "Apoptosis":                         0.80,
    "BCL-2 family signaling":            0.82,
    "Intrinsic apoptosis pathway":       0.80,
    "Synthetic lethality":               0.88,

    # ── mTOR/PI3K ─────────────────────────────────────────────────────────────
    # mTOR upregulated in ATRT; CDK4/6-i resistance often via PI3K/mTOR.
    # Source: Chi 2019 AACR Abstract; Torchia 2015
    "PI3K-Akt signaling":                0.75,
    "mTOR signaling":                    0.75,
    "PTEN signaling":                    0.72,

    # ── DNA damage / PARP ─────────────────────────────────────────────────────
    # HDAC-i induces BRCAness → PARP sensitisation (Adimoolam 2007).
    "PARP signaling":                    0.72,
    "DNA damage response":               0.70,
    "DNA repair":                        0.68,
    "ATR signaling":                     0.65,
    "Homologous recombination":          0.65,

    # ── p53 signalling ────────────────────────────────────────────────────────
    "p53 signaling":                     0.78,
    "MDM2-p53 interaction":              0.75,
    "Tumour suppressor loss":            0.75,
    "Oncogene addiction":                0.80,

    # ── Immune/TME ────────────────────────────────────────────────────────────
    # ATRT has low mutational burden; immune context is less central than GBM.
    # Source: Frühwald 2020
    "T-cell checkpoint signaling":       0.52,
    "PD-1/PD-L1 signaling":             0.50,
    "Tumour microenvironment":           0.55,
    "TGF-beta signaling":                0.55,

    # ── RTK / MAPK / angiogenesis ─────────────────────────────────────────────
    "Receptor tyrosine kinase signaling": 0.72,
    "MAPK signaling":                    0.68,
    "RAS signaling":                     0.65,
    "EGFR signaling":                    0.65,
    "PDGFRA signaling":                  0.65,
    "JAK-STAT signaling":                0.65,
    "STAT3 signaling":                   0.68,
    "NF-kB signaling":                   0.62,
    "Angiogenesis":                      0.55,
    "VEGF signaling":                    0.52,
    "Hypoxia response":                  0.52,
    "Cancer metabolism":                 0.55,

    # ── Stemness ──────────────────────────────────────────────────────────────
    # SOX2/LIN28A/SALL4 are highly expressed in ATRT (Torchia 2015).
    "Neural stem cell":                  0.62,
    "Stemness":                          0.62,

    # ── Generic/broad pathways ─────────────────────────────────────────────────
    "Histone H3 methylation":            0.78,
    "SMAD signaling":                    0.50,   # low — ATRT lacks ACVR1 mutations
    "BMP signaling pathway":             0.35,   # low — DIPG-specific biology
    "BMP-SMAD signaling":                0.35,   # low — DIPG-specific biology
    "ACVR1 signaling":                   0.20,   # NOT an ATRT pathway — DIPG only
    "ALK2 signaling":                    0.20,   # NOT an ATRT pathway — DIPG only
}

# ── DIPG-specific pathways that must NOT be high in ATRT ─────────────────────
# These remain in the dict with low weights so they don't accidentally boost
# drugs like ACVR1 inhibitors. They are not removed entirely in case a drug's
# pathway annotation mentions BMP for general cancer biology reasons.
_DIPG_SPECIFIC_LOW_WEIGHT_NOTE = """
ACVR1/BMP pathway weights are intentionally LOW in ATRT scoring (0.20-0.35).
ACVR1 gain-of-function mutations occur in ~20-25% of DIPG (Wu 2014 Nat Genet),
NOT in ATRT. SMAD/BMP signalling is a DIPG-specific biology.

If a drug is labelled as an "ACVR1 inhibitor" or "BMP signalling inhibitor",
it should NOT receive a high pathway score in ATRT — this was a bug in v1.0.

Source: Buczkowicz 2014 Nat Genet PMID 24705254 (ACVR1 in DIPG);
        Hasselblatt 2011 Acta Neuropathol PMID 20625942 (ATRT has no ACVR1 mutations)
"""


def get_pathway_weight(
    pathway_name: str,
    disease_weights: Optional[Dict[str, float]] = None,
) -> float:
    """
    Return the pathway weight for ATRT scoring.

    Checks disease_weights override first (e.g. subgroup-specific weights),
    then falls back to PATHWAY_WEIGHTS, then to 0.50 default.
    """
    if disease_weights:
        if pathway_name in disease_weights:
            return disease_weights[pathway_name]
        for key, weight in disease_weights.items():
            if key.lower() in pathway_name.lower() or pathway_name.lower() in key.lower():
                return weight

    if pathway_name in PATHWAY_WEIGHTS:
        return PATHWAY_WEIGHTS[pathway_name]

    pathway_lower = pathway_name.lower()
    best_weight  = None
    best_match_len = 0
    for key, weight in PATHWAY_WEIGHTS.items():
        key_lower = key.lower()
        if key_lower in pathway_lower or pathway_lower in key_lower:
            match_len = len(min(key_lower, pathway_lower, key=len))
            if match_len > best_match_len:
                best_match_len = match_len
                best_weight    = weight

    return best_weight if best_weight is not None else 0.50


# ─────────────────────────────────────────────────────────────────────────────
# Component scoring functions — inputs to hypothesis_generator
# ─────────────────────────────────────────────────────────────────────────────

def score_gene_overlap(
    drug_targets:  List[str],
    disease_genes: List[str],
    gene_weights:  Optional[Dict[str, float]] = None,
) -> Tuple[float, Dict]:
    if not drug_targets or not disease_genes:
        return 0.0, {"overlap": [], "n_overlap": 0, "n_disease": len(disease_genes)}

    drug_set    = set(t.upper() for t in drug_targets)
    disease_set = set(g.upper() for g in disease_genes)
    overlap     = drug_set & disease_set

    if not overlap:
        return 0.0, {"overlap": [], "n_overlap": 0, "n_disease": len(disease_set)}

    if gene_weights:
        weighted_hits  = sum(gene_weights.get(g, 0.5) for g in overlap)
        weighted_total = sum(gene_weights.get(g, 0.5) for g in disease_set)
        score = weighted_hits / max(weighted_total, 1.0)
    else:
        score = len(overlap) / len(disease_set)

    # Multi-hit bonus: hitting several ATRT vulnerabilities simultaneously
    if len(overlap) >= 3:
        score = min(1.0, score * 1.15)
    elif len(overlap) >= 2:
        score = min(1.0, score * 1.08)

    return round(min(score, 1.0), 4), {
        "overlap": sorted(overlap),
        "n_overlap": len(overlap),
        "n_disease": len(disease_set),
        "n_drug_targets": len(drug_set),
    }


def score_pathway_overlap(
    drug_pathways:    List[str],
    disease_pathways: List[str],
    disease_weights:  Optional[Dict[str, float]] = None,
) -> Tuple[float, Dict]:
    if not drug_pathways or not disease_pathways:
        return 0.0, {"matched_pathways": [], "weighted_score": 0.0}

    matched   = []
    total_w   = 0.0
    max_total = sum(
        get_pathway_weight(dp, disease_weights) for dp in disease_pathways
    )

    for drug_path in drug_pathways:
        best_w, best_match = 0.0, None
        for disease_path in disease_pathways:
            dp_lower  = drug_path.lower()
            dis_lower = disease_path.lower()
            if (dp_lower == dis_lower
                    or dp_lower in dis_lower
                    or dis_lower in dp_lower):
                w = get_pathway_weight(disease_path, disease_weights)
                if w > best_w:
                    best_w, best_match = w, disease_path
        if best_match:
            matched.append({
                "drug_pathway": drug_path,
                "matched_to":   best_match,
                "weight":       best_w,
            })
            total_w += best_w

    score = total_w / max(max_total, 1.0) if max_total > 0 else 0.0
    return round(min(score, 1.0), 4), {
        "matched_pathways": matched,
        "n_matched":        len(matched),
        "weighted_score":   round(total_w, 4),
        "max_possible":     round(max_total, 4),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Weight constants — exported for other modules
# ─────────────────────────────────────────────────────────────────────────────

WEIGHT_GENE       = 0.40
WEIGHT_PATHWAY    = 0.30
WEIGHT_BBB        = 0.20
WEIGHT_LITERATURE = 0.10
WEIGHT_PPI        = 0.15
WEIGHT_SIMILARITY = 0.10
WEIGHT_MECHANISM  = 0.20


# ─────────────────────────────────────────────────────────────────────────────
# DrugScorer — produces component scores for hypothesis_generator
# ─────────────────────────────────────────────────────────────────────────────

class DrugScorer:
    """
    Produces component evidence scores for drug candidates in ATRT context.

    Loads ATRT-specific gene weights and pathway weights when disease="atrt".
    DIPG-specific pathway weights (ACVR1/BMP) are intentionally low.

    Source: ATRT_GENE_SCORE_WEIGHTS from atrt_specialization.py;
            PATHWAY_WEIGHTS from this module (ATRT-corrected v2.0)
    """

    SUPPORTED_CNS_DISEASES = {
        "glioblastoma", "gbm", "dipg", "glioma",
        "medulloblastoma", "ependymoma", "atrt",
        "rhabdoid", "malignant rhabdoid tumor",
    }

    def __init__(
        self,
        disease:                str = "atrt",
        disease_genes:          Optional[List[str]] = None,
        disease_pathways:       Optional[List[str]] = None,
        custom_gene_weights:    Optional[Dict[str, float]] = None,
        custom_pathway_weights: Optional[Dict[str, float]] = None,
    ):
        self.disease          = disease.lower().strip()
        self.disease_genes    = disease_genes or []
        self.disease_pathways = disease_pathways or []
        self.gene_weights     = custom_gene_weights or {}
        self.pathway_weights  = dict(PATHWAY_WEIGHTS)

        if custom_pathway_weights:
            self.pathway_weights.update(custom_pathway_weights)

        # Load ATRT-specific weights
        if self.disease in ("atrt", "rhabdoid", "malignant rhabdoid tumor"):
            self._load_atrt_weights()

        # DIPG weights loaded only when disease is explicitly DIPG
        if self.disease in ("dipg", "diffuse intrinsic pontine glioma",
                            "h3k27m", "h3k27m glioma"):
            self._load_dipg_weights()

        logger.info(
            "DrugScorer v2.0 (ATRT): disease='%s', genes=%d",
            self.disease, len(self.disease_genes),
        )

    def _load_atrt_weights(self) -> None:
        """Load ATRT-specific gene weights from atrt_specialization module."""
        try:
            try:
                from .atrt_specialization import (
                    ATRT_CORE_GENES,
                    ATRT_GENE_SCORE_WEIGHTS,
                    ATRT_PATHWAY_WEIGHTS,
                )
            except ImportError:
                from atrt_specialization import (
                    ATRT_CORE_GENES,
                    ATRT_GENE_SCORE_WEIGHTS,
                    ATRT_PATHWAY_WEIGHTS,
                )

            if not self.disease_genes:
                self.disease_genes = ATRT_CORE_GENES
            if not self.gene_weights:
                self.gene_weights = ATRT_GENE_SCORE_WEIGHTS
            # Merge ATRT-specific pathway weights on top
            self.pathway_weights.update(ATRT_PATHWAY_WEIGHTS)
            logger.info("  → Loaded ATRT-specific weights (SMARCB1-loss context)")

        except ImportError:
            logger.warning(
                "  ⚠️  atrt_specialization module not found — "
                "using default ATRT PATHWAY_WEIGHTS from scorer.py"
            )

    def _load_dipg_weights(self) -> None:
        """Load DIPG-specific weights — ONLY for DIPG disease context."""
        try:
            try:
                from .dipg_specialization import (
                    DIPG_PATHWAY_WEIGHTS, DIPG_CORE_GENES,
                    get_dipg_gene_score_weights,
                )
            except ImportError:
                from dipg_specialization import (
                    DIPG_PATHWAY_WEIGHTS, DIPG_CORE_GENES,
                    get_dipg_gene_score_weights,
                )
            self.pathway_weights.update(DIPG_PATHWAY_WEIGHTS)
            if not self.disease_genes:
                self.disease_genes = DIPG_CORE_GENES
            if not self.gene_weights:
                self.gene_weights = get_dipg_gene_score_weights()
            logger.info("  → Loaded DIPG-specific weights (H3K27M context)")
        except ImportError:
            logger.warning("  ⚠️  dipg_specialization module not found")

    def score(self, candidate: Dict) -> Dict:
        """
        Compute component evidence scores for a single drug candidate.

        Returns candidate with gene_score, pathway_score, and preliminary_score
        added. These are inputs to hypothesis_generator — NOT final rankings.
        """
        targets   = candidate.get("targets") or candidate.get("drug_targets") or []
        pathways  = candidate.get("pathways") or []
        bbb       = candidate.get("bbb_score", 0.65)

        gene_sc, gene_details = score_gene_overlap(
            targets, self.disease_genes, self.gene_weights
        )
        path_sc, path_details = score_pathway_overlap(
            pathways, self.disease_pathways, self.pathway_weights
        )

        candidate["gene_score"]      = gene_sc
        candidate["pathway_score"]   = path_sc
        candidate["bbb_score"]       = bbb
        candidate["gene_overlap"]    = gene_details
        candidate["pathway_overlap"] = path_details

        # Lightweight preliminary score for filtering only
        preliminary = gene_sc * 0.50 + path_sc * 0.30 + bbb * 0.20
        candidate["preliminary_score"] = round(preliminary, 4)

        return candidate

    def score_batch(self, candidates: List[Dict]) -> List[Dict]:
        """Score all candidates (component scores for hypothesis_generator)."""
        for c in candidates:
            self.score(c)
        candidates.sort(
            key=lambda x: x.get("preliminary_score", 0), reverse=True
        )
        logger.info(
            "Scored %d candidates (ATRT component scores for hypothesis_generator)",
            len(candidates),
        )
        return candidates


def sensitivity_analysis(
    candidates: list, perturbation: float = 0.10
) -> dict:
    """
    Weight sensitivity analysis — confirms top-4 ranking stability.

    Source: Methods section — top-4 ranking stable under ±10% weight
    perturbations (Spearman ρ > 0.90).
    """
    base_weights = {
        "gene_score":    WEIGHT_GENE,
        "pathway_score": WEIGHT_PATHWAY,
        "bbb_score":     WEIGHT_BBB,
        "literature_score": WEIGHT_LITERATURE,
    }

    def _score(comp, w):
        return sum(comp.get(k, 0.0) * v for k, v in w.items())

    def _spearman(ranks_a, ranks_b):
        n = len(ranks_a)
        if n < 2:
            return 1.0
        d2 = sum((a - b) ** 2 for a, b in zip(ranks_a, ranks_b))
        return 1.0 - (6 * d2) / (n * (n * n - 1))

    def _rank(scores):
        indexed = sorted(enumerate(scores), key=lambda x: x[1], reverse=True)
        ranks   = [0] * len(scores)
        for r, (i, _) in enumerate(indexed, 1):
            ranks[i] = r
        return ranks

    if not candidates:
        return {
            "rank_correlation_min": 1.0,
            "stable": True,
            "paper_statement": "Sensitivity analysis skipped — no candidates.",
            "perturbation_results": [],
        }

    base_scores = [_score(c, base_weights) for c in candidates]
    base_ranks  = _rank(base_scores)
    results     = []
    min_rho     = 1.0

    for key in base_weights:
        for sign in (+1, -1):
            delta     = sign * perturbation
            w_perturb = {
                k: max(0.01, v + (delta if k == key else 0.0))
                for k, v in base_weights.items()
            }
            total     = sum(w_perturb.values())
            w_perturb = {k: v / total for k, v in w_perturb.items()}
            p_scores  = [_score(c, w_perturb) for c in candidates]
            p_ranks   = _rank(p_scores)
            rho       = _spearman(base_ranks, p_ranks)
            results.append({
                "weight_changed": key,
                "direction":      "up" if sign > 0 else "down",
                "new_weight":     round(w_perturb[key], 3),
                "rho":            round(rho, 4),
            })
            if rho < min_rho:
                min_rho = rho

    stable = min_rho > 0.90
    return {
        "rank_correlation_min": round(min_rho, 4),
        "stable":               stable,
        "paper_statement": (
            f"A ±{perturbation:.0%} perturbation of all scoring weights yielded "
            f"minimum Spearman ρ = {min_rho:.3f} "
            f"({'stable' if stable else 'unstable'})."
        ),
        "perturbation_results": results,
    }


# Backward compatibility alias
ProductionScorer = DrugScorer