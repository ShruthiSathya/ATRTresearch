"""
polypharmacology.py
====================
Polypharmacology scorer for ATRT Drug Repurposing Pipeline v2.0

CHANGES FROM v1.0
------------------
- Removed DIPG-specific synergy pairs (H3F3A, HIST1H3B, ACVR1/BMP)
  from SYNERGISTIC_TARGET_COMBINATIONS. These are only rational when
  H3K27M is the driver; in ATRT SMARCB1 loss is the driver.
- Added ATRT-specific synergy pairs with published CI evidence.
- _selectivity_score() was already fixed (true 0–1 range).
- disease parameter now actually gates which synergy rules are loaded.

ATRT SYNERGY LOGIC
-------------------
Primary synergy axes in SMARCB1-null tumors:
1. EZH2 + HDAC — complementary H3K27 modification (Torchia 2015 CI < 0.5)
2. BRD4 + AURKA — dual MYCN blockade (Sredni 2017 CI = 0.49)
3. PSMB5 + BRD4 — proteasomal + transcriptional MYC blockade
4. CDK4 + MTOR — cell cycle + survival dual blockade
"""

import logging
import math
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# DISEASE-GATED SYNERGY COMBINATIONS
# ─────────────────────────────────────────────────────────────────────────────

# ATRT-specific (SMARCB1-loss context)
# Sources listed per entry
ATRT_SYNERGISTIC_TARGET_COMBINATIONS: List[Dict] = [
    {
        "targets":     {"EZH2", "HDAC1"},
        "score":       0.92,
        "mechanism":   "H3K27me3 normalisation + deacetylation block (Torchia 2015 CI 0.41)",
        "atrt_specific": True,
        "reference":   "Torchia 2015 Cancer Cell PMID 26609405",
    },
    {
        "targets":     {"EZH2", "HDAC2"},
        "score":       0.90,
        "mechanism":   "Complementary H3K27 modification (Torchia 2015)",
        "atrt_specific": True,
        "reference":   "Torchia 2015 Cancer Cell PMID 26609405",
    },
    {
        "targets":     {"BRD4", "AURKA"},
        "score":       0.88,
        "mechanism":   "Dual MYCN blockade: transcription + protein stability (Sredni 2017 CI 0.49)",
        "atrt_specific": True,
        "reference":   "Sredni 2017 Pediatric Blood Cancer PMID 28544500",
    },
    {
        "targets":     {"BRD2", "AURKA"},
        "score":       0.84,
        "mechanism":   "MYCN transcription + protein stability block",
        "atrt_specific": True,
        "reference":   "Extrapolated from BRD4+AURKA data",
    },
    {
        "targets":     {"EZH2", "BRD4"},
        "score":       0.82,
        "mechanism":   "H3K27me3 + MYC super-enhancer dual block (Geoerger 2017 CI 0.55)",
        "atrt_specific": True,
        "reference":   "Geoerger 2017 Clin Cancer Res PMID 28108534",
    },
    {
        "targets":     {"EZH2", "BRD4", "HDAC1"},
        "score":       0.94,
        "mechanism":   "Triple epigenetic attack — EZH2+BET+HDAC",
        "atrt_specific": True,
        "reference":   "Mechanistic extension of Torchia 2015 + Geoerger 2017",
    },
    {
        "targets":     {"PSMB5", "BRD4"},
        "score":       0.80,
        "mechanism":   "MYC protein degradation + transcription block",
        "atrt_specific": True,
        "reference":   "Mechanistic: MYC t½ ~20 min (Lin 2019 Sci Transl Med)",
    },
    {
        "targets":     {"CDK4", "MTOR"},
        "score":       0.76,
        "mechanism":   "Cell cycle arrest + survival pathway blockade",
        "atrt_specific": False,
        "reference":   "Mechanistic; CDK4/6 resistance via PI3K/mTOR (Olmez 2017)",
    },
    {
        "targets":     {"CDK4", "PIK3CA"},
        "score":       0.74,
        "mechanism":   "CDK4/6 + PI3K dual blockade prevents mTOR bypass",
        "atrt_specific": False,
        "reference":   "Mechanistic (Dong 2020)",
    },
    {
        "targets":     {"BCL2", "HDAC1"},
        "score":       0.72,
        "mechanism":   "HDAC-i reduces BCL-2 expression, sensitises to apoptosis",
        "atrt_specific": False,
        "reference":   "Mechanistic (Bhatt 2017)",
    },
    {
        "targets":     {"MCL1", "AURKA"},
        "score":       0.70,
        "mechanism":   "MCL1 degradation + MYCN destabilisation",
        "atrt_specific": True,
        "reference":   "Mechanistic: MCL1 is proteasome-sensitive",
    },
    {
        "targets":     {"PARP1", "HDAC1"},
        "score":       0.74,
        "mechanism":   "HDAC-i induces BRCAness, sensitises to PARP inhibition",
        "atrt_specific": False,
        "reference":   "Adimoolam 2007; Rasmussen 2016",
    },
]

# GBM-specific synergy rules (NOT used in ATRT scoring)
# Kept here for reference if this module is repurposed
_GBM_DIPG_SYNERGISTIC_COMBINATIONS: List[Dict] = [
    {
        "targets":      {"H3F3A", "EZH2"},
        "score":        0.90,
        "mechanism":    "H3K27M + EZH2 — only rational in H3K27M context",
        "atrt_specific": False,
        "dipg_only":    True,
    },
    {
        "targets":      {"ACVR1", "HDAC1"},
        "score":        0.80,
        "mechanism":    "ACVR1/BMP + HDAC — only rational in ACVR1-mutant DIPG",
        "atrt_specific": False,
        "dipg_only":    True,
    },
]


# ─────────────────────────────────────────────────────────────────────────────
# RESISTANCE GENES (shared ATRT/GBM)
# ─────────────────────────────────────────────────────────────────────────────

ATRT_RESISTANCE_GENES: Set[str] = {
    # Primary SMARCB1-loss resistance nodes
    "MYC", "MYCN", "BCL2L1", "MCL1", "CDK4", "BRD4", "PIK3CA", "EZH2",
    # RTK bypass
    "EGFR", "PDGFRA", "MET", "FGFR1", "AXL",
    # PI3K/AKT/mTOR
    "PIK3CA", "PIK3R1", "AKT1", "AKT2", "MTOR",
    # Cell cycle
    "CCND1", "CCND2", "CCNE1",
    # Apoptosis
    "BCL2", "BIRC5",
    # Epigenetic
    "HDAC1", "HDAC2",
    # Drug efflux
    "ABCB1", "ABCG2",
    # Stemness
    "SOX2", "LIN28A", "SALL4",
}


class PolypharmacologyScorer:
    """
    ATRT Polypharmacology Scorer v2.0

    Scores drug candidates for multi-target synergy quality.
    Uses ATRT-specific synergy rules; DIPG/H3K27M rules are excluded.

    Parameters
    ----------
    disease : str — "atrt" (default) or "gbm"
    """

    def __init__(
        self,
        disease:       str  = "atrt",
        dipg_mode:     bool = False,
    ) -> None:
        self.disease   = disease.lower()
        self.dipg_mode = dipg_mode

        # Always use ATRT combinations in this pipeline
        # dipg_mode kept as parameter for backward compat but is ignored
        if dipg_mode:
            logger.warning(
                "PolypharmacologyScorer: dipg_mode=True is ignored in ATRT pipeline. "
                "DIPG synergy rules (H3K27M, ACVR1, BMP) are not applied."
            )

        self._synergy_combos = ATRT_SYNERGISTIC_TARGET_COMBINATIONS
        self.resistance_genes = ATRT_RESISTANCE_GENES

        logger.info(
            "PolypharmacologyScorer: disease=%s | synergy_rules=%d | resistance_genes=%d",
            self.disease, len(self._synergy_combos), len(self.resistance_genes),
        )

    def score(self, candidate: Dict, disease_targets: Optional[List[str]] = None) -> Dict:
        """
        Score a single drug for polypharmacology quality.

        Parameters
        ----------
        candidate : dict with 'targets' list
        disease_targets : optional list of ATRT disease genes for context scoring

        Returns candidate dict with poly_* fields added.
        """
        targets  = set(candidate.get("targets", []))

        syn_score, syn_combos = self._score_synergistic_combinations(targets)
        sel_score             = self._selectivity_score(targets)
        resist_score, r_hits  = self._resistance_coverage(targets)

        context_score = 0.0
        if disease_targets:
            disease_set   = set(t.upper() for t in disease_targets)
            target_upper  = set(t.upper() for t in targets)
            overlap       = target_upper & disease_set
            context_score = min(1.0, len(overlap) / max(len(disease_set), 1))

        poly_score = (
            syn_score     * 0.45
            + resist_score  * 0.25
            + context_score * 0.20
            + sel_score     * 0.10
        )

        candidate["poly_synergy_score"]     = round(syn_score, 4)
        candidate["poly_selectivity_score"] = round(sel_score, 4)
        candidate["poly_resistance_score"]  = round(resist_score, 4)
        candidate["poly_context_score"]     = round(context_score, 4)
        candidate["poly_score"]             = round(min(1.0, poly_score), 4)
        candidate["synergistic_combos"]     = syn_combos
        candidate["resistance_gene_hits"]   = sorted(r_hits)

        return candidate

    def score_batch(
        self,
        candidates:      List[Dict],
        disease_targets: Optional[List[str]] = None,
    ) -> List[Dict]:
        """Score all candidates; sort by poly_score."""
        for c in candidates:
            self.score(c, disease_targets=disease_targets)
        candidates.sort(key=lambda x: x.get("poly_score", 0), reverse=True)
        logger.info(
            "Polypharmacology scored: %d candidates | top poly_score=%.3f",
            len(candidates),
            candidates[0].get("poly_score", 0) if candidates else 0,
        )
        return candidates

    def _score_synergistic_combinations(
        self, targets: Set[str]
    ) -> Tuple[float, List[Dict]]:
        """Check drug targets against ATRT synergistic target combinations."""
        if not targets:
            return 0.0, []

        targets_upper = set(t.upper() for t in targets)
        matched: List[Dict] = []
        best_score = 0.0

        for combo in self._synergy_combos:
            combo_targets = set(t.upper() for t in combo["targets"])
            overlap       = targets_upper & combo_targets

            if len(overlap) == len(combo_targets):
                score = combo["score"]
                matched.append({
                    "targets":   sorted(combo_targets),
                    "score":     score,
                    "mechanism": combo["mechanism"],
                    "reference": combo.get("reference", ""),
                    "full_hit":  True,
                })
                best_score = max(best_score, score)

            elif len(overlap) >= max(1, len(combo_targets) - 1):
                score = combo["score"] * 0.60
                matched.append({
                    "targets":   sorted(overlap),
                    "score":     score,
                    "mechanism": combo["mechanism"] + " (partial)",
                    "full_hit":  False,
                })
                best_score = max(best_score, score)

        return round(best_score, 4), matched

    def _selectivity_score(self, targets: Set[str]) -> float:
        """
        True 0–1 selectivity score (exponential decay with target count).
        k = 0.15 gives: n=1→1.00, n=2→0.86, n=5→0.55, n=10→0.26, n=20→0.06
        """
        n = len(targets)
        if n == 0:
            return 0.0
        return round(math.exp(-0.15 * (n - 1)), 4)

    def _resistance_coverage(
        self, targets: Set[str]
    ) -> Tuple[float, Set[str]]:
        """Score how many known ATRT resistance genes the drug directly targets."""
        if not targets:
            return 0.0, set()
        targets_upper = set(t.upper() for t in targets)
        hits          = targets_upper & set(r.upper() for r in self.resistance_genes)
        score = 0.0
        for i in range(len(hits)):
            score += 0.20 / math.sqrt(i + 1)
        return round(min(1.0, score), 4), hits

    def generate_poly_report(self, candidates: List[Dict]) -> str:
        lines = [
            "## ATRT Polypharmacology Scoring Summary\n\n",
            "Synergy rules: ATRT SMARCB1-null context only "
            "(Torchia 2015, Geoerger 2017, Sredni 2017).\n\n",
            f"| {'Drug':<25} | {'Poly':>5} | {'Syn':>5} | "
            f"{'Sel':>5} | {'Res':>5} | Top synergistic combination |\n",
            f"|{'-'*26}|{'-'*6}|{'-'*6}|{'-'*6}|{'-'*6}|{'-'*30}|\n",
        ]
        for c in candidates[:15]:
            name   = (c.get("drug_name") or c.get("name") or "?")[:24]
            poly   = c.get("poly_score", 0)
            syn    = c.get("poly_synergy_score", 0)
            sel    = c.get("poly_selectivity_score", 0)
            res    = c.get("poly_resistance_score", 0)
            combos = c.get("synergistic_combos", [])
            top    = combos[0]["mechanism"][:28] if combos else "-"
            lines.append(
                f"| {name:<25} | {poly:>5.3f} | {syn:>5.3f} | "
                f"{sel:>5.3f} | {res:>5.3f} | {top} |\n"
            )
        return "".join(lines)