import logging
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# DIPG core gene sets
# ─────────────────────────────────────────────────────────────────────────────

DIPG_CORE_GENES: List[str] = [
    "H3F3A", "HIST1H3B", "HIST1H3C",
    "EZH2", "EED", "SUZ12",
    "KDM6A", "KDM6B",
    "ACVR1", "BMPR1A", "BMPR2", "SMAD1", "SMAD5", "ID1", "ID2",
    "PDGFRA", "PDGFA", "PDGFB",
    "PIK3CA", "PIK3R1", "PTEN",
    "CDKN2A", "CDK4", "CDK6", "CCND1", "RB1",
    "TP53", "MDM2", "MDM4",
    "MYC", "MYCN",
    "KMT2A", "KMT2D", "SETD2", "ATRX", "DAXX",
    "BRD4", "BRD2", "BRD3",
    "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC6",
    "PARP1", "ATM", "ATR",
    "EGFR", "MET", "FGFR1",
    "SOX2", "OLIG2", "NKX2-1",
    "CD274", "PDCD1", "TIGIT",
    "IDH1", "IDH2",
]

DIPG_ACVR1_SUBTYPE_GENES: List[str] = [
    "ACVR1", "BMPR1A", "BMPR2",
    "SMAD1", "SMAD4", "SMAD5", "SMAD6", "SMAD7", "SMAD9",
    "ID1", "ID2", "ID3", "ID4",
    "NOGGIN", "CHORDIN", "GREM1",
]

DIPG_PATHWAY_WEIGHTS: Dict[str, float] = {
    "H3K27 methylation":                    1.00,
    "PRC2 complex":                         1.00,
    "Epigenetic regulation of gene expression": 0.95,
    "H3K27 demethylation":                  0.95,
    "Histone modification":                 0.90,
    "Chromatin remodeling":                 0.90,
    "HDAC deacetylase activity":            0.95,
    "HDAC inhibition":                      0.95,
    "Histone deacetylation":                0.95,
    "BRD4 signaling":                       0.95,
    "BET bromodomain":                      0.95,
    "BMP signaling pathway":                1.00,
    "ACVR1 signaling":                      1.00,
    "BMP-SMAD signaling":                   0.95,
    "PDGFRA signaling":                     0.95,
    "PI3K-Akt signaling":                   0.85,
    "mTOR signaling":                       0.80,
    "CDK4/6 signaling":                     0.90,
    "Cell cycle regulation":                0.80,
    "Rb-E2F signaling":                     0.85,
    "p53 signaling":                        0.80,
    "PARP signaling":                       0.85,
    "Proteasome inhibition":                0.80,
    "T-cell checkpoint signaling":          0.55,
    "VEGF signaling":                       0.60,
}

DIPG_RESISTANCE_GENES: Set[str] = {
    "H3F3A", "HIST1H3B", "PDGFRA", "EGFR", "PTEN", "PIK3CA",
    "CDKN2A", "CDK4", "MDM2", "TP53", "ATRX", "MYC", "MYCN",
    "KDM6A", "KDM6B", "BRD4", "ABCB1", "ABCG2", "ABCC1",
    "BCL2", "BCL2L1", "MCL1", "SOX2", "OLIG2",
}

DIPG_NOVEL_TARGETS: Dict[str, Dict] = {
    "ACVR1": {
        "rationale": "Gain-of-function mutations in ~25% DIPG; absent in adult GBM",
        "drugs_to_test": ["Dorsomorphin", "LDN-193189", "K02288", "INCB000928"],
        "novelty_score": 1.0,
        "clinical_status": "Preclinical only — no approved ACVR1 inhibitor",
    },
    "BRD4": {
        "rationale": "H3K27M creates hyperacetylated enhancers; BRD4 reads H3K27ac",
        "drugs_to_test": ["JQ1", "OTX015", "Molibresib", "ABBV-744"],
        "novelty_score": 0.85,
        "clinical_status": "BET inhibitors in adult GBM trials; DIPG data sparse",
    },
    "HDAC1/2": {
        "rationale": "H3K27M global hypomethylation balanced by hyperacetylation",
        "drugs_to_test": ["Panobinostat", "Vorinostat", "Entinostat"],
        "novelty_score": 0.75,
        "clinical_status": "Panobinostat in DIPG trials (PBTC-047)",
    },
    "CDK4/6": {
        "rationale": "CDK4 amplification in subset; Abemaciclib has best CNS penetrance",
        "drugs_to_test": ["Abemaciclib", "Ribociclib"],
        "novelty_score": 0.70,
        "clinical_status": "DIPG-specific trials needed",
    },
}

DIPG_DISEASE_PARAMS: Dict = {
    "baseline_orr":                   0.05,
    "baseline_pfs6":                  0.10,
    "bbb_barrier":                    0.80,
    "h3k27m_prevalence":              0.80,
    "acvr1_mutation_fraction":        0.25,
    "immune_desert_fraction":         0.85,
    "resistance_genes":               DIPG_RESISTANCE_GENES,
}


# ─────────────────────────────────────────────────────────────────────────────
# EZH2 inhibitor identification
# ─────────────────────────────────────────────────────────────────────────────

def _is_ezh2_inhibitor(drug: Dict) -> bool:
    """
    Returns True if the drug is an EZH2 inhibitor.

    FIX v5.5: EZH2 inhibitors are PENALISED in H3K27M DIPG because H3K27M
    dominant-negatively inhibits PRC2/EZH2. The enzyme is already suppressed —
    inhibiting it further has no mechanistic rationale.

    Sources: Bender et al. 2014 (Cancer Cell); Mohammad et al. 2017 (Nat Med).
    EZH2 score set to 0.22 in tissue_expression.py for this reason.
    This function ensures the DIPG specialization module is CONSISTENT
    with the tissue scorer — previously these were contradictory.
    """
    from .pipeline_config import EZH2_INHIBITOR

    name_lower = (drug.get("name") or drug.get("drug_name") or "").lower().strip()
    mech_lower = (drug.get("mechanism") or "").lower()
    targets    = [t.upper() for t in (drug.get("targets") or [])]

    # Check curated drug name list
    for suffix in (" hydrochloride", " hcl", " sodium", " mesylate", " phosphate"):
        name_lower = name_lower.replace(suffix, "")
    if name_lower.strip() in EZH2_INHIBITOR["known_inhibitors"]:
        return True

    # Check mechanism string
    if any(kw in mech_lower for kw in EZH2_INHIBITOR["mechanism_keywords"]):
        return True

    # Check if EZH2 is the ONLY target (single EZH2 inhibitor vs. pan-HDAC
    # that also hits EZH2 as an off-target — the latter is not penalised)
    if targets == ["EZH2"]:
        return True

    return False


# ─────────────────────────────────────────────────────────────────────────────
# DIPG-specific scoring class
# ─────────────────────────────────────────────────────────────────────────────

class DIPGSpecializedScorer:
    """
    Augments base pipeline scoring with DIPG/H3K27M-specific logic.

    v5.5 fixes:
    1. EZH2 inhibitors now PENALISED (not boosted) — consistent with tissue scorer.
       H3K27M dominant-negatively inhibits PRC2 → EZH2 inhibitors are non-rational.
    2. ACVR1 subgroup stratification: reports ACVR1-relevant candidates separately.
    3. Proteasome inhibitors correctly boosted (H3K27M creates proteotoxic stress).
    4. EZH2 removed from h3k27m_vulnerability_genes (was internally contradictory).
    """

    UNTESTED_IN_DIPG: Set[str] = {
        "tazemetostat", "abemaciclib", "ribociclib", "molibresib",
        "entinostat", "tucidinostat", "navitoclax", "selinexor",
        "alisertib", "olaparib", "niraparib", "veliparib",
        "alpelisib", "copanlisib", "afatinib",
    }

    DIPG_MECHANISM_KEYWORDS: Dict[str, float] = {
        "hdac inhibitor":           0.40,
        "histone deacetylase":      0.40,
        "pan-hdac":                 0.45,
        "class i hdac":             0.35,
        # EZH2 inhibitor: REMOVED — penalised instead (see _ezh2_penalty below)
        "bet inhibitor":            0.40,
        "brd4 inhibitor":           0.40,
        "bromodomain":              0.35,
        "cdk4":                     0.30,
        "cdk4/6":                   0.35,
        "cdk 4":                    0.30,
        "pdgfr":                    0.30,
        "pdgfra":                   0.35,
        "pi3k inhibitor":           0.25,
        "mtor inhibitor":           0.25,
        "bmp inhibitor":            0.45,
        "acvr1 inhibitor":          0.50,
        "parp inhibitor":           0.25,
        "dna damage":               0.20,
        "autophagy":                0.20,
        "proteasome inhibitor":     0.35,
        "proteasome":               0.30,
        "marizomib":                0.38,
    }

    def __init__(
        self,
        apply_bbb_penalty: bool  = True,
        novelty_bonus:     float = 0.08,
        h3k27m_bonus:      float = 0.12,
        acvr1_bonus:       float = 0.10,
        cns_drug_boost:    float = 0.05,
    ):
        self.apply_bbb_penalty = apply_bbb_penalty
        self.novelty_bonus     = novelty_bonus
        self.h3k27m_bonus      = h3k27m_bonus
        self.acvr1_bonus       = acvr1_bonus
        self.cns_drug_boost    = cns_drug_boost

        from .pipeline_config import EZH2_INHIBITOR
        self._ezh2_penalty = EZH2_INHIBITOR["composite_penalty"]
        self._ezh2_rationale = EZH2_INHIBITOR["rationale"]

        logger.info(
            "✅ DIPGSpecializedScorer v5.5: h3k27m_bonus=%.2f, acvr1_bonus=%.2f, "
            "novelty_bonus=%.2f, cns_boost=%.2f, ezh2_penalty=%.2f",
            h3k27m_bonus, acvr1_bonus, novelty_bonus, cns_drug_boost, self._ezh2_penalty,
        )

    def _h3k27m_vulnerability_score(self, drug: Dict) -> float:
        mechanism  = (drug.get("mechanism", "") or "").lower()
        name_lower = (drug.get("name", drug.get("drug_name", "")) or "").lower()
        targets    = [t.upper() for t in (drug.get("targets") or [])]

        score = 0.0
        for keyword, weight in self.DIPG_MECHANISM_KEYWORDS.items():
            if keyword in mechanism:
                score = max(score, weight)

        # FIX v5.5: EZH2 and EED/SUZ12 removed from vulnerability genes.
        # H3K27M inhibits PRC2 — drugs targeting EZH2 are non-rational.
        # PSMB5/PSMB2 (proteasome) added: H3K27M creates UPR/proteotoxic stress.
        h3k27m_vulnerability_genes = {
            "BRD4", "BRD2", "BRD3",         # BET bromodomain — H3K27M super-enhancers
            "HDAC1", "HDAC2", "HDAC3",       # HDAC — H3K27M hyperacetylation
            "HDAC4", "HDAC6",
            "PSMB5", "PSMB2", "PSMB8",       # Proteasome — H3K27M proteotoxic stress
            "PSMB1", "PSMA1", "PSMA3",       # Source: Lin et al. 2019 Sci Transl Med
            # NOTE: EZH2, EED, SUZ12 intentionally EXCLUDED.
            # They were here in v5.4 — contradicted the curated tissue score of 0.22.
        }
        target_hits = set(targets) & h3k27m_vulnerability_genes
        if target_hits:
            score = max(score, 0.35 + len(target_hits) * 0.05)

        if "valproic" in name_lower or "valproate" in name_lower:
            score = max(score, 0.40)

        return min(score * self.h3k27m_bonus / 0.40, self.h3k27m_bonus)

    def _acvr1_score(self, drug: Dict) -> float:
        mechanism = (drug.get("mechanism", "") or "").lower()
        targets   = [t.upper() for t in (drug.get("targets") or [])]

        acvr1_genes  = {"ACVR1", "BMPR1A", "BMPR2", "SMAD1", "SMAD5", "SMAD4", "ID1", "ID2"}
        bmp_keywords = ("bmp", "acvr1", "activin", "bone morphogenetic", "alk2", "type i bmp receptor")

        target_hits = set(targets) & acvr1_genes
        mech_hit    = any(kw in mechanism for kw in bmp_keywords)

        if target_hits or mech_hit:
            return self.acvr1_bonus * 0.25 if not target_hits else self.acvr1_bonus
        return 0.0

    def _novelty_bonus(self, drug: Dict) -> float:
        name_lower = (drug.get("name", drug.get("drug_name", "")) or "").lower()
        for suffix in (" hydrochloride", " hcl", " sodium", " mesylate"):
            name_lower = name_lower.replace(suffix, "")
        name_lower = name_lower.strip()

        if name_lower in self.UNTESTED_IN_DIPG:
            return self.novelty_bonus
        for untested in self.UNTESTED_IN_DIPG:
            if untested in name_lower or name_lower in untested:
                return self.novelty_bonus * 0.75
        return 0.0

    def _ezh2_inhibitor_penalty(self, drug: Dict) -> Tuple[bool, float, str]:
        """
        FIX v5.5: returns (is_ezh2_inhibitor, penalty_multiplier, rationale).

        EZH2 inhibitors receive a 0.50 composite score multiplier in H3K27M DIPG.
        This is CONSISTENT with the tissue_expression.py curated score of 0.22
        for EZH2 — which was already penalising EZH2 at the tissue level but
        the specialization module was previously BOOSTING these drugs, creating
        a contradiction that inflated EZH2 inhibitor rankings.
        """
        if _is_ezh2_inhibitor(drug):
            return True, self._ezh2_penalty, self._ezh2_rationale
        return False, 1.0, ""

    def score_candidate(self, candidate: Dict) -> Dict:
        base_score   = candidate.get("score", 0.0)

        # FIX v5.5: Check EZH2 inhibitor BEFORE computing bonuses
        is_ezh2, ezh2_mult, ezh2_rationale = self._ezh2_inhibitor_penalty(candidate)

        if is_ezh2:
            # Apply penalty to base score and skip h3k27m bonus
            # (EZH2 inhibitors don't get H3K27M vulnerability bonus)
            h3k27m_bonus = 0.0
            acvr1_bonus  = 0.0
            novelty      = 0.0
            bbb_pen      = candidate.get("bbb_penetrance", "UNKNOWN")
            cns_boost    = self.cns_drug_boost if bbb_pen == "HIGH" else 0.0
            adjusted     = round(min(1.0, base_score * ezh2_mult + cns_boost), 4)
            logger.info(
                "EZH2 inhibitor penalty applied to %s: %.3f → %.3f (×%.2f). %s",
                candidate.get("name", "?"), base_score, adjusted,
                ezh2_mult, ezh2_rationale[:80],
            )
        else:
            h3k27m_bonus = self._h3k27m_vulnerability_score(candidate)
            acvr1_bonus  = self._acvr1_score(candidate)
            novelty      = self._novelty_bonus(candidate)
            bbb_pen      = candidate.get("bbb_penetrance", "UNKNOWN")
            cns_boost    = self.cns_drug_boost if bbb_pen == "HIGH" else 0.0
            total_bonus  = h3k27m_bonus + acvr1_bonus + novelty + cns_boost
            adjusted     = min(1.0, base_score + total_bonus)
            ezh2_mult    = 1.0

        candidate["dipg_score"]       = round(adjusted, 4)
        candidate["dipg_bonus_total"] = round(adjusted - base_score, 4)
        candidate["dipg_components"]  = {
            "base_score":         round(base_score, 4),
            "h3k27m_bonus":       round(h3k27m_bonus if not is_ezh2 else 0.0, 4),
            "acvr1_bonus":        round(acvr1_bonus if not is_ezh2 else 0.0, 4),
            "novelty_bonus":      round(novelty if not is_ezh2 else 0.0, 4),
            "cns_boost":          round(cns_boost, 4),
            "is_untested_dipg":   (novelty > 0) if not is_ezh2 else False,
            "h3k27m_relevant":    (h3k27m_bonus > 0) if not is_ezh2 else False,
            "acvr1_relevant":     (acvr1_bonus > 0) if not is_ezh2 else False,
            "is_ezh2_inhibitor":  is_ezh2,
            "ezh2_penalty_applied": is_ezh2,
            "ezh2_penalty_multiplier": ezh2_mult,
            "ezh2_rationale":     ezh2_rationale if is_ezh2 else "",
        }
        candidate["score"] = adjusted
        return candidate

    def score_batch(self, candidates: List[Dict]) -> List[Dict]:
        logger.info("🧬 DIPG specialization v5.5: scoring %d candidates", len(candidates))
        for c in candidates:
            self.score_candidate(c)
        candidates.sort(key=lambda x: x.get("score", 0), reverse=True)

        n_h3k27m  = sum(1 for c in candidates if c.get("dipg_components", {}).get("h3k27m_relevant"))
        n_acvr1   = sum(1 for c in candidates if c.get("dipg_components", {}).get("acvr1_relevant"))
        n_novel   = sum(1 for c in candidates if c.get("dipg_components", {}).get("is_untested_dipg"))
        n_ezh2pen = sum(1 for c in candidates if c.get("dipg_components", {}).get("is_ezh2_inhibitor"))

        logger.info(
            "   H3K27M-relevant: %d | ACVR1-relevant: %d | "
            "Untested in DIPG: %d | EZH2 inhibitors penalised: %d",
            n_h3k27m, n_acvr1, n_novel, n_ezh2pen,
        )
        return candidates

    def generate_acvr1_subgroup_report(self, candidates: List[Dict]) -> str:
        """
        FIX v5.5: New method — report candidates most relevant to the ACVR1-mutant
        subgroup (~25% of DIPG). Helps stratify treatment planning.
        """
        acvr1_relevant = [
            c for c in candidates
            if c.get("dipg_components", {}).get("acvr1_relevant")
            and c.get("score", 0) > 0.40
        ]
        acvr1_relevant.sort(key=lambda x: x.get("score", 0), reverse=True)

        lines = [
            "# ACVR1-Mutant DIPG Subgroup Candidates (~25% of DIPG)",
            "",
            "These candidates have mechanistic rationale specifically for the ACVR1-mutant",
            "subgroup. ACVR1 mutations are near-absent in adult GBM, making these",
            "DIPG-specific vulnerabilities.",
            "",
            "| Rank | Drug | Score | BBB | ACVR1 Target | Mechanism |",
            "|------|------|-------|-----|--------------|-----------|",
        ]
        for i, c in enumerate(acvr1_relevant[:10], 1):
            targets = [t for t in (c.get("targets") or [])
                       if t.upper() in {"ACVR1", "BMPR1A", "BMPR2", "SMAD1", "SMAD5"}]
            bbb  = c.get("bbb_penetrance", "?")
            mech = (c.get("mechanism") or "")[:40]
            lines.append(
                f"| {i} | **{c.get('name', '?')}** | "
                f"{c.get('score', 0):.3f} | {bbb} | "
                f"{', '.join(targets) or 'indirect'} | {mech} |"
            )

        if not acvr1_relevant:
            lines.append("*No ACVR1-relevant candidates found above score threshold.*")

        lines += [
            "",
            f"**Estimated ACVR1-mutant patients in cohort**: ~{round(95 * 0.25):.0f} "
            f"of 95 H3K27M+ samples (~25% prevalence)",
            "",
            "**Clinical note**: ACVR1 inhibitors (dorsomorphin, LDN-193189) lack",
            "approved drug status. BMP pathway inhibition in this subgroup is an",
            "unmet clinical need.",
        ]
        return "\n".join(lines)

    def generate_novelty_report(self, candidates: List[Dict], top_n: int = 10) -> str:
        novel_candidates = [
            c for c in candidates
            if c.get("dipg_components", {}).get("is_untested_dipg")
        ]
        novel_candidates.sort(key=lambda x: x.get("score", 0), reverse=True)

        lines = [
            "# Novel DIPG Repurposing Candidates (Untested in DIPG Clinical Trials)",
            "",
            f"| Rank | Drug | Score | H3K27M | ACVR1 | EZH2-pen | CNS | Mechanism |",
            f"|------|------|-------|--------|-------|----------|-----|-----------|",
        ]
        for i, c in enumerate(novel_candidates[:top_n], 1):
            comp  = c.get("dipg_components", {})
            bbb   = c.get("bbb_penetrance", "?")
            mech  = (c.get("mechanism", "") or "")[:40]
            h3    = "✓" if comp.get("h3k27m_relevant") else "-"
            ac    = "✓" if comp.get("acvr1_relevant") else "-"
            ez    = "⚠" if comp.get("is_ezh2_inhibitor") else "-"
            lines.append(
                f"| {i} | **{c.get('drug_name', c.get('name', '?'))}** | "
                f"{c.get('score', 0):.3f} | {h3} | {ac} | {ez} | {bbb} | {mech} |"
            )

        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Helper functions (unchanged API)
# ─────────────────────────────────────────────────────────────────────────────

def get_dipg_gene_score_weights() -> Dict[str, float]:
    return {
        "H3F3A": 1.00, "HIST1H3B": 1.00, "ACVR1": 0.95, "PDGFRA": 0.90,
        "BRD4": 0.90, "KDM6A": 0.85, "KDM6B": 0.80,
        "HDAC1": 0.80, "HDAC2": 0.80, "ATRX": 0.75, "SETD2": 0.75,
        "CDK4": 0.85, "CDK6": 0.75, "CDKN2A": 0.80,
        "PIK3CA": 0.75, "PTEN": 0.80, "TP53": 0.75,
        "MYC": 0.70, "MYCN": 0.70, "EGFR": 0.65,
        "ABCB1": 0.40, "ABCG2": 0.40,
        # FIX v5.5: EZH2 weight reduced — it's already suppressed by H3K27M
        "EZH2": 0.20,
    }


def get_dipg_disease_data_supplement() -> Dict:
    gene_weights = get_dipg_gene_score_weights()
    return {
        "name":           "Diffuse intrinsic pontine glioma",
        "aliases":        ["DIPG", "H3K27M-mutant glioma", "diffuse midline glioma H3K27M"],
        "genes":          DIPG_CORE_GENES,
        "gene_scores":    gene_weights,
        "pathways":       list(DIPG_PATHWAY_WEIGHTS.keys()),
        "is_rare":        True,
        "is_pediatric":   True,
        "h3k27m_mutant":  True,
        "bbb_relevant":   True,
        "disease_params": DIPG_DISEASE_PARAMS,
        "source":         "DIPG specialization module v5.5",
        "notes": (
            "DIPG is a pediatric brainstem glioma driven by H3K27M histone mutation. "
            "~80% carry H3K27M. ACVR1 mutations in ~25%. Median OS 9-11 months. "
            "EZH2 inhibitors are non-rational in H3K27M DIPG (Bender et al. 2014)."
        ),
    }


def augment_disease_data_for_dipg(disease_data: Dict) -> Dict:
    supplement = get_dipg_disease_data_supplement()
    augmented  = dict(disease_data)

    existing_genes = set(disease_data.get("genes", []))
    dipg_genes     = set(supplement["genes"])
    augmented["genes"]         = list(dipg_genes | existing_genes)
    augmented["gene_scores"]   = supplement["gene_scores"]
    augmented["pathways"]      = supplement["pathways"]
    augmented["disease_params"] = supplement["disease_params"]
    augmented["is_pediatric"]  = True
    augmented["h3k27m_mutant"] = True
    augmented["bbb_relevant"]  = True
    augmented["source"]        = supplement["source"]
    augmented["dipg_augmented"] = True

    logger.info(
        "Disease data augmented for DIPG v5.5: %d total genes",
        len(augmented["genes"]),
    )
    return augmented