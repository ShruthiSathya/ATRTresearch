"""
run_atrt_pipeline.py
====================
Full ATRT pipeline orchestrator. Analogous to run_dipg_pipeline.py.

SETUP CHECKLIST (run before first execution)
---------------------------------------------
1. Download GSE70678 expression matrix to data/raw_omics/GSE70678_ATRT_expression.txt
   GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70678

2. Download GSE106982 to data/raw_omics/GSE106982_ATRT_methylation.txt (optional)
   GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106982

3. DepMap: already downloaded. No new files needed.
   Update depmap_essentiality.py to include ATRT_SUBTYPE_TERMS (see atrt_pipeline_config.py).

4. (Optional) CBTN ATRT genomics from Cavatica → data/validation/cbtn_genomics/atrt/

5. Run CMap query with ATRT signature:
   python -m backend.pipeline.prepare_cmap_query --disease atrt --output data/cmap_query/
   Then submit to clue.io and run integrate_cmap_results.py with --output atrt_cmap_scores.json

Run:
   python -m backend.pipeline.run_atrt_pipeline --disease atrt --top_n 20
"""

import asyncio
import json
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def _is_atrt(disease_name: str) -> bool:
    keywords = (
        "atrt", "atypical teratoid", "rhabdoid tumor", "rhabdoid tumour",
        "smarcb1", "smarcb1-deficient",
    )
    return any(k in disease_name.lower() for k in keywords)


class ATRTPipelineValidator:
    """
    Validates ATRT-specific genomic data from CBTN or GSE70678.
    Analogous to PedcBioPortalValidator in discovery_pipeline.py.

    Outputs:
    - smarcb1_loss_count: samples with SMARCB1 homozygous deletion or mutation
    - subgroup_counts: {TYR: n, SHH: n, MYC: n}
    - upregulated_genes: set of RNA-upregulated genes vs normal
    """

    def __init__(self, data_dir: str = "data/validation/cbtn_genomics/atrt/"):
        self.data_dir = Path(data_dir)
        # Lazy import to avoid circular dependency
        from .atrt_pipeline_config import ATRT_GENOMICS
        from .pipeline_config import GENOMICS as ATRT_GENOMICS
        self.config = ATRT_GENOMICS

    def validate_atrt_cohort(self) -> Dict:
        """
        Load ATRT genomic data and compute co-occurrence / subgroup statistics.
        Falls back gracefully if files not found.
        """
        import pandas as pd

        result = {
            "smarcb1_loss_count": 0,
            "smarca4_loss_count": 0,
            "subgroup_counts":    {"TYR": 0, "SHH": 0, "MYC": 0},
            "upregulated_genes":  set(),
            "total_samples":      0,
            "has_rna_data":       False,
        }

        # Try CBTN ATRT genomics
        cna_path = self.data_dir / "cna.txt"
        rna_path = self.data_dir / "rna_zscores.txt"
        mut_path = self.data_dir / "mutations.txt"

        if not cna_path.exists() and not rna_path.exists():
            logger.info(
                "ATRT genomic files not found at %s. "
                "Falling back to GSE70678 if available. "
                "This is non-fatal — pipeline will use prevalence priors.",
                self.data_dir,
            )
            return self._load_gse70678_fallback(result)

        # Load CNA for SMARCB1 deletion
        if cna_path.exists():
            try:
                cna = pd.read_csv(cna_path, sep="\t", index_col=0)
                cna.index = cna.index.astype(str).str.upper().str.strip()
                result["total_samples"] = len(cna.columns)

                for alias in self.config["smarcb1_aliases"]:
                    if alias in cna.index:
                        row = pd.to_numeric(cna.loc[alias], errors="coerce")
                        result["smarcb1_loss_count"] = int(
                            (row <= self.config["smarcb1_del_threshold"]).sum()
                        )
                        logger.info(
                            "SMARCB1 loss: %d/%d samples (CNA ≤ %d)",
                            result["smarcb1_loss_count"],
                            result["total_samples"],
                            self.config["smarcb1_del_threshold"],
                        )
                        break

                for alias in self.config["smarca4_aliases"]:
                    if alias in cna.index:
                        row = pd.to_numeric(cna.loc[alias], errors="coerce")
                        result["smarca4_loss_count"] = int(
                            (row <= self.config["smarcb1_del_threshold"]).sum()
                        )
                        break

            except Exception as e:
                logger.warning("CNA loading failed: %s", e)

        # Load RNA for upregulated genes
        if rna_path.exists():
            try:
                rna = pd.read_csv(rna_path, sep="\t", index_col=0)
                threshold = self.config["rna_upregulation_zscore"]
                atrt_indicators = self.config["rna_atrt_col_indicators"]
                normal_indicators = self.config["rna_normal_col_indicators"]
                metadata_rows = self.config["rna_metadata_rows"]

                gene_rows = [r for r in rna.index
                             if str(r).upper().strip() not in metadata_rows]

                if len(gene_rows) >= self.config["rna_min_genes_required"]:
                    rna_expr = rna.loc[gene_rows]
                    atrt_cols = [c for c in rna_expr.columns
                                 if any(ind.lower() in c.lower() for ind in atrt_indicators)
                                 and not any(n.lower() in c.lower() for n in normal_indicators)]

                    if atrt_cols:
                        mean_expr = rna_expr[atrt_cols].apply(
                            pd.to_numeric, errors="coerce"
                        ).mean(axis=1)
                        result["upregulated_genes"] = set(
                            mean_expr[mean_expr > threshold].index
                        )
                        result["has_rna_data"] = (
                            len(result["upregulated_genes"]) >= 10
                        )
                        logger.info(
                            "ATRT RNA: %d upregulated genes (z > %.1f) in %d samples",
                            len(result["upregulated_genes"]), threshold, len(atrt_cols),
                        )

            except Exception as e:
                logger.warning("RNA loading failed: %s", e)

        return result

    def _load_gse70678_fallback(self, result: Dict) -> Dict:
        """Load GSE70678 as RNA fallback when CBTN data unavailable."""
        import pandas as pd

        gse_path = Path("data/raw_omics/GSE70678_ATRT_expression.txt")
        if not gse_path.exists():
            logger.info(
                "GSE70678 not found. Using prevalence priors only. "
                "Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70678"
            )
            # Use published prevalence priors (Johann 2016)
            result["smarcb1_loss_count"] = 0
            result["total_samples"] = 0
            result["subgroup_counts"] = {"TYR": 18, "SHH": 18, "MYC": 13}  # Johann 2016 n=49
            return result

        try:
            df = pd.read_csv(gse_path, sep="\t", index_col=0)
            df.index = df.index.astype(str).str.upper().str.strip()

            # Separate ATRT from normal columns using GSE70678 naming
            atrt_indicators = self.config["rna_atrt_col_indicators"]
            normal_indicators = self.config["rna_normal_col_indicators"]

            atrt_cols = [c for c in df.columns
                         if any(ind.lower() in c.lower() for ind in atrt_indicators)
                         and not any(n.lower() in c.lower() for n in normal_indicators)]
            normal_cols = [c for c in df.columns
                           if any(n.lower() in c.lower() for n in normal_indicators)]

            if not atrt_cols:
                logger.warning("No ATRT columns identified in GSE70678. Check column names.")
                return result

            result["total_samples"] = len(atrt_cols)
            threshold = self.config["rna_upregulation_zscore"]

            df_numeric = df.apply(pd.to_numeric, errors="coerce")
            atrt_mean = df_numeric[atrt_cols].mean(axis=1)

            if normal_cols:
                normal_mean = df_numeric[normal_cols].mean(axis=1)
                diff = atrt_mean - normal_mean
            else:
                diff = atrt_mean

            result["upregulated_genes"] = set(diff[diff > threshold].index)
            result["has_rna_data"] = len(result["upregulated_genes"]) >= 10
            result["smarcb1_loss_count"] = len(atrt_cols)  # All ATRT are SMARCB1-null

            logger.info(
                "GSE70678 fallback: %d ATRT samples, %d upregulated genes",
                len(atrt_cols), len(result["upregulated_genes"]),
            )

        except Exception as e:
            logger.warning("GSE70678 loading failed: %s", e)

        return result


async def run_atrt_pipeline(
    disease_name: str = "atrt",
    subgroup: Optional[str] = None,
    top_n: int = 20,
    location: str = "unknown",  # "infratentorial", "supratentorial", "unknown"
    predict_combinations: bool = True,
) -> Dict:
    """
    Run the full ATRT scoring pipeline.

    Parameters
    ----------
    disease_name : str
        Should be "atrt" or similar.
    subgroup : str or None
        "TYR", "SHH", "MYC", or None (pan-ATRT scoring).
    top_n : int
        Number of top candidates to return.
    location : str
        Tumour location — affects BBB penalty severity.
    predict_combinations : bool
        Whether to run synergy prediction.
    """
    from backend.pipeline.discovery_pipeline import ProductionPipeline
    from backend.pipeline.atrt_specialization import (
        ATRTSpecializedScorer,
        ATRT_CORE_GENES,
        get_atrt_disease_data_supplement,
        augment_disease_data_for_atrt,
    )
    from backend.pipeline.polypharmacology import PolypharmacologyScorer
    from backend.pipeline.synergy_predictor import SynergyPredictor
    from backend.pipeline.atrt_pipeline_config import (
        ATRT_BBB_PENALTIES,
        ATRT_COMPOSITE_WEIGHTS,
    )

    logger.info("=" * 70)
    logger.info("ATRT Drug Repurposing Pipeline v1.0")
    logger.info("  Subgroup: %s | Location: %s", subgroup or "pan-ATRT", location)
    logger.info("=" * 70)

    # ── Step 1: Base pipeline ─────────────────────────────────────────────────
    logger.info("[1/5] Running base ProductionPipeline for ATRT...")
    pipeline = ProductionPipeline()
    await pipeline.initialize(disease=disease_name)
    results = await pipeline.run(disease_name=disease_name, top_k=top_n)

    hypotheses = results.get("hypotheses", [])
    stats = results.get("stats", {})

    # ── Step 2: ATRT genomic validation ────────────────────────────────────────
    logger.info("[2/5] Running ATRT genomic validation...")
    validator = ATRTPipelineValidator()
    genomic_stats = validator.validate_atrt_cohort()

    logger.info(
        "  SMARCB1 loss: %d samples | RNA upregulated: %d genes | has_rna: %s",
        genomic_stats.get("smarcb1_loss_count", 0),
        len(genomic_stats.get("upregulated_genes", set())),
        genomic_stats.get("has_rna_data", False),
    )

    # ── Step 3: ATRT specialization ───────────────────────────────────────────
    logger.info("[3/5] Applying ATRT SMARCB1/EZH2 specialization...")
    candidates = await pipeline._data_fetcher.fetch_approved_drugs()

    if not candidates:
        logger.warning("No candidates — using fallback library")
        candidates = [
            {"name": "Tazemetostat",  "targets": ["EZH2"]},
            {"name": "Panobinostat",  "targets": ["HDAC1", "HDAC2", "HDAC3"]},
            {"name": "Alisertib",     "targets": ["AURKA"]},
            {"name": "Birabresib",    "targets": ["BRD4", "BRD2", "BRD3"]},
            {"name": "Abemaciclib",   "targets": ["CDK4", "CDK6"]},
            {"name": "Marizomib",     "targets": ["PSMB5", "PSMB2", "PSMB1"]},
        ]

    atrt_scorer = ATRTSpecializedScorer(
        subgroup=subgroup,
        smarcb1_synleth_bonus=0.15,
        novelty_bonus=0.08,
    )
    candidates = atrt_scorer.score_batch(candidates)

    # ── Step 4: ATRT BBB penalty (location-aware) ─────────────────────────────
    bbb_penalties = ATRT_BBB_PENALTIES.get(location, ATRT_BBB_PENALTIES["unknown_location"])
    n_penalised = 0
    for c in candidates:
        bbb = c.get("bbb_penetrance", "UNKNOWN")
        if bbb in bbb_penalties:
            c["score_before_atrt_bbb"] = c["score"]
            c["score"] = round(c["score"] * bbb_penalties[bbb], 4)
            c["atrt_bbb_penalty"] = bbb_penalties[bbb]
            n_penalised += 1

    if n_penalised:
        logger.info(
            "ATRT BBB penalty applied (%s location): %d candidates penalised",
            location, n_penalised,
        )

    # ── Step 5: IC50 annotation ────────────────────────────────────────────────
    try:
        from backend.pipeline.published_ic50_atrt_validation import (
            annotate_atrt_candidates_with_ic50,
        )
        candidates = annotate_atrt_candidates_with_ic50(candidates)
    except Exception as e:
        logger.debug("ATRT IC50 annotation skipped: %s", e)

    # Sort final
    top_candidates = sorted(
        candidates, key=lambda x: x.get("score", 0), reverse=True
    )[:top_n]

    # ── Subgroup report ────────────────────────────────────────────────────────
    subgroup_report = atrt_scorer.generate_subgroup_report(top_candidates)

    # ── Combination prediction ─────────────────────────────────────────────────
    top_combinations = []
    if predict_combinations:
        syn_predictor = SynergyPredictor()
        top_combinations = syn_predictor.predict_top_combinations(top_candidates[:30])

    # ── Summary logging ────────────────────────────────────────────────────────
    logger.info("\n" + "=" * 70)
    logger.info("ATRT RESULTS — Top 10")
    logger.info("%-30s %6s %8s %8s %8s", "Drug", "Score", "EZH2↑", "AURKA", "BBB")
    logger.info("-" * 65)
    for c in top_candidates[:10]:
        name   = (c.get("name") or "?")[:29]
        score  = c.get("score", 0)
        ezh2b  = "YES" if c.get("atrt_components", {}).get("ezh2_boosted") else "-"
        aurka  = "YES" if c.get("atrt_components", {}).get("is_aurka_inhibitor") else "-"
        bbb    = c.get("bbb_penetrance", "?")
        logger.info("%-30s %6.3f %8s %8s %8s", name, score, ezh2b, aurka, bbb)

    return {
        "hypotheses":         hypotheses,
        "top_candidates":     top_candidates,
        "subgroup_report":    subgroup_report,
        "top_combinations":   top_combinations,
        "genomic_stats":      {
            k: v for k, v in genomic_stats.items()
            if k != "upregulated_genes"
        },
        "stats": {
            **stats,
            "smarcb1_loss_count":   genomic_stats.get("smarcb1_loss_count", 0),
            "rna_upregulated_genes": len(genomic_stats.get("upregulated_genes", set())),
            "subgroup":             subgroup or "pan-ATRT",
            "location":             location,
            "escape_bypass_mode":   "RNA-confirmed" if genomic_stats.get("has_rna_data") else "curated fallback",
            "n_ezh2_boosted":       sum(
                1 for c in top_candidates
                if c.get("atrt_components", {}).get("ezh2_boosted")
            ),
        },
        "pipeline_stats": {
            "disease":           "ATRT",
            "subgroup":          subgroup or "pan-ATRT",
            "total_candidates":  len(candidates),
            "location":          location,
            "v1_features": [
                "EZH2 inhibitor BOOST (synthetic lethality with SMARCB1 loss)",
                "AURKA inhibitor boost (MYCN stabilisation, MYC subgroup)",
                "SMO/GLI inhibitor boost (SHH subgroup)",
                "Location-aware BBB penalty (infratentorial vs supratentorial)",
                "GSE70678 bulk RNA tissue scoring",
                "Published IC50 validation (BT16, BT37, G401, A204, CHLA lines)",
                "Subgroup stratification report (TYR/SHH/MYC — Johann 2016)",
            ],
            "key_difference_from_dipg": (
                "EZH2 inhibitors (tazemetostat) are BOOSTED in ATRT due to SMARCB1 "
                "synthetic lethality — the OPPOSITE of DIPG where EZH2 inhibitors "
                "are penalised. This is the most important disease-specific difference."
            ),
        },
    }


async def _cli_main(
    disease: str,
    subgroup: Optional[str],
    location: str,
    output: Optional[str],
    top_n: int,
    combinations: bool,
):
    result = await run_atrt_pipeline(
        disease_name=disease,
        subgroup=subgroup,
        top_n=top_n,
        location=location,
        predict_combinations=combinations,
    )

    print("\n" + "=" * 70)
    print("ATRT TOP CANDIDATES")
    print("=" * 70)
    for i, c in enumerate(result["top_candidates"][:10], 1):
        ezh2 = "⬆ EZH2-BOOST" if c.get("atrt_components", {}).get("ezh2_boosted") else ""
        ic50 = f"  IC50={c.get('ic50_um', 'N/A')} µM ({c.get('ic50_cell_line', '')})" if c.get("ic50_validated") else ""
        print(f"  {i:2d}. {c.get('name','?'):<28} score={c.get('score',0):.3f}  "
              f"BBB={c.get('bbb_penetrance','?'):<10}{ezh2}{ic50}")

    print("\n" + "=" * 70)
    print("SUBGROUP STRATIFICATION REPORT")
    print("=" * 70)
    print(result["subgroup_report"])

    print("\n" + "=" * 70)
    print("KEY PIPELINE STATS")
    print("=" * 70)
    s = result["stats"]
    print(f"  Subgroup:            {s.get('subgroup', 'pan-ATRT')}")
    print(f"  EZH2 inhibitors boosted: {s.get('n_ezh2_boosted', 0)}")
    print(f"  RNA upregulated genes:   {s.get('rna_upregulated_genes', 0)}")
    print(f"  Escape bypass mode:      {s.get('escape_bypass_mode', 'N/A')}")

    if output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w") as f:
            json.dump(result, f, indent=2, default=str)
        print(f"\nResults saved to: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ATRT Drug Repurposing Pipeline v1.0")
    parser.add_argument("--disease",   default="atrt")
    parser.add_argument("--subgroup",  default=None,
                        choices=["TYR", "SHH", "MYC", None],
                        help="ATRT molecular subgroup (omit for pan-ATRT)")
    parser.add_argument("--location",  default="unknown",
                        choices=["infratentorial", "supratentorial", "unknown"],
                        help="Tumour location for BBB penalty")
    parser.add_argument("--output",    default="results/atrt_pipeline_results.json")
    parser.add_argument("--top_n",     default=20, type=int)
    parser.add_argument("--combinations", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )
    asyncio.run(_cli_main(
        args.disease, args.subgroup, args.location,
        args.output, args.top_n, args.combinations,
    ))