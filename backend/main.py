"""
main.py — ATRT Drug Repurposing API
FastAPI entry point. Imports ProductionPipeline from discovery_pipeline
(not production_pipeline which does not exist).
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
import logging
import traceback

# ── correct import path ───────────────────────────────────────────────────────
from pipeline.discovery_pipeline import ProductionPipeline
from pipeline.clinical_validator import ClinicalValidator
from pipeline.drug_filter import DrugSafetyFilter

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="ATRT Drug Repurposing API",
    description=(
        "Multi-omic computational pipeline for prioritising drug combinations "
        "in SMARCB1-deficient Atypical Teratoid/Rhabdoid Tumor (ATRT). "
        "Scoring: tissue (GSE70678) × 0.40 + DepMap (BT16/BT37/G401/A204) × 0.35 "
        "+ escape bypass × 0.20 + PPI × 0.05. "
        "EZH2 inhibitors boosted ×1.40 (Knutson 2013 PNAS synthetic lethality)."
    ),
    version="1.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global pipeline — initialised once at startup
_pipeline: ProductionPipeline | None = None


@app.on_event("startup")
async def startup_event() -> None:
    global _pipeline
    logger.info("🚀 Starting ATRT Drug Repurposing API v1.0 ...")
    logger.info(
        "📊 Data sources: OpenTargets (EFO_0002915), DepMap 24Q4, "
        "GSE70678, STRING-DB, clue.io CMap"
    )
    try:
        _pipeline = ProductionPipeline()
        logger.info("✅ Pipeline ready (disease=atrt)")
    except Exception:
        # Log the full traceback so the researcher knows exactly what failed
        logger.error("❌ Pipeline initialisation failed:\n%s", traceback.format_exc())
        # Re-raise so uvicorn marks the service as unhealthy rather than
        # silently accepting requests that will always return success=False
        raise


@app.on_event("shutdown")
async def shutdown_event() -> None:
    logger.info("👋 ATRT API shutdown")


# ─────────────────────────────────────────────────────────────────────────────
# Health
# ─────────────────────────────────────────────────────────────────────────────

@app.get("/", tags=["Health"])
async def root() -> dict:
    return {
        "status": "online",
        "service": "ATRT Drug Repurposing API",
        "version": "1.0.0",
        "disease": "ATRT (SMARCB1-deficient)",
        "pipeline_ready": _pipeline is not None,
    }


# ─────────────────────────────────────────────────────────────────────────────
# /analyze  — main scoring endpoint
# ─────────────────────────────────────────────────────────────────────────────

@app.post("/analyze", tags=["Analysis"])
async def analyze_disease(request: dict) -> dict:
    """
    Score drug candidates against the ATRT SMARCB1-loss molecular context.

    Request body
    ------------
    disease_name : str   (default "atrt")
    min_score    : float (default 0.20)
    max_results  : int   (default 10)
    subgroup     : str   (optional — "TYR" | "SHH" | "MYC")
    location     : str   (optional — "infratentorial" | "supratentorial" | "unknown_location")

    Returns
    -------
    success       : bool
    candidates    : list[dict]  — scored and filtered drug candidates
    filtered_count: int
    filtered_drugs: list[dict]  — safety-filtered drugs with reason
    stats         : dict        — pipeline run statistics
    """
    if _pipeline is None:
        return {
            "success": False,
            "error": (
                "Pipeline not initialised. Check startup logs for the root cause. "
                "Common causes: DepMap CSVs missing from data/depmap/, "
                "or a Python import error in discovery_pipeline.py."
            ),
        }

    disease_name = request.get("disease_name", "atrt")
    min_score    = float(request.get("min_score", 0.20))
    max_results  = int(request.get("max_results", 10))
    subgroup     = request.get("subgroup")           # TYR / SHH / MYC / None
    location     = request.get("location", "unknown_location")

    logger.info(
        "Analysis request: disease=%s subgroup=%s location=%s min_score=%.2f",
        disease_name, subgroup or "pan-ATRT", location, min_score,
    )

    try:
        result = await _pipeline.run(
            disease_name=disease_name,
            top_k=max_results * 2,   # get extra before safety filtering
            subgroup=subgroup,
            location=location,
        )
        if not result:
            return {"success": False, "error": "Pipeline returned empty result"}

        candidates = result.get("top_candidates", [])

        # Normalise field names expected by DrugSafetyFilter
        for c in candidates:
            if "indication" not in c:
                c["indication"] = "ATRT / rhabdoid tumor (SMARCB1-deficient)"
            if "mechanism" not in c:
                c["mechanism"] = c.get("drug_class", "")
            if "drug_name" not in c:
                c["drug_name"] = c.get("name", "")

        # Safety filter
        safety_filter = DrugSafetyFilter()
        try:
            safe_candidates, filtered_out = await safety_filter.filter_candidates(
                candidates=candidates,
                disease_name=disease_name,
                remove_absolute=True,
                remove_relative=True,
            )
        except Exception as filter_err:
            logger.error("Safety filter error (returning unfiltered): %s", filter_err)
            safe_candidates = candidates
            filtered_out = []

        # Apply min_score threshold and cap at max_results
        safe_candidates = [
            c for c in safe_candidates if c.get("score", 0) >= min_score
        ][:max_results]

        logger.info(
            "Result: %d safe candidates (filtered %d)",
            len(safe_candidates), len(filtered_out),
        )

        return {
            "success": True,
            "candidates": safe_candidates,
            "filtered_count": len(filtered_out),
            "filtered_drugs": [
                {
                    "drug_name": c.get("drug_name", c.get("name", "?")),
                    "reason": c.get("contraindication", {}).get("reason", "Unknown"),
                    "severity": c.get("contraindication", {}).get("severity", "unknown"),
                }
                for c in filtered_out
            ],
            "stats": result.get("stats", {}),
        }

    except Exception:
        logger.error("Analysis error:\n%s", traceback.format_exc())
        return {"success": False, "error": traceback.format_exc()}


# ─────────────────────────────────────────────────────────────────────────────
# /validate_clinical  — per-drug clinical evidence lookup
# ─────────────────────────────────────────────────────────────────────────────

@app.post("/validate_clinical", tags=["Analysis"])
async def validate_clinical(request: dict) -> dict:
    """
    Validate a drug candidate using ClinicalTrials.gov, PubMed, and OpenFDA.

    Request body
    ------------
    drug_name    : str (required)
    disease_name : str (required)
    drug_data    : dict (optional — additional drug metadata)
    disease_data : dict (optional — additional disease metadata)
    """
    drug_name    = request.get("drug_name")
    disease_name = request.get("disease_name")
    drug_data    = request.get("drug_data", {})
    disease_data = request.get("disease_data", {})

    if not drug_name or not disease_name:
        return {"success": False, "error": "Missing drug_name or disease_name"}

    logger.info("Clinical validation: %s for %s", drug_name, disease_name)

    validator = ClinicalValidator()
    try:
        validation_result = await validator.validate_candidate(
            drug_name=drug_name,
            disease_name=disease_name,
            drug_data=drug_data,
            disease_data=disease_data,
        )
        return {"success": True, "validation": validation_result}
    except Exception:
        logger.error("Clinical validation error:\n%s", traceback.format_exc())
        return {"success": False, "error": traceback.format_exc()}
    finally:
        await validator.close()


# ─────────────────────────────────────────────────────────────────────────────
# /subgroup_report  — subgroup-stratified scoring
# ─────────────────────────────────────────────────────────────────────────────

@app.post("/subgroup_report", tags=["Analysis"])
async def subgroup_report(request: dict) -> dict:
    """
    Run pan-ATRT scoring stratified by molecular subgroup (TYR / SHH / MYC).
    Returns separate candidate lists per subgroup.

    Request body
    ------------
    max_results : int (default 8 per subgroup)
    """
    if _pipeline is None:
        return {"success": False, "error": "Pipeline not initialised"}

    max_results = int(request.get("max_results", 8))
    results_by_subgroup: dict = {}

    for subgroup in ("TYR", "SHH", "MYC"):
        try:
            result = await _pipeline.run(
                disease_name="atrt",
                top_k=max_results,
                subgroup=subgroup,
                location="unknown_location",
            )
            results_by_subgroup[subgroup] = {
                "candidates": result.get("top_candidates", [])[:max_results],
                "stats":      result.get("stats", {}),
            }
        except Exception:
            logger.error(
                "Subgroup %s failed:\n%s", subgroup, traceback.format_exc()
            )
            results_by_subgroup[subgroup] = {"error": traceback.format_exc()}

    return {"success": True, "subgroups": results_by_subgroup}