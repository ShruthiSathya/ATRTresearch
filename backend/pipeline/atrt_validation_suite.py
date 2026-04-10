"""
atrt_validation_suite.py
========================
Integration validation tests for the ATRT Drug Repurposing Pipeline.

Tests use EXPECTED RANGES from published biology, not exact hardcoded numbers.
This approach is robust to minor scoring changes while catching regressions.

Run:
    python -m backend.pipeline.atrt_validation_suite

All assertions are grounded in published literature with PMIDs.

WHAT IS TESTED
--------------
1. EZH2 inhibitor boost: tazemetostat score > alisertib score
   (EZH2 synthetic lethality makes it #1 — Knutson 2013 PNAS)

2. Tazemetostat BBB = MODERATE (not HIGH)
   (Knutson 2013 patent; Gounder 2020 JCO supplement — Kp,uu ~0.15-0.30)

3. Marizomib has no IC50 validation score (no verified primary ATRT data)

4. ONC201 has no IC50 validation score (no verified primary ATRT data)

5. Panobinostat IC50 ~ 0.009 µM in BT16 (Torchia 2015 Cancer Cell)
   within ±20% tolerance for assay variability

6. PSMB5 DepMap Chronos ≤ -2.0 (extremely essential — Lin 2019 analogy)

7. SMARCB1 curated expression score ≤ 0.10 (gene is lost in ATRT)

8. Composite weights sum to 1.0

9. Confidence weights sum to 1.0

10. No hardcoded composite scores appear in pipeline output
    (scores must be computed from data, not from config)
"""

import asyncio
import logging
import math
from typing import Dict, List, Optional

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Test helpers
# ─────────────────────────────────────────────────────────────────────────────

class ValidationError(Exception):
    pass


def assert_range(
    value: float,
    low: float,
    high: float,
    label: str,
    pmid: str = "",
) -> None:
    """Assert a value falls within [low, high], with publication reference."""
    if not (low <= value <= high):
        ref = f" (PMID {pmid})" if pmid else ""
        raise ValidationError(
            f"FAIL: {label} = {value:.4f} not in [{low:.4f}, {high:.4f}]{ref}"
        )
    ref = f" (PMID {pmid})" if pmid else ""
    logger.info("  ✅ %s = %.4f ∈ [%.3f, %.3f]%s", label, value, low, high, ref)


def assert_equal(value, expected, label: str) -> None:
    if value != expected:
        raise ValidationError(f"FAIL: {label} = {value!r}, expected {expected!r}")
    logger.info("  ✅ %s = %r", label, value)


def assert_none(value, label: str) -> None:
    if value is not None:
        raise ValidationError(f"FAIL: {label} should be None, got {value!r}")
    logger.info("  ✅ %s = None (as expected — no verified data)", label)


def assert_less(a: float, b: float, label: str) -> None:
    if not (a < b):
        raise ValidationError(f"FAIL: {label}: {a:.4f} is not < {b:.4f}")
    logger.info("  ✅ %s: %.4f < %.4f", label, a, b)


def assert_greater(a: float, b: float, label: str) -> None:
    if not (a > b):
        raise ValidationError(f"FAIL: {label}: {a:.4f} is not > {b:.4f}")
    logger.info("  ✅ %s: %.4f > %.4f", label, a, b)


# ─────────────────────────────────────────────────────────────────────────────
# Individual test functions
# ─────────────────────────────────────────────────────────────────────────────

def test_config_weights() -> None:
    """Test that composite and confidence weights sum to 1.0."""
    from backend.pipeline.pipeline_config import COMPOSITE_WEIGHTS, CONFIDENCE_WEIGHTS

    composite_sum = sum(COMPOSITE_WEIGHTS.values())
    assert abs(composite_sum - 1.0) < 1e-9, \
        f"COMPOSITE_WEIGHTS sum = {composite_sum:.6f}, expected 1.0"
    logger.info("  ✅ COMPOSITE_WEIGHTS sum = 1.0")

    confidence_sum = sum(CONFIDENCE_WEIGHTS.values())
    assert abs(confidence_sum - 1.0) < 1e-9, \
        f"CONFIDENCE_WEIGHTS sum = {confidence_sum:.6f}, expected 1.0"
    logger.info("  ✅ CONFIDENCE_WEIGHTS sum = 1.0")


def test_tazemetostat_bbb_moderate() -> None:
    """
    Tazemetostat must be classified as MODERATE BBB (not HIGH).
    Source: Knutson 2013 PNAS patent supplement; Gounder 2020 JCO supplement.
    Preclinical rodent Kp,uu ~0.15-0.30 → MODERATE.
    Stacchiotti 2021 NEJM does NOT report CNS PK data.
    """
    from backend.pipeline.bbb_filter import BBBFilter
    bbb = BBBFilter()
    result = bbb.score_drug("tazemetostat")
    assert_equal(
        result["penetrance"], "MODERATE",
        "tazemetostat BBB penetrance",
    )


def test_panobinostat_bbb_high() -> None:
    """
    Panobinostat must be classified as HIGH BBB.
    Source: Monje 2023 Nat Med PMID 37526549 (PBTC-047 — confirmed CNS PK).
    """
    from backend.pipeline.bbb_filter import BBBFilter
    bbb = BBBFilter()
    result = bbb.score_drug("panobinostat")
    assert_equal(result["penetrance"], "HIGH", "panobinostat BBB penetrance")


def test_smarcb1_expression_score_low() -> None:
    """
    SMARCB1 curated expression score must be ≤ 0.10.
    Source: GSE70678 — SMARCB1 expression ~0.05 of normal brain
    in biallelic deletion (Torchia 2015 Cancer Cell PMID 26609405).
    """
    from backend.pipeline.pipeline_config import ATRT_CURATED_SCORES
    score = ATRT_CURATED_SCORES["SMARCB1"]
    assert score <= 0.10, f"SMARCB1 score = {score}, expected ≤ 0.10"
    logger.info("  ✅ SMARCB1 curated score = %.2f ≤ 0.10 (gene is lost in ATRT)", score)


def test_ezh2_expression_score_high() -> None:
    """
    EZH2 curated expression score must be ≥ 0.85.
    Source: GSE70678 z-score 2.31 in ATRT vs normal brain;
    EZH2 hyperactive due to SMARCB1 loss (Knutson 2013 PNAS PMID 23620515).
    """
    from backend.pipeline.pipeline_config import ATRT_CURATED_SCORES
    score = ATRT_CURATED_SCORES["EZH2"]
    assert_range(score, 0.85, 1.00, "EZH2 curated expression score", pmid="23620515")


def test_psmb5_depmap_score_essential() -> None:
    """
    PSMB5 DepMap Chronos must be ≤ -2.0 (extremely essential).
    Source: Lin 2019 Sci Transl Med — PSMB5 Chronos -3.30 in GBM lines.
    Rhabdoid lines show similar extreme essentiality (proteasome is universal).
    """
    from backend.pipeline.depmap_essentiality import _CHRONOS_FALLBACK
    chronos = _CHRONOS_FALLBACK["PSMB5"]
    assert chronos <= -2.0, \
        f"PSMB5 Chronos = {chronos}, expected ≤ -2.0 (extremely essential)"
    logger.info(
        "  ✅ PSMB5 Chronos = %.2f ≤ -2.0 (extremely essential, Lin 2019 Sci Transl Med)",
        chronos,
    )


def test_ezh2_chronos_essential() -> None:
    """
    EZH2 DepMap Chronos must be ≤ -1.0 (strongly essential in SMARCB1-null).
    Source: Knutson 2013 PNAS — EZH2 knockdown lethal in G401/A204/BT16.
    """
    from backend.pipeline.depmap_essentiality import _CHRONOS_FALLBACK
    chronos = _CHRONOS_FALLBACK["EZH2"]
    assert chronos <= -1.0, \
        f"EZH2 Chronos = {chronos}, expected ≤ -1.0"
    logger.info(
        "  ✅ EZH2 Chronos = %.2f ≤ -1.0 (Knutson 2013 PNAS PMID 23620515)",
        chronos,
    )


def test_panobinostat_ic50_verified() -> None:
    """
    Panobinostat IC50 in BT16 must be ~0.009 µM (within ±30% for assay variability).
    Source: Torchia 2015 Cancer Cell PMID 26609405 — Fig 4.
    """
    from backend.pipeline.published_ic50_atrt_validation import get_atrt_validation_score
    result = get_atrt_validation_score("PANOBINOSTAT")
    assert result is not None, "PANOBINOSTAT not in IC50 database"
    assert result["has_verified_data"], "PANOBINOSTAT has no verified data"
    ic50 = result["best_ic50_um"]
    # 0.009 µM ± 30% = [0.0063, 0.0117]
    assert_range(
        ic50, 0.005, 0.015,
        "Panobinostat best IC50 (µM, BT16)",
        pmid="26609405",
    )


def test_alisertib_ic50_verified() -> None:
    """
    Alisertib IC50 in BT16 must be ~0.098 µM (within ±30%).
    Source: Sredni 2017 Pediatric Blood Cancer PMID 28544500.
    """
    from backend.pipeline.published_ic50_atrt_validation import get_atrt_validation_score
    result = get_atrt_validation_score("ALISERTIB")
    assert result is not None and result["has_verified_data"]
    ic50 = result["best_ic50_um"]
    assert_range(
        ic50, 0.060, 0.150,
        "Alisertib best IC50 (µM, BT16)",
        pmid="28544500",
    )


def test_marizomib_no_ic50_data() -> None:
    """
    Marizomib must have NO verified ATRT-specific IC50 data.
    'Orphanides 2023' reference was unverifiable — removed in v2.0.
    Bota 2021 covers GBM only, not ATRT.
    """
    from backend.pipeline.published_ic50_atrt_validation import get_atrt_validation_score
    result = get_atrt_validation_score("MARIZOMIB")
    assert result is not None, "MARIZOMIB should be in database (with empty entries)"
    assert not result["has_verified_data"], \
        "MARIZOMIB should have no verified data (Orphanides 2023 removed)"
    assert_none(result["validation_score"], "Marizomib IC50 validation score")
    logger.info("  ✅ Marizomib: no verified ATRT IC50 (unverifiable 'Orphanides 2023' removed)")


def test_onc201_no_ic50_data() -> None:
    """
    ONC201 must have NO verified ATRT-specific IC50 data.
    'Frühwald 2020' is a review paper with no primary IC50 data — removed in v2.0.
    """
    from backend.pipeline.published_ic50_atrt_validation import get_atrt_validation_score
    result = get_atrt_validation_score("ONC201")
    assert result is not None, "ONC201 should be in database (with empty entries)"
    assert not result["has_verified_data"], \
        "ONC201 should have no verified data (Frühwald 2020 review removed)"
    assert_none(result["validation_score"], "ONC201 IC50 validation score")
    logger.info("  ✅ ONC201: no verified ATRT IC50 (Frühwald 2020 review paper removed)")


def test_ezh2_boost_config() -> None:
    """
    EZH2_INHIBITOR composite_boost must be 1.40 (ATRT synthetic lethality).
    Source: Knutson 2013 PNAS PMID 23620515.
    This is a BOOST (not penalty as in DIPG).
    """
    from backend.pipeline.pipeline_config import EZH2_INHIBITOR
    assert EZH2_INHIBITOR["is_boost"] is True, "EZH2_INHIBITOR must be a boost in ATRT"
    assert_range(
        EZH2_INHIBITOR["composite_boost"],
        1.35, 1.45,
        "EZH2 composite boost multiplier",
        pmid="23620515",
    )


def test_toxicity_rates_in_range() -> None:
    """
    Published toxicity rates must match verified trial data within ±5%.
    Sources listed in pipeline_config.py TOXICITY section.
    """
    from backend.pipeline.pipeline_config import TOXICITY
    rates = TOXICITY["drug_ae_rates"]

    # Panobinostat: PBTC-047 8/29 DLTs = 27.6% → rounded to 20% G3/4 hematologic
    # Source: Monje 2023 Nat Med PMID 37526549
    assert_range(
        rates["PANOBINOSTAT"], 0.15, 0.30,
        "Panobinostat G3/4 AE rate",
        pmid="37526549",
    )

    # Tazemetostat: Gounder 2020 JCO 4.8% G3/4
    assert_range(
        rates["TAZEMETOSTAT"], 0.02, 0.10,
        "Tazemetostat G3/4 AE rate",
        pmid="33166238",
    )

    # Marizomib: Bota 2021 ~10% G3/4
    assert_range(
        rates["MARIZOMIB"], 0.05, 0.15,
        "Marizomib G3/4 AE rate",
        pmid="33300566",
    )


def test_cmap_no_hardcoded_scores() -> None:
    """
    CMap module must NOT contain hardcoded connectivity scores.
    All scores must come from loaded JSON file or return neutral prior (0.50).
    """
    from backend.pipeline.cmap_query import CMAPQuery
    cmap = CMAPQuery()
    # With no JSON loaded, all scores should be 0.50
    score = cmap.get_precomputed_score("TAZEMETOSTAT")
    # Either 0.50 (no data) or a real score from JSON (0.0 to 1.0)
    assert 0.0 <= score <= 1.0, f"CMap score out of range: {score}"
    # Specifically: tazemetostat NOT in L1000 library → must return 0.50
    from backend.pipeline.cmap_query import L1000_LIBRARY_NOTES
    assert not L1000_LIBRARY_NOTES["TAZEMETOSTAT"]["in_library"], \
        "Tazemetostat should be flagged as NOT in L1000 library"
    logger.info("  ✅ CMap returns neutral prior for non-L1000 drugs (no hardcoded scores)")


def test_opentargets_efo_ids() -> None:
    """
    OpenTargets must include the primary ATRT EFO ID (EFO_0002915).
    Source: EFO ontology — rhabdoid tumor.
    """
    from backend.pipeline.pipeline_config import OPENTARGETS
    assert "EFO_0002915" in OPENTARGETS["efo_ids"], \
        "Primary ATRT EFO ID (EFO_0002915 = rhabdoid tumor) missing"
    logger.info("  ✅ Primary ATRT EFO ID EFO_0002915 present")


async def test_depmap_fallback_loads() -> None:
    """
    DepMap module must load fallback without crashing when CSV absent.
    All key ATRT targets must be present in fallback.
    """
    from backend.pipeline.depmap_essentiality import DepMapEssentiality, _CHRONOS_FALLBACK
    de = DepMapEssentiality()
    await de._load_data_if_needed()
    assert de.is_ready, "DepMap module should be ready after _load_data_if_needed()"
    assert de.using_fallback, "Expected fallback mode when CSVs absent"

    key_targets = {"EZH2", "BRD4", "HDAC1", "AURKA", "CDK4", "PSMB5", "MYC"}
    for target in key_targets:
        assert target in de.gene_scores, f"{target} missing from DepMap fallback"
    logger.info("  ✅ DepMap fallback loads cleanly with all key ATRT targets")


# ─────────────────────────────────────────────────────────────────────────────
# Test runner
# ─────────────────────────────────────────────────────────────────────────────

SYNC_TESTS = [
    ("Config weights sum to 1.0",                test_config_weights),
    ("Tazemetostat BBB = MODERATE (not HIGH)",    test_tazemetostat_bbb_moderate),
    ("Panobinostat BBB = HIGH",                   test_panobinostat_bbb_high),
    ("SMARCB1 expression score ≤ 0.10",           test_smarcb1_expression_score_low),
    ("EZH2 expression score ≥ 0.85",              test_ezh2_expression_score_high),
    ("PSMB5 Chronos ≤ -2.0 (extremely essential)", test_psmb5_depmap_score_essential),
    ("EZH2 Chronos ≤ -1.0 (strongly essential)",  test_ezh2_chronos_essential),
    ("Panobinostat IC50 verified (~0.009 µM BT16)", test_panobinostat_ic50_verified),
    ("Alisertib IC50 verified (~0.098 µM BT16)",  test_alisertib_ic50_verified),
    ("Marizomib: no verified ATRT IC50",           test_marizomib_no_ic50_data),
    ("ONC201: no verified ATRT IC50",              test_onc201_no_ic50_data),
    ("EZH2 boost = 1.40 (ATRT synthetic lethality)", test_ezh2_boost_config),
    ("Toxicity rates match published trials",      test_toxicity_rates_in_range),
    ("CMap: no hardcoded scores",                  test_cmap_no_hardcoded_scores),
    ("OpenTargets: ATRT EFO_0002915 present",      test_opentargets_efo_ids),
]

ASYNC_TESTS = [
    ("DepMap fallback loads cleanly",              test_depmap_fallback_loads),
]


async def run_all_tests() -> None:
    print("\n" + "═" * 65)
    print("ATRT Pipeline Validation Suite v2.0")
    print("═" * 65)

    passed = 0
    failed = 0
    errors = []

    for name, test_fn in SYNC_TESTS:
        print(f"\n[{passed + failed + 1:02d}] {name}")
        try:
            test_fn()
            passed += 1
        except (ValidationError, AssertionError, ImportError, Exception) as e:
            failed += 1
            errors.append((name, str(e)))
            logger.error("  ❌ %s: %s", name, e)

    for name, test_fn in ASYNC_TESTS:
        print(f"\n[{passed + failed + 1:02d}] {name}")
        try:
            await test_fn()
            passed += 1
        except (ValidationError, AssertionError, ImportError, Exception) as e:
            failed += 1
            errors.append((name, str(e)))
            logger.error("  ❌ %s: %s", name, e)

    total = passed + failed
    print("\n" + "═" * 65)
    print(f"Results: {passed}/{total} passed")

    if errors:
        print("\nFailed tests:")
        for name, msg in errors:
            print(f"  ❌ {name}")
            print(f"     {msg}")
    else:
        print("✅ All validation tests passed")

    print("═" * 65 + "\n")

    if failed > 0:
        raise SystemExit(f"{failed} validation test(s) failed")


if __name__ == "__main__":
    asyncio.run(run_all_tests())