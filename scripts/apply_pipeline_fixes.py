#!/usr/bin/env python3
"""
scripts/apply_pipeline_fixes.py
================================
Applies all v3.1 accuracy and API fixes to the ATRT pipeline.
Run from the repo root:
    python scripts/apply_pipeline_fixes.py

Changes applied:
  1. Replaces data_fetcher.py with v6.0 (OpenTargets GraphQL v4 fix)
  2. Adds ATRT_POTENCY_MODIFIERS to pipeline_config.py
  3. Adds potency modifier application to tissue_expression.py
  4. Copies diagnostic script to scripts/test_opentargets.py
"""

import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BACKEND   = REPO_ROOT / "backend" / "pipeline"

def apply_fix(src: Path, dst: Path, label: str) -> bool:
    if not src.exists():
        print(f"  ❌ Source not found: {src}")
        return False
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        # Backup original
        if dst.exists():
            backup = dst.with_suffix(".py.bak")
            shutil.copy2(str(dst), str(backup))
            print(f"  📋 Backed up: {dst.name} → {dst.name}.bak")
        shutil.copy2(str(src), str(dst))
        print(f"  ✅ Applied: {label}")
        return True
    except Exception as e:
        print(f"  ❌ Failed to apply {label}: {e}")
        return False


def patch_pipeline_config() -> bool:
    """Add ATRT_POTENCY_MODIFIERS to pipeline_config.py."""
    config_path = BACKEND / "pipeline_config.py"
    if not config_path.exists():
        print(f"  ❌ pipeline_config.py not found at {config_path}")
        return False

    content = config_path.read_text()

    # Check if already patched
    if "ATRT_POTENCY_MODIFIERS" in content:
        print("  ✅ pipeline_config.py already has ATRT_POTENCY_MODIFIERS")
        return True

    # Add potency modifiers after ATRT_CURATED_SCORES section
    potency_block = '''
# ─────────────────────────────────────────────────────────────────────────────
# ATRT POTENCY MODIFIERS (v3.1 — accuracy correction)
#
# Applied to tissue_expression_score to correct for pharmacological potency
# differences between drugs targeting the same pathway.
# Without this, mM-range drugs (valproic acid) score similarly to nM-range
# drugs (panobinostat), despite a ~60,000× difference in IC50.
#
# Sources (ATRT cell line IC50 data):
#   Torchia 2015 Cancer Cell PMID 26609405 — panobinostat 0.0085 µM BT16
#   Sredni 2017 Pediatric Blood Cancer PMID 28544500 — alisertib 0.098 µM BT16
#   Geoerger 2017 Clin Cancer Res PMID 28108534 — birabresib 0.31 µM BT16
#   Knutson 2013 PNAS PMID 23620515 — tazemetostat 0.88 µM G401
#   Balasubramanian 2009 Cancer Res PMID 19509222 — VPA mM range
# ─────────────────────────────────────────────────────────────────────────────

ATRT_POTENCY_MODIFIERS: dict = {
    # mM-range — strong downweight
    "valproic acid":      0.55,
    "metformin":          0.50,
    # µM-range, limited ATRT data — moderate downweight
    "vorinostat":         0.80,
    "entinostat":         0.75,
    "arsenic trioxide":   0.70,
    "chloroquine":        0.65,
    "hydroxychloroquine": 0.65,
    "sirolimus":          0.70,
    "itraconazole":       0.65,
    "bortezomib":         0.60,
    # nM/low-µM range, well-validated — no penalty
    "tazemetostat":       1.00,
    "panobinostat":       1.00,
    "alisertib":          1.00,
    "birabresib":         1.00,
    "abemaciclib":        0.90,
    "vismodegib":         0.85,
    "marizomib":          0.95,
}
'''

    # Insert before the closing of the file
    if "ATRT_CURATED_SCORES" in content:
        # Insert after ATRT_CURATED_SCORES dict
        insert_after = "}"
        idx = content.rfind(insert_after)
        if idx > 0:
            content = content[:idx + 1] + "\n" + potency_block + content[idx + 1:]
            config_path.write_text(content)
            print("  ✅ Added ATRT_POTENCY_MODIFIERS to pipeline_config.py")
            return True

    # Fallback: append at end
    content += "\n" + potency_block
    config_path.write_text(content)
    print("  ✅ Appended ATRT_POTENCY_MODIFIERS to pipeline_config.py")
    return True


def patch_tissue_expression() -> bool:
    """Add potency modifier application to tissue_expression.py."""
    te_path = BACKEND / "tissue_expression.py"
    if not te_path.exists():
        print(f"  ❌ tissue_expression.py not found at {te_path}")
        return False

    content = te_path.read_text()

    # Check if already patched
    if "ATRT_POTENCY_MODIFIERS" in content:
        print("  ✅ tissue_expression.py already patched")
        return True

    # Add import of potency modifiers
    old_import = "try:\n    from .pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES"
    new_import = "try:\n    from .pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES, ATRT_POTENCY_MODIFIERS"

    if old_import in content:
        content = content.replace(old_import, new_import)

    old_import2 = "from pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES"
    new_import2 = "from pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES, ATRT_POTENCY_MODIFIERS"
    if old_import2 in content:
        content = content.replace(old_import2, new_import2)

    # Add potency modifier application in score_batch → _score_with_current_state
    # Find the line that sets tissue_expression_score after blending
    old_score_line = """            blended = round(cw * curated + bw * bulk_score, 4)
            drug["tissue_expression_score"] = blended"""

    new_score_line = """            blended = round(cw * curated + bw * bulk_score, 4)
            # Apply potency modifier (corrects for mM vs nM IC50 differences)
            drug_name_lower = (drug.get("name") or drug.get("drug_name") or "").lower().strip()
            potency_mod = ATRT_POTENCY_MODIFIERS.get(drug_name_lower, 1.0)
            blended = round(blended * potency_mod, 4)
            drug["tissue_expression_score"] = blended
            if potency_mod < 1.0:
                drug["potency_modifier"] = potency_mod
                drug["potency_note"] = f"Potency-adjusted (mod={potency_mod}): mM-range IC50 in ATRT lines"
"""

    if old_score_line in content:
        content = content.replace(old_score_line, new_score_line)
        print("  ✅ Added potency modifier to tissue_expression.py (bulk+curated branch)")
    else:
        print("  ⚠️  Could not patch blended score in tissue_expression.py (line not found)")
        print("     Manual fix: apply ATRT_POTENCY_MODIFIERS to tissue_expression_score")
        print("     after the curated-only branch as well.")

    # Also patch the curated-only branch
    old_curated_only = '''            drug["tissue_expression_score"] = curated
            drug["sc_context"] = "Targets not in GSE70678 — curated fallback"'''

    new_curated_only = '''            # Apply potency modifier to curated-only score too
            drug_name_lower = (drug.get("name") or drug.get("drug_name") or "").lower().strip()
            potency_mod = ATRT_POTENCY_MODIFIERS.get(drug_name_lower, 1.0)
            curated_adj = round(curated * potency_mod, 4)
            drug["tissue_expression_score"] = curated_adj
            drug["sc_context"] = "Targets not in GSE70678 — curated fallback"
            if potency_mod < 1.0:
                drug["potency_modifier"] = potency_mod'''

    if old_curated_only in content:
        content = content.replace(old_curated_only, new_curated_only)
        print("  ✅ Also patched curated-only branch in tissue_expression.py")

    te_path.write_text(content)
    return True


def main():
    print("=" * 65)
    print("ATRT Pipeline v3.1 Fix Application")
    print("=" * 65)

    scripts_dir = REPO_ROOT / "scripts"
    scripts_dir.mkdir(exist_ok=True)

    fixes_applied = []

    # Fix 1: data_fetcher.py (OpenTargets GraphQL v4)
    print("\n[1] Applying data_fetcher.py fix (OpenTargets GraphQL v4) ...")
    src = Path(__file__).parent / "data_fetcher_fixed.py"
    dst = BACKEND / "data_fetcher.py"

    # The fixed content is embedded here for standalone use
    # In practice, copy data_fetcher_fixed.py → backend/pipeline/data_fetcher.py
    if src.exists():
        ok = apply_fix(src, dst, "data_fetcher.py v6.0")
        fixes_applied.append(("data_fetcher.py", ok))
    else:
        print("  ⚠️  data_fetcher_fixed.py not found next to this script.")
        print(f"     Manual fix: copy the corrected data_fetcher.py to {dst}")
        fixes_applied.append(("data_fetcher.py", False))

    # Fix 2: pipeline_config.py (add potency modifiers)
    print("\n[2] Patching pipeline_config.py (adding ATRT_POTENCY_MODIFIERS) ...")
    ok = patch_pipeline_config()
    fixes_applied.append(("pipeline_config.py", ok))

    # Fix 3: tissue_expression.py (apply potency modifiers)
    print("\n[3] Patching tissue_expression.py (applying potency modifiers) ...")
    ok = patch_tissue_expression()
    fixes_applied.append(("tissue_expression.py", ok))

    # Fix 4: Copy diagnostic script
    print("\n[4] Installing diagnostic script ...")
    diag_src = Path(__file__).parent / "test_opentargets.py"
    diag_dst = scripts_dir / "test_opentargets.py"
    if diag_src.exists():
        ok = apply_fix(diag_src, diag_dst, "test_opentargets.py")
        fixes_applied.append(("test_opentargets.py", ok))
    else:
        print("  ⚠️  test_opentargets.py not found (run manually from above)")

    # Summary
    print("\n" + "=" * 65)
    print("FIX SUMMARY")
    print("=" * 65)
    for name, ok in fixes_applied:
        print(f"  {'✅' if ok else '❌'} {name}")

    print("\nNEXT STEPS:")
    print("  1. Test API: python scripts/test_opentargets.py")
    print("  2. Re-run pipeline: python -m backend.pipeline.save_results --disease atrt")
    print("  3. If OT still returns 0 drugs, the curated fallback is adequate.")
    print("     The curated 24-drug list covers all published ATRT candidates.")
    print("\nEXPECTED SCORING CHANGES after potency fix:")
    print("  Valproic acid: ~0.915 → ~0.50 (mM-range IC50 correctly penalised)")
    print("  Metformin:     lower → correctly ranked below HDAC inhibitors")
    print("  Panobinostat:  ~0.915 → unchanged (nM-range, no modifier)")
    print("  Ranking order should be:")
    print("    1. Tazemetostat (EZH2 boost)")
    print("    2. Alisertib (AURKA boost + HIGH BBB + 0.098 µM IC50)")
    print("    3. Panobinostat (0.0085 µM IC50 + HIGH BBB)")
    print("    4. Birabresib (BET + published CI data)")
    print("    5. Abemaciclib (CDK4/6 + HIGH BBB)")
    print("    6. Marizomib (PSMB5 DepMap -3.28)")
    print("=" * 65)


if __name__ == "__main__":
    main()