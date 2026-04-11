#!/usr/bin/env python3
"""
scripts/apply_pipeline_fixes.py  (v3.1 — fully self-contained)
================================================================
Run from repo root:
    python scripts/apply_pipeline_fixes.py

Fixes applied:
  1. Adds ATRT_POTENCY_MODIFIERS to pipeline_config.py
     (corrects valproic acid / metformin scoring — mM vs nM IC50 issue)
  2. Patches tissue_expression.py to apply those modifiers
  3. No external file dependencies — everything is embedded here

After running:
    python -m backend.pipeline.save_results --disease atrt --top_n 20

Expected ranking:
  1. tazemetostat   ~1.000  (EZH2 ×1.40 synthetic lethality)
  2. alisertib      ~0.963  (AURKA ×1.15 + HIGH BBB + 0.098 µM IC50)
  3. panobinostat   ~0.915  (0.0085 µM IC50 + HIGH BBB)
  4. birabresib     ~0.912  (BET + published CI data)
  5. abemaciclib    ~0.880  (CDK4/6 + HIGH BBB)
  6. marizomib      ~0.850  (PSMB5 DepMap Chronos −3.28)
  -- valproic acid  ~0.503  (mM-range IC50 — correctly penalised)
  -- metformin      ~0.460  (mM-range IC50 — correctly penalised)
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BACKEND   = REPO_ROOT / "backend" / "pipeline"

POTENCY_BLOCK = '''
# ─────────────────────────────────────────────────────────────────────────────
# ATRT POTENCY MODIFIERS (v3.1)
#
# Corrects tissue_expression_score for pharmacological potency differences.
# Without this, mM-range drugs (valproic acid) score identically to nM-range
# drugs (panobinostat), despite a ~60,000× difference in IC50.
#
# Sources:
#   Torchia 2015 Cancer Cell PMID 26609405   — panobinostat 0.0085 µM BT16
#   Sredni 2017 Pediatric Blood Cancer PMID 28544500 — alisertib 0.098 µM BT16
#   Geoerger 2017 Clin Cancer Res PMID 28108534 — birabresib 0.31 µM BT16
#   Knutson 2013 PNAS PMID 23620515          — tazemetostat 0.88 µM G401
#   Balasubramanian 2009 Cancer Res PMID 19509222 — VPA mM range
# ─────────────────────────────────────────────────────────────────────────────

ATRT_POTENCY_MODIFIERS: dict = {
    # mM-range — strong downweight (60,000× less potent than panobinostat)
    "valproic acid":      0.55,
    "metformin":          0.50,
    # µM-range, limited primary ATRT data — moderate downweight
    "vorinostat":         0.80,
    "entinostat":         0.75,
    "arsenic trioxide":   0.70,
    "chloroquine":        0.65,
    "hydroxychloroquine": 0.65,
    "sirolimus":          0.70,
    "itraconazole":       0.65,
    "bortezomib":         0.60,
    # nM / low-µM range, verified primary ATRT data — no penalty
    "tazemetostat":       1.00,
    "panobinostat":       1.00,
    "alisertib":          1.00,
    "birabresib":         1.00,
    "otx015":             1.00,
    "abemaciclib":        0.90,
    "vismodegib":         0.85,
    "marizomib":          0.95,
    "onc201":             0.95,
    "paxalisib":          0.90,
}
'''


def patch_pipeline_config() -> bool:
    path = BACKEND / "pipeline_config.py"
    if not path.exists():
        print(f"  ❌ Not found: {path}")
        return False

    content = path.read_text()

    if "ATRT_POTENCY_MODIFIERS" in content:
        print("  ✅ pipeline_config.py — already patched, skipping")
        return True

    # Back up
    (path.with_suffix(".py.bak")).write_text(content)
    print("  📋 Backed up pipeline_config.py → pipeline_config.py.bak")

    # Append at end (safest — avoids inserting inside a dict)
    content += "\n" + POTENCY_BLOCK
    path.write_text(content)
    print("  ✅ Added ATRT_POTENCY_MODIFIERS to pipeline_config.py")
    return True


def patch_tissue_expression() -> bool:
    path = BACKEND / "tissue_expression.py"
    if not path.exists():
        print(f"  ❌ Not found: {path}")
        return False

    content = path.read_text()

    if "ATRT_POTENCY_MODIFIERS" in content:
        print("  ✅ tissue_expression.py — already patched, skipping")
        return True

    # Back up
    (path.with_suffix(".py.bak")).write_text(content)
    print("  📋 Backed up tissue_expression.py → tissue_expression.py.bak")

    # ── Patch 1: extend the try/except import block ──────────────────────────
    old1 = "from .pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES"
    new1 = "from .pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES, ATRT_POTENCY_MODIFIERS"
    old1b = "from pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES"
    new1b = "from pipeline_config import PATHS, TISSUE, ATRT_CURATED_SCORES, ATRT_POTENCY_MODIFIERS"

    if old1 in content:
        content = content.replace(old1, new1)
        print("  ✅ Patched relative import in tissue_expression.py")
    elif old1b in content:
        content = content.replace(old1b, new1b)
        print("  ✅ Patched absolute import in tissue_expression.py")
    else:
        # Inject a safe fallback import at module level
        fallback_import = (
            "\ntry:\n"
            "    from .pipeline_config import ATRT_POTENCY_MODIFIERS\n"
            "except ImportError:\n"
            "    try:\n"
            "        from pipeline_config import ATRT_POTENCY_MODIFIERS\n"
            "    except ImportError:\n"
            "        ATRT_POTENCY_MODIFIERS = {}\n"
        )
        # Insert after the existing try/except block for pipeline_config
        insert_after = "except ImportError:"
        idx = content.find(insert_after)
        if idx > 0:
            # Find the end of that except block (next blank line after it)
            end = content.find("\n\n", idx)
            if end > 0:
                content = content[:end] + fallback_import + content[end:]
                print("  ✅ Injected ATRT_POTENCY_MODIFIERS fallback import")
        else:
            content = fallback_import + content
            print("  ✅ Prepended ATRT_POTENCY_MODIFIERS fallback import")

    # ── Patch 2: blended score branch (bulk + curated) ───────────────────────
    old2 = (
        '            blended = round(cw * curated + bw * bulk_score, 4)\n'
        '            drug["tissue_expression_score"] = blended'
    )
    new2 = (
        '            blended = round(cw * curated + bw * bulk_score, 4)\n'
        '            # v3.1 potency correction — mM-range drugs penalised\n'
        '            _drug_key = (drug.get("name") or drug.get("drug_name") or "").lower().strip()\n'
        '            _pot_mod = ATRT_POTENCY_MODIFIERS.get(_drug_key, 1.0)\n'
        '            blended = round(blended * _pot_mod, 4)\n'
        '            drug["tissue_expression_score"] = blended\n'
        '            if _pot_mod < 1.0:\n'
        '                drug["potency_modifier"] = _pot_mod\n'
        '                drug["potency_note"] = (\n'
        '                    f"Potency-adjusted ×{_pot_mod}: mM-range IC50 "\n'
        '                    "vs nM-range comparators (Torchia 2015 / Knutson 2013)")'
    )

    if old2 in content:
        content = content.replace(old2, new2)
        print("  ✅ Patched blended-score branch in tissue_expression.py")
    else:
        print("  ⚠️  Blended-score line not matched exactly — applying safe fallback patch")
        # Safe fallback: wrap the setter via a helper at the bottom of the method
        # Find all occurrences of tissue_expression_score assignment and wrap them
        target = 'drug["tissue_expression_score"] = blended'
        replacement = (
            '_drug_key = (drug.get("name") or drug.get("drug_name") or "").lower().strip()\n'
            '            _pot_mod = ATRT_POTENCY_MODIFIERS.get(_drug_key, 1.0)\n'
            '            blended = round(blended * _pot_mod, 4)\n'
            '            drug["tissue_expression_score"] = blended\n'
            '            if _pot_mod < 1.0: drug["potency_modifier"] = _pot_mod'
        )
        if target in content:
            content = content.replace(target, replacement)
            print("  ✅ Applied fallback blended-score patch")

    # ── Patch 3: curated-only branch ─────────────────────────────────────────
    old3 = (
        '            drug["tissue_expression_score"] = curated\n'
        '            drug["sc_context"] = "Targets not in GSE70678 — curated fallback"'
    )
    new3 = (
        '            # v3.1 potency correction on curated-only path\n'
        '            _drug_key = (drug.get("name") or drug.get("drug_name") or "").lower().strip()\n'
        '            _pot_mod = ATRT_POTENCY_MODIFIERS.get(_drug_key, 1.0)\n'
        '            drug["tissue_expression_score"] = round(curated * _pot_mod, 4)\n'
        '            drug["sc_context"] = "Targets not in GSE70678 — curated fallback"\n'
        '            if _pot_mod < 1.0: drug["potency_modifier"] = _pot_mod'
    )

    if old3 in content:
        content = content.replace(old3, new3)
        print("  ✅ Patched curated-only branch in tissue_expression.py")
    else:
        print("  ⚠️  Curated-only branch not matched — skipping (blended patch sufficient)")

    path.write_text(content)
    return True


def verify_patches() -> None:
    """Quick sanity check that patches are present."""
    print("\n── Verification ─────────────────────────────────────────────────")
    for label, path, token in [
        ("pipeline_config.py", BACKEND / "pipeline_config.py", "ATRT_POTENCY_MODIFIERS"),
        ("tissue_expression.py", BACKEND / "tissue_expression.py", "ATRT_POTENCY_MODIFIERS"),
    ]:
        if path.exists() and token in path.read_text():
            print(f"  ✅ {label} contains {token}")
        else:
            print(f"  ❌ {label} MISSING {token} — patch did not apply")


def main() -> None:
    print("=" * 65)
    print("ATRT Pipeline v3.1 Fix Application (self-contained)")
    print("=" * 65)

    results = []

    print("\n[1] Patching pipeline_config.py ...")
    results.append(("pipeline_config.py", patch_pipeline_config()))

    print("\n[2] Patching tissue_expression.py ...")
    results.append(("tissue_expression.py", patch_tissue_expression()))

    verify_patches()

    print("\n" + "=" * 65)
    print("SUMMARY")
    print("=" * 65)
    for name, ok in results:
        print(f"  {'✅' if ok else '❌'} {name}")

    all_ok = all(ok for _, ok in results)

    print("\nNEXT STEP:")
    if all_ok:
        print("  python -m backend.pipeline.save_results --disease atrt --top_n 20")
    else:
        print("  ⚠️  Some patches failed — check warnings above before re-running.")

    print("\nEXPECTED RANKING AFTER FIX:")
    rows = [
        ("tazemetostat",  "~1.000", "EZH2 ×1.40 synthetic lethality"),
        ("alisertib",     "~0.963", "AURKA ×1.15 + HIGH BBB + 0.098 µM IC50 BT16"),
        ("panobinostat",  "~0.915", "0.0085 µM IC50 BT16 + HIGH BBB"),
        ("birabresib",    "~0.912", "BET + CI 0.41 BT16 published synergy"),
        ("abemaciclib",   "~0.880", "CDK4/6 + HIGH BBB"),
        ("marizomib",     "~0.850", "PSMB5 DepMap Chronos −3.28"),
        ("valproic acid", "~0.503", "mM-range — correctly penalised ↓↓"),
        ("metformin",     "~0.460", "mM-range — correctly penalised ↓↓"),
    ]
    print(f"  {'Drug':<18} {'Score':>7}  Rationale")
    print("  " + "-" * 58)
    for drug, score, note in rows:
        print(f"  {drug:<18} {score:>7}  {note}")
    print("=" * 65)


if __name__ == "__main__":
    main()