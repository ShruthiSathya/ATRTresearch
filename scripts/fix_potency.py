#!/usr/bin/env python3
"""
scripts/fix_potency.py
======================
Direct fix for the VPA/metformin potency scoring bug.
Run from repo root:
    python scripts/fix_potency.py
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TE_PATH   = REPO_ROOT / "backend" / "pipeline" / "tissue_expression.py"
PC_PATH   = REPO_ROOT / "backend" / "pipeline" / "pipeline_config.py"

POTENCY_MAP = {
    "valproic acid":      0.55,
    "metformin":          0.50,
    "metformin hcl":      0.50,
    "vorinostat":         0.80,
    "entinostat":         0.75,
    "arsenic trioxide":   0.70,
    "chloroquine":        0.65,
    "hydroxychloroquine": 0.65,
    "sirolimus":          0.70,
    "itraconazole":       0.65,
    "bortezomib":         0.60,
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


# ── Step 1: add/update ATRT_POTENCY_MODIFIERS in pipeline_config.py ──────────

def fix_pipeline_config() -> bool:
    content = PC_PATH.read_text()

    if "ATRT_POTENCY_MODIFIERS" in content:
        print("  pipeline_config.py — already has ATRT_POTENCY_MODIFIERS, updating values...")
        # Remove old block and re-append fresh one
        start = content.find("\nATRT_POTENCY_MODIFIERS")
        if start != -1:
            content = content[:start]

    PC_PATH.with_suffix(".py.bak").write_text(PC_PATH.read_text())

    lines = ["", "ATRT_POTENCY_MODIFIERS: dict = {"]
    for k, v in POTENCY_MAP.items():
        lines.append(f'    "{k}": {v},')
    lines.append("}")
    content += "\n".join(lines) + "\n"
    PC_PATH.write_text(content)
    print("  ✅ pipeline_config.py — ATRT_POTENCY_MODIFIERS written")
    return True


# ── Step 2: replace score_batch in tissue_expression.py with patched version ─

def fix_tissue_expression() -> bool:
    content = TE_PATH.read_text()

    if "POTENCY_MAP_INLINE" in content:
        print("  tissue_expression.py — already patched (POTENCY_MAP_INLINE found)")
        return True

    TE_PATH.with_suffix(".py.bak").write_text(content)

    # Build inline potency map string for injection (no import needed — avoids
    # any import-path issues entirely)
    map_lines = ["    _POTENCY = {"]
    for k, v in POTENCY_MAP.items():
        map_lines.append(f'        "{k}": {v},')
    map_lines.append("    }  # POTENCY_MAP_INLINE")
    potency_map_str = "\n".join(map_lines)

    # The patch: add a _apply_potency() helper method to TissueExpressionScorer
    # and call it from _score_with_current_state.
    # We inject it by appending to the class definition before the last method.

    # Find _score_with_current_state and wrap every tissue_expression_score assignment
    # Strategy: find all lines setting drug["tissue_expression_score"] and wrap them.

    lines = content.split("\n")
    new_lines = []
    patched_count = 0

    i = 0
    while i < len(lines):
        line = lines[i]

        # Inject potency map once at the start of _score_with_current_state
        if "def _score_with_current_state(self" in line:
            new_lines.append(line)
            # Find the first line of the method body (after docstring if any)
            i += 1
            while i < len(lines) and lines[i].strip() == "":
                new_lines.append(lines[i])
                i += 1
            # Insert inline potency map
            new_lines.append(potency_map_str)
            new_lines.append("")
            continue

        # Wrap every assignment to tissue_expression_score
        if 'drug["tissue_expression_score"]' in line and "=" in line and "potency" not in line:
            indent = len(line) - len(line.lstrip())
            spaces = " " * indent
            # Extract the value being assigned
            # e.g.  drug["tissue_expression_score"] = blended
            # or    drug["tissue_expression_score"] = curated
            rhs = line.split("=", 1)[1].strip()
            new_lines.append(f'{spaces}_drug_key = (drug.get("name") or drug.get("drug_name") or "").lower().strip()')
            new_lines.append(f'{spaces}_pot = _POTENCY.get(_drug_key, 1.0)')
            new_lines.append(f'{spaces}drug["tissue_expression_score"] = round({rhs} * _pot, 4)')
            new_lines.append(f'{spaces}if _pot < 1.0: drug["potency_modifier"] = _pot')
            patched_count += 1
            i += 1
            continue

        new_lines.append(line)
        i += 1

    if patched_count == 0:
        print("  ❌ No tissue_expression_score assignments found — check file manually")
        return False

    TE_PATH.write_text("\n".join(new_lines))
    print(f"  ✅ tissue_expression.py — patched {patched_count} score assignment(s)")
    return True


# ── Step 3: verify ────────────────────────────────────────────────────────────

def verify() -> None:
    print("\n── Verification ──────────────────────────────────────────────────")

    # Check pipeline_config
    pc = PC_PATH.read_text()
    if "ATRT_POTENCY_MODIFIERS" in pc and '"valproic acid": 0.55' in pc:
        print("  ✅ pipeline_config.py — ATRT_POTENCY_MODIFIERS present with correct VPA value")
    else:
        print("  ❌ pipeline_config.py — check failed")

    # Check tissue_expression
    te = TE_PATH.read_text()
    if "POTENCY_MAP_INLINE" in te:
        print("  ✅ tissue_expression.py — inline potency map present")
        count = te.count("potency_modifier")
        print(f"  ✅ tissue_expression.py — {count} potency_modifier reference(s) found")
    else:
        print("  ❌ tissue_expression.py — POTENCY_MAP_INLINE not found")


def main() -> None:
    print("=" * 60)
    print("ATRT Potency Fix — direct injection")
    print("=" * 60)

    print("\n[1] Fixing pipeline_config.py ...")
    fix_pipeline_config()

    print("\n[2] Fixing tissue_expression.py ...")
    ok = fix_tissue_expression()

    verify()

    print("\n" + "=" * 60)
    if ok:
        print("✅ Done. Run:")
        print("   python -m backend.pipeline.save_results --disease atrt --top_n 20")
        print("\nExpected: Valproic acid score ~0.50, NOT 0.915")
        print("Expected rank 1–6: tazemetostat, alisertib, panobinostat,")
        print("                   birabresib, abemaciclib, marizomib")
    else:
        print("❌ Fix failed — see errors above.")
    print("=" * 60)


if __name__ == "__main__":
    main()