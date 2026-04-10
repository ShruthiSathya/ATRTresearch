"""
generate_figures.py
===================
ATRT Drug Repurposing Pipeline — Publication Figures

Generates 4 figures from results/atrt_pipeline_results.json:
  fig1_smarcb1_biology.png  — SMARCB1 loss + EZH2 synthetic lethality diagram
  fig2_drug_rankings.png    — Top 12 drugs with stacked score components
  fig3_score_scatter.png    — DepMap vs tissue, coloured by BBB
  fig4_confidence.png       — Top hypothesis confidence breakdown

Run:
  python -m backend.pipeline.generate_figures

Requires:
  results/atrt_pipeline_results.json  (run save_results.py first)
"""

import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

try:
    from .pipeline_config import COMPOSITE_WEIGHTS, CONFIDENCE_WEIGHTS
except ImportError:
    from pipeline_config import COMPOSITE_WEIGHTS, CONFIDENCE_WEIGHTS

RESULTS_FILE = Path("results/atrt_pipeline_results.json")
FIGURES_DIR  = Path("figures/atrt")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Colour palette
ATRT_BLUE    = "#1565C0"
ATRT_GREEN   = "#2E7D32"
ATRT_GOLD    = "#F9A825"
DARK_GREY    = "#2C2C2C"
LIGHT_GREY   = "#E0E0E0"
RED_ACCENT   = "#C5383B"
PURPLE       = "#6A1B9A"
TEAL         = "#00695C"

BBB_COLORS = {
    "HIGH":     ATRT_GREEN,
    "MODERATE": ATRT_GOLD,
    "LOW":      RED_ACCENT,
    "UNKNOWN":  "grey",
}

plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         11,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
})


def load() -> dict:
    if not RESULTS_FILE.exists():
        print(f"ERROR: {RESULTS_FILE} not found.")
        print("Run: python -m backend.pipeline.save_results first.")
        sys.exit(1)
    with open(RESULTS_FILE) as f:
        data = json.load(f)
    s = data.get("stats", {})
    print(f"Loaded {RESULTS_FILE}")
    print(f"  timestamp      : {data.get('run_timestamp','?')}")
    print(f"  drugs screened : {s.get('n_drugs_screened','?')}")
    print(f"  SMARCB1 loss   : {s.get('smarcb1_loss_count', 0)}/{s.get('total_atrt_samples', 0)}")
    print(f"  candidates     : {len(data.get('top_candidates', []))}")
    return data


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — SMARCB1 Biology + Subgroup Context
# ─────────────────────────────────────────────────────────────────────────────

def fig1_smarcb1_biology(data: dict) -> None:
    """
    Figure 1: SMARCB1 loss statistics and EZH2 synthetic lethality context.
    Shows ATRT cohort breakdown and drug boost summary.
    """
    s       = data.get("stats", {})
    smarcb1 = s.get("smarcb1_loss_count", 0)
    smarca4 = s.get("smarca4_loss_count", 0)
    total   = s.get("total_atrt_samples", 0)
    n_screen = s.get("n_drugs_screened", 0)
    n_ezh2  = s.get("n_ezh2_boosted", 0)
    n_aurka = s.get("n_aurka_boosted", 0)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.suptitle(
        "ATRT Molecular Context: SMARCB1 Loss → EZH2 Synthetic Lethality\n"
        f"n={total} ATRT samples | {n_screen} drugs screened",
        fontsize=13, fontweight="bold", y=1.02,
    )

    # Panel A: SMARCB1/SMARCA4 loss breakdown
    ax = axes[0]
    categories = ["SMARCB1\nloss (~95%)", "SMARCA4\nloss (~5%)", "Other"]
    counts     = [smarcb1, smarca4, max(0, total - smarcb1 - smarca4)]
    colors     = [RED_ACCENT, ATRT_GOLD, LIGHT_GREY]
    bars = ax.bar(categories, counts, color=colors,
                  edgecolor=DARK_GREY, linewidth=0.8, width=0.55)

    ymax = max(counts) if counts else 1
    for bar, cnt in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                cnt + ymax * 0.02, str(cnt),
                ha="center", va="bottom", fontweight="bold", fontsize=12)

    ax.set_ylabel("Samples", fontsize=12)
    ax.set_title(
        "SWI/SNF Subunit Loss in ATRT Cohort\n"
        "(Hasselblatt 2011; CBTN data)",
        fontsize=11, fontweight="bold"
    )

    ax.text(0.5, -0.22,
            "SMARCB1 loss = biallelic deletion/mutation (defines disease)\n"
            "SMARCA4 loss = alternate subunit (~5% of ATRT)",
            ha="center", transform=ax.transAxes,
            fontsize=9, color=DARK_GREY, style="italic")

    # Panel B: Drug boost summary + EZH2 rationale
    ax2 = axes[1]

    boost_labels = [
        "EZH2 inhibitors\nboosted (×1.40)",
        "AURKA inhibitors\nboosted (×1.15–1.30)",
        "SMO/GLI inhibitors\nboosted (SHH subgroup)",
        "Other candidates\n(no special boost)",
    ]
    cands  = data.get("top_candidates", [])
    n_smo  = sum(1 for c in cands
                 if any("SMO" in b for b in c.get("atrt_boosts_applied", [])))
    other  = max(0, len(cands) - n_ezh2 - n_aurka - n_smo)
    boost_counts = [n_ezh2, n_aurka, n_smo, other]
    boost_colors = [RED_ACCENT, PURPLE, TEAL, LIGHT_GREY]

    bars2 = ax2.barh(range(4), boost_counts, color=boost_colors,
                     edgecolor=DARK_GREY, linewidth=0.6, height=0.5)
    ax2.set_yticks(range(4))
    ax2.set_yticklabels(boost_labels, fontsize=10)
    ax2.set_xlabel("Number of candidates in top screening results", fontsize=10)
    ax2.set_title(
        "ATRT-Specific Drug Score Adjustments\n"
        "(Based on SMARCB1-loss biology)",
        fontsize=11, fontweight="bold"
    )

    for bar, cnt in zip(bars2, boost_counts):
        if cnt > 0:
            ax2.text(cnt + 0.1, bar.get_y() + bar.get_height() / 2,
                     str(cnt), va="center", fontweight="bold", fontsize=11)

    # Add EZH2 rationale text box
    ax2.text(0.02, -0.18,
             "EZH2 boost: SMARCB1 normally opposes PRC2. When SMARCB1 is lost,\n"
             "EZH2 becomes hyperactive and ESSENTIAL → synthetic lethality\n"
             "(Knutson 2013 PNAS; FDA Breakthrough Therapy for tazemetostat)",
             ha="left", transform=ax2.transAxes,
             fontsize=8.5, color=DARK_GREY, style="italic",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF3E0",
                       edgecolor=RED_ACCENT, alpha=0.85))

    plt.tight_layout()
    out = FIGURES_DIR / "fig1_smarcb1_biology.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  fig1 ✅ {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — Top Drug Rankings
# ─────────────────────────────────────────────────────────────────────────────

def fig2_drug_rankings(data: dict, top_n: int = 12) -> None:
    candidates = data.get("top_candidates", [])
    if not candidates:
        print("  fig2 skipped — no candidates")
        return

    top = candidates[:top_n]
    drugs  = [c["name"] for c in top]
    comps  = {
        f"Tissue/GSE70678 ({COMPOSITE_WEIGHTS['tissue']*100:.0f}%)": [
            c.get("tissue_expression_score", 0) * COMPOSITE_WEIGHTS["tissue"]
            for c in top
        ],
        f"DepMap CRISPR ({COMPOSITE_WEIGHTS['depmap']*100:.0f}%)": [
            c.get("depmap_score", 0) * COMPOSITE_WEIGHTS["depmap"]
            for c in top
        ],
        f"Escape bypass ({COMPOSITE_WEIGHTS['escape']*100:.0f}%)": [
            c.get("escape_bypass_score", 0) * COMPOSITE_WEIGHTS["escape"]
            for c in top
        ],
        f"PPI network ({COMPOSITE_WEIGHTS['ppi']*100:.0f}%)": [
            c.get("ppi_score", 0) * COMPOSITE_WEIGHTS["ppi"]
            for c in top
        ],
    }
    composite  = [c.get("score", 0) for c in top]
    bbb_status = [c.get("bbb_penetrance", "UNKNOWN") for c in top]
    ezh2       = [c.get("ezh2_boosted", False) for c in top]
    aurka      = [c.get("aurka_boosted", False) for c in top]
    failed     = [c.get("clinical_failure", False) for c in top]

    s        = data.get("stats", {})
    n_screen = s.get("n_drugs_screened", len(candidates))

    fig, ax = plt.subplots(figsize=(max(12, len(drugs) * 1.1), 6.5))
    x      = np.arange(len(drugs))
    width  = 0.18
    offsets= [-1.5, -0.5, 0.5, 1.5]
    colors = [ATRT_BLUE, ATRT_GOLD, ATRT_GREEN, PURPLE]

    for (label, vals), offset, color in zip(comps.items(), offsets, colors):
        ax.bar(x + offset * width, vals, width,
               label=label, color=color, alpha=0.85,
               edgecolor="white", linewidth=0.4)

    ax2 = ax.twinx()
    ax2.plot(x, composite, "o-", color=RED_ACCENT,
             linewidth=2.5, markersize=8, label="Composite score", zorder=5)
    ax2.set_ylim(0, 1.15)
    ax2.set_ylabel("Composite score", color=RED_ACCENT, fontsize=12)
    ax2.tick_params(axis="y", colors=RED_ACCENT)
    ax2.spines["right"].set_color(RED_ACCENT)
    ax2.spines["top"].set_visible(False)

    for i in range(len(drugs)):
        markers = []
        if ezh2[i]:
            markers.append(("E", RED_ACCENT))
        if aurka[i]:
            markers.append(("A", PURPLE))
        if failed[i]:
            markers.append(("X", "grey"))
        if bbb_status[i] == "HIGH":
            markers.append(("★", ATRT_GREEN))
        for j, (m, mc) in enumerate(markers):
            ax2.text(i, composite[i] + 0.04 + j * 0.07, m,
                     ha="center", fontsize=9, color=mc, fontweight="bold")

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    extra  = [
        mpatches.Patch(color=RED_ACCENT,   label="E EZH2 inhibitor (boosted ×1.40)"),
        mpatches.Patch(color=PURPLE,       label="A AURKA inhibitor (boosted ×1.15–1.30)"),
        mpatches.Patch(color=ATRT_GREEN,   label="★ HIGH BBB penetrance"),
        mpatches.Patch(color="grey",       label="X Known GBM clinical failure"),
    ]
    ax.legend(h1 + h2 + extra, l1 + l2 + [p.get_label() for p in extra],
              loc="upper right", fontsize=8, framealpha=0.9, ncol=2)

    ax.set_xticks(x)
    ax.set_xticklabels(drugs, rotation=30, ha="right", fontsize=11)
    ax.set_ylabel("Weighted component score", fontsize=12)
    ax.set_title(
        f"Top {len(drugs)} ATRT Drug Candidates — Multi-Omic Composite Scoring\n"
        f"{n_screen} drugs screened | GSE70678 tissue + DepMap rhabdoid CRISPR",
        fontsize=12, fontweight="bold",
    )

    max_stacked = max(
        sum(comps[k][i] for k in comps) for i in range(len(drugs))
    ) if drugs else 0.6
    ax.set_ylim(0, max_stacked * 1.35)

    plt.tight_layout()
    out = FIGURES_DIR / "fig2_drug_rankings.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  fig2 ✅ {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — Score Scatter
# ─────────────────────────────────────────────────────────────────────────────

def fig3_score_scatter(data: dict) -> None:
    candidates = data.get("top_candidates", [])
    if not candidates:
        print("  fig3 skipped — no candidates")
        return

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    fig.suptitle(
        "ATRT Score Component Relationships",
        fontsize=13, fontweight="bold"
    )

    pairs = [
        ("depmap_score", "tissue_expression_score", "DepMap CRISPR (ATRT lines)",
         "GSE70678 Expression"),
        ("ppi_score", "escape_bypass_score", "PPI Network", "Escape Bypass"),
        ("depmap_score", "score", "DepMap CRISPR", "Composite Score"),
    ]

    for ax, (xkey, ykey, xlabel, ylabel) in zip(axes, pairs):
        xs     = [c.get(xkey, 0) for c in candidates]
        ys     = [c.get(ykey, 0) for c in candidates]
        colors = [BBB_COLORS.get(c.get("bbb_penetrance", "UNKNOWN"), "grey")
                  for c in candidates]
        names  = [c.get("name", "?") for c in candidates]

        # Mark EZH2-boosted drugs differently
        for i, c in enumerate(candidates):
            if c.get("ezh2_boosted"):
                ax.scatter([xs[i]], [ys[i]], c=RED_ACCENT, s=120, alpha=0.9,
                           edgecolors=DARK_GREY, linewidths=0.8,
                           marker="*", zorder=5)
            else:
                ax.scatter([xs[i]], [ys[i]], c=[colors[i]], s=65, alpha=0.75,
                           edgecolors=DARK_GREY, linewidths=0.4, zorder=3)

        top5 = sorted(range(len(candidates)),
                      key=lambda i: candidates[i].get("score", 0),
                      reverse=True)[:5]
        for i in top5:
            ax.annotate(names[i], (xs[i], ys[i]),
                        textcoords="offset points", xytext=(5, 3),
                        fontsize=8, color=DARK_GREY, fontweight="bold")

        if len(xs) > 3:
            try:
                z  = np.polyfit(xs, ys, 1)
                xr = np.linspace(min(xs), max(xs), 100)
                ax.plot(xr, np.poly1d(z)(xr), "--",
                        color="grey", alpha=0.45, linewidth=1.2)
            except Exception:
                pass

        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_xlim(-0.05, 1.1)
        ax.set_ylim(-0.05, 1.1)

    legend_elems = [
        mpatches.Patch(facecolor=ATRT_GREEN, label="HIGH BBB"),
        mpatches.Patch(facecolor=ATRT_GOLD,  label="MODERATE BBB"),
        mpatches.Patch(facecolor=RED_ACCENT,  label="LOW BBB"),
        mpatches.Patch(facecolor="grey",      label="UNKNOWN BBB"),
        plt.scatter([], [], marker="*", c=RED_ACCENT, s=100,
                    label="EZH2 inhibitor (boosted ×1.40)"),
    ]
    axes[2].legend(handles=legend_elems, loc="lower right",
                   fontsize=9, framealpha=0.9)

    plt.tight_layout()
    out = FIGURES_DIR / "fig3_score_scatter.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  fig3 ✅ {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 — Top Hypothesis Confidence Breakdown
# ─────────────────────────────────────────────────────────────────────────────

def fig4_confidence(data: dict) -> None:
    hyps   = data.get("hypotheses", [])
    cb_top = data.get("confidence_breakdown", {})
    cb     = hyps[0].get("confidence_breakdown", cb_top) if hyps else cb_top
    if not cb:
        print("  fig4 skipped — no confidence_breakdown")
        return

    dep_raw  = cb.get("depmap_essentiality", 0)
    bbb_raw  = cb.get("bbb_penetrance", 0)
    div_raw  = cb.get("mechanistic_diversity", 0)
    conf_adj = cb.get("confidence_adjusted", cb.get("confidence", 0))
    conf_raw = cb.get("confidence_raw", conf_adj)
    conf_opt = cb.get("confidence_adjusted_optimistic", conf_adj)
    tox_mult = cb.get("toxicity_multiplier", 1.0)
    tox_flag = cb.get("toxicity_flag", "")
    hyp0     = hyps[0] if hyps else {}
    combo    = cb_top.get("drug_combo") or hyp0.get("drug_or_combo", "Top hypothesis")

    dep_w = dep_raw * CONFIDENCE_WEIGHTS["depmap"]
    bbb_w = bbb_raw * CONFIDENCE_WEIGHTS["bbb"]
    div_w = div_raw * CONFIDENCE_WEIGHTS["diversity"]
    total = dep_w + bbb_w + div_w

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    fig.suptitle(
        f"{combo}\nAdjusted confidence range: [{conf_adj:.2f}, {conf_opt:.2f}]  "
        f"(raw {conf_raw:.2f} × tox {tox_mult:.2f} — {tox_flag})",
        fontsize=12, fontweight="bold", y=1.02,
    )

    # Panel A: Confidence components
    ax = axes[0]
    comp_labels = [
        f"DepMap Essentiality\n(ATRT lines, w={CONFIDENCE_WEIGHTS['depmap']})",
        f"BBB Penetrance\n(curated PK, w={CONFIDENCE_WEIGHTS['bbb']})",
        f"Target Diversity\n(Jaccard, w={CONFIDENCE_WEIGHTS['diversity']})",
    ]
    weighted   = [dep_w, bbb_w, div_w]
    bar_colors = [ATRT_BLUE, ATRT_GOLD, ATRT_GREEN]
    sources    = ["Broad CRISPR\n(external)", "PK literature\n(external)",
                  "Jaccard overlap\n(computed)"]

    bars = ax.bar(comp_labels, weighted, color=bar_colors,
                  edgecolor=DARK_GREY, linewidth=0.8, width=0.5)

    ymax = max(weighted) if weighted else 0.5
    for bar, raw, w, src, col in zip(bars, [dep_raw, bbb_raw, div_raw],
                                     weighted, sources, bar_colors):
        ax.text(bar.get_x() + bar.get_width() / 2,
                w + ymax * 0.04, f"raw={raw:.2f} → {w:.3f}",
                ha="center", va="bottom", fontsize=9.5, fontweight="bold")
        if w > ymax * 0.1:
            tc = "white" if col != ATRT_GOLD else DARK_GREY
            ax.text(bar.get_x() + bar.get_width() / 2,
                    w / 2, src, ha="center", va="center",
                    fontsize=8, color=tc, fontweight="bold",
                    multialignment="center")

    ax.axhline(total, color=RED_ACCENT, linestyle="--", linewidth=2, alpha=0.8)
    ax.text(len(comp_labels) - 0.6, total + ymax * 0.04,
            f"Total = {total:.2f}",
            color=RED_ACCENT, fontsize=11, fontweight="bold")
    ax.set_ylabel("Weighted contribution to confidence", fontsize=11)
    ax.set_title(
        "Confidence Components\n(externally grounded — not self-referential)",
        fontsize=11, fontweight="bold"
    )
    ax.set_ylim(0, ymax * 1.8)

    # Panel B: Top-8 score composition
    ax2   = axes[1]
    top8  = data.get("top_candidates", [])[:8]
    if top8:
        names = [c["name"] for c in top8]
        dep   = [c.get("depmap_score", 0) for c in top8]
        tis   = [c.get("tissue_expression_score", 0) for c in top8]
        esc   = [c.get("escape_bypass_score", 0) for c in top8]
        ppi   = [c.get("ppi_score", 0) for c in top8]
        ezh2b = [c.get("ezh2_boosted", False) for c in top8]

        y = np.arange(len(names))
        h = 0.55
        ax2.barh(y, [d * COMPOSITE_WEIGHTS["depmap"] for d in dep], h,
                 color=ATRT_GOLD, label=f'DepMap ({COMPOSITE_WEIGHTS["depmap"]*100:.0f}%)')
        ax2.barh(y, [t * COMPOSITE_WEIGHTS["tissue"] for t in tis], h,
                 left=[d * COMPOSITE_WEIGHTS["depmap"] for d in dep],
                 color=ATRT_BLUE,
                 label=f'Tissue/GSE70678 ({COMPOSITE_WEIGHTS["tissue"]*100:.0f}%)')
        left2 = [d * COMPOSITE_WEIGHTS["depmap"] + t * COMPOSITE_WEIGHTS["tissue"]
                 for d, t in zip(dep, tis)]
        ax2.barh(y, [e * COMPOSITE_WEIGHTS["escape"] for e in esc], h,
                 left=left2, color=ATRT_GREEN, label="Escape (20%)")
        left3 = [l + e * COMPOSITE_WEIGHTS["escape"] for l, e in zip(left2, esc)]
        ax2.barh(y, [p * COMPOSITE_WEIGHTS["ppi"] for p in ppi], h,
                 left=left3, color=PURPLE, label="PPI (5%)")

        # Mark EZH2-boosted bars
        for i, boosted in enumerate(ezh2b):
            if boosted:
                ax2.text(0.02, i, "★EZH2↑", va="center",
                         fontsize=7.5, color=RED_ACCENT, fontweight="bold")

        ax2.set_yticks(y)
        ax2.set_yticklabels(names, fontsize=10)
        ax2.set_xlabel("Composite score", fontsize=11)
        ax2.set_xlim(0, 1.05)
        ax2.set_title("Score Composition — Top 8 Candidates\n(★ = EZH2 inhibitor, boosted ×1.40)",
                      fontsize=11, fontweight="bold")
        ax2.legend(loc="lower right", fontsize=9, framealpha=0.9)

    plt.tight_layout()
    out = FIGURES_DIR / "fig4_confidence.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  fig4 ✅ {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 55)
    print("generate_figures.py — ATRT pipeline results")
    print("=" * 55)

    results = load()
    print()

    fig1_smarcb1_biology(results)
    fig2_drug_rankings(results, top_n=12)
    fig3_score_scatter(results)
    fig4_confidence(results)

    print("\n" + "=" * 55)
    print(f"Figures saved to: {FIGURES_DIR}/")
    print("Add to README.md:")
    for name in ["fig1_smarcb1_biology", "fig2_drug_rankings",
                 "fig3_score_scatter", "fig4_confidence"]:
        print(f"  ![{name}](figures/atrt/{name}.png)")