#!/usr/bin/env python3
"""
make_figures.py

Step 5 / 6 figures for the DC-TAP-seq x ENCODE overlap pipeline.

Reads the pre-computed matrices and enrichment tables and regenerates:
  heatmaps/heatmap_{CELL}.png       positive elements x informative targets
  enrichment/volcano_neg_powered.png       log2 OR vs -log10 FDR (both cells)
  enrichment/dotplot_neg_powered.png       top-18 enriched targets w/ 95% CI (both cells)
  enrichment/dotplot_neg_powered_all.png   ALL FDR<0.1 targets w/ 95% CI (reference)

Inputs (produced by compute_overlap_enrichment.py):
  overlap/overlap_matrix_{CELL}.csv
  overlap/target_assay_{CELL}.csv
  enrichment/enrichment_{CELL}.csv
  elements/element_manifest.csv   (for significant-gene row labels)

Usage:
  python make_figures.py
"""
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
# Editable text in vector output: embed TrueType (type-42) glyphs as real
# text objects so Illustrator/Inkscape can select and edit them, rather than
# converting characters to outline shapes.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, leaves_list

warnings.filterwarnings("ignore")

CELLS = ["K562", "WTC11"]
FLAG_COLS = ["chr", "start", "end", "is_pos_neg", "is_pos_pos", "is_bg_powered", "is_bg_all"]

# ChIP targets pulled to the far-left of every heatmap so the expected
# chromatin trends (active mark / repressive mark / insulator) read at a glance.
KEY_MARKS = ["H3K27ac", "H3K27me3", "CTCF"]
# element-to-TSS proximity cutoff for the heatmap annotation strip
TSS_NEAR_BP = 100_000


def save_fig(fig, path_noext, dpi=300):
    """Write both an editable-text PDF (for Illustrator) and a PNG (for the
    HTML report) from a single figure. `path_noext` has no extension."""
    fig.savefig(f"{path_noext}.pdf", bbox_inches="tight")
    fig.savefig(f"{path_noext}.png", dpi=dpi, bbox_inches="tight")

# Manuscript element-category color scheme (Table S3 element_category)
CATEGORY_COLORS = {
    "H3K27me3 element": "#429130",
    "CTCF element":     "#49bcbc",
    "High H3K27ac":     "#c5373d",
    "H3K27ac":          "#d9694a",
    "No H3K27ac":       "#c5cad7",
}
CATEGORY_ORDER = ["H3K27me3 element", "CTCF element", "High H3K27ac", "H3K27ac", "No H3K27ac"]
# Row ordering requested by the reviewer reply: active -> insulator -> repressive
ROW_CAT_ORDER = ["High H3K27ac", "H3K27ac", "No H3K27ac", "CTCF element", "H3K27me3 element"]


def load(cell):
    full = pd.read_csv(f"overlap/overlap_matrix_{cell}.csv", index_col=0)
    for c in ["is_pos_neg", "is_pos_pos", "is_bg_powered", "is_bg_all"]:
        full[c] = full[c].astype(bool)
    mat = full.drop(columns=FLAG_COLS)
    ta = pd.read_csv(f"overlap/target_assay_{cell}.csv", index_col=0)["assay"]
    enr = pd.read_csv(f"enrichment/enrichment_{cell}.csv")
    return full, mat, ta, enr


def cluster_order(binmat, axis=0):
    X = binmat.values if axis == 0 else binmat.values.T
    if X.shape[0] < 3:
        return list(range(X.shape[0]))
    try:
        return list(leaves_list(linkage(X, method="average", metric="jaccard")))
    except Exception:
        return list(range(X.shape[0]))


def _hex2rgb(hexc):
    hexc = hexc.lstrip("#")
    return [int(hexc[j:j + 2], 16) / 255 for j in (0, 2, 4)]


def make_heatmap(cell, full, mat, ta, enr, man):
    posneg = set(full.index[full["is_pos_neg"]])
    pospos = full.index[full["is_pos_pos"]]
    pos_all = list(full.index[full["is_pos_neg"]]) + [k for k in pospos if k not in posneg]

    enr_targets = set(enr[(enr["background"] == "powered") & (enr["FDR"] < 0.1)]["target"])
    submat = mat.loc[pos_all]
    frac_targets = set(submat.columns[submat.mean(axis=0) >= 0.25])
    keep_t = sorted(t for t in (enr_targets | frac_targets) if submat[t].sum() > 0)
    # always show the three key marks so the expected trends (and the absence
    # of the repressive mark) are visible, even at zero overlap
    for km in KEY_MARKS:
        if km in submat.columns and km not in keep_t:
            keep_t.append(km)
    M = submat[keep_t]

    minfo = man[man["cell"] == cell].drop_duplicates("elem_key").set_index("elem_key")
    sg = minfo["sig_genes"]
    cat = minfo["element_category"]
    dist = minfo["min_sig_dist_tss"]        # min element-to-TSS distance over sig pairs
    neg_keys = [k for k in M.index if k in posneg]
    pos_keys = [k for k in M.index if k not in posneg]

    # --- rows: ordered BY CHROMATIN CATEGORY (not clustered), split neg/pos ---
    cat_rank = {c: i for i, c in enumerate(ROW_CAT_ORDER)}

    def order_block(keys):
        # sort by category (High H3K27ac -> H3K27ac -> No H3K27ac -> CTCF -> H3K27me3),
        # then by element coordinate for a stable within-category order
        return sorted(keys, key=lambda k: (cat_rank.get(cat.get(k), 99), k))

    neg_ord = order_block(neg_keys)
    pos_ord = order_block(pos_keys)
    row_order = neg_ord + pos_ord
    n_neg = len(neg_ord)

    # --- columns: KEY MARKS (H3K27ac, H3K27me3, CTCF) pinned far-left ---
    key_cols = [c for c in KEY_MARKS if c in keep_t]
    rest = [c for c in keep_t if c not in key_cols]
    rest_ord = [rest[i] for i in cluster_order(M[rest], axis=1)] if len(rest) >= 3 else rest
    # within the remaining columns keep histone marks together, then TFs
    rest_hist = [c for c in rest_ord if ta.get(c) == "Histone ChIP-seq"]
    rest_tf = [c for c in rest_ord if ta.get(c) != "Histone ChIP-seq"]
    cols = key_cols + rest_hist + rest_tf
    D = M.loc[row_order, cols]
    n_rows = len(row_order)

    # left panels: [category strip | <100kb strip | overlap matrix]
    main_w = max(10, len(cols) * 0.14)
    fig, (axc, axd, ax) = plt.subplots(
        1, 3, figsize=(main_w + 1.8, max(5, n_rows * 0.17)),
        gridspec_kw={"width_ratios": [0.55, 0.4, main_w], "wspace": 0.03})

    # category strip
    strip = np.zeros((n_rows, 1, 3))
    for i, k in enumerate(row_order):
        strip[i, 0] = _hex2rgb(CATEGORY_COLORS.get(cat.get(k), "#ffffff"))
    axc.imshow(strip, aspect="auto", interpolation="none")
    axc.set_xticks([]); axc.set_xlabel("category", fontsize=6)
    ylabels = [f"{sg.get(k, '?')} | {k}  ({'-' if k in posneg else '+'})" for k in row_order]
    axc.set_yticks(range(n_rows)); axc.set_yticklabels(ylabels, fontsize=5)

    # <100 kb to TSS strip: filled = near (<100kb), open = far
    dstrip = np.ones((n_rows, 1, 3))
    near_flags = []
    for i, k in enumerate(row_order):
        d = dist.get(k, np.nan)
        near = (not pd.isna(d)) and (d < TSS_NEAR_BP)
        near_flags.append(near)
        dstrip[i, 0] = _hex2rgb("#2c3e50") if near else _hex2rgb("#ecf0f1")
    axd.imshow(dstrip, aspect="auto", interpolation="none")
    axd.set_xticks([]); axd.set_yticks([])
    axd.set_xlabel("<100kb\nto TSS", fontsize=5.5)

    # main overlap matrix
    ax.imshow(D.values, aspect="auto", cmap="Greys", vmin=0, vmax=1, interpolation="none")
    ax.set_yticks([])
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, fontsize=5)
    for t, lab in zip(ax.get_xticklabels(), cols):
        if lab in key_cols:
            t.set_color("#c0392b"); t.set_fontweight("bold")
        elif ta.get(lab) == "Histone ChIP-seq":
            t.set_color("#e67e22")
        else:
            t.set_color("#333333")

    # neg/pos divider on all three panels
    for a in (axc, axd, ax):
        a.axhline(n_neg - 0.5, color="#c0392b", lw=1.5)
    # divider after the pinned key-mark columns
    ax.axvline(len(key_cols) - 0.5, color="#c0392b", lw=1.2)
    # divider between the remaining histone block and the TF block
    ax.axvline(len(key_cols) + len(rest_hist) - 0.5, color="#2980b9", lw=1.0, ls="--")

    ax.set_title(f"{cell}: DC-TAP-seq positive elements x ChIP-seq target overlap\n"
                 f"rows above red line = negative effect (-), below = positive (+); "
                 f"rows ordered by chromatin category; "
                 f"red bold columns = key marks (pinned left)", fontsize=8)
    ax.set_xlabel(f"ChIP-seq targets (n={len(cols)}: {len(key_cols)} key + "
                  f"{len(rest_hist)} other histone + {len(rest_tf)} TF)", fontsize=7)

    # legends: element category + <100kb strip
    present = [c for c in ROW_CAT_ORDER if c in set(cat.get(k) for k in row_order)]
    handles = [Line2D([0], [0], marker="s", color="w", markerfacecolor=CATEGORY_COLORS[c],
                      markersize=7, label=c) for c in present]
    handles += [
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#2c3e50", markersize=7,
               label="<100 kb to TSS"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#ecf0f1",
               markeredgecolor="#bbb", markersize=7, label=u"\u2265100 kb to TSS"),
    ]
    ax.legend(handles=handles, title="Element category / proximity", loc="upper left",
              bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=6, title_fontsize=6.5)

    fig.tight_layout()
    os.makedirs("heatmaps", exist_ok=True)
    save_fig(fig, f"heatmaps/heatmap_{cell}")
    plt.close(fig)
    print(f"heatmap_{cell}: {D.shape[0]} elements x {D.shape[1]} targets "
          f"({len(key_cols)} key marks pinned; {sum(near_flags)} rows <100kb to TSS)")


def make_volcano(data):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, cell in zip(axes, CELLS):
        enr = data[cell]
        r = enr[(enr["effect_sign"] == "neg") & (enr["background"] == "powered")].copy()
        r["neglog10fdr"] = -np.log10(r["FDR"].clip(lower=1e-300))
        hist = r["assay"] == "Histone ChIP-seq"
        ax.scatter(r.loc[~hist, "log2_OR"], r.loc[~hist, "neglog10fdr"], s=18,
                   c="#7f8c8d", alpha=0.5, label="TF", edgecolors="none")
        ax.scatter(r.loc[hist, "log2_OR"], r.loc[hist, "neglog10fdr"], s=34,
                   c="#e67e22", alpha=0.9, label="Histone", edgecolors="k", linewidths=0.3)
        ax.axhline(-np.log10(0.1), ls="--", lw=0.8, c="#c0392b")
        ax.text(ax.get_xlim()[1], -np.log10(0.1), " FDR 0.1", color="#c0392b",
                va="bottom", ha="right", fontsize=7)
        for _, row in r[r["FDR"] < 0.1].nsmallest(10, "FDR").iterrows():
            ax.annotate(row["target"], (row["log2_OR"], row["neglog10fdr"]), fontsize=6.5,
                        xytext=(3, 2), textcoords="offset points")
        ax.set_xlabel("log2 odds ratio (positives vs background)")
        ax.set_ylabel("-log10 FDR")
        ax.set_title(f"{cell}: negative-effect positives")
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
        ax.legend(frameon=False, fontsize=7, loc="upper left")
    fig.suptitle("ChIP-seq target enrichment in DC-TAP-seq negative-effect positives "
                 "(vs well-powered background)", fontsize=9.5)
    fig.tight_layout()
    save_fig(fig, "enrichment/volcano_neg_powered")
    plt.close(fig)
    print("volcano_neg_powered.png/.pdf")


def _log2or_ci(row, z=1.96):
    """95% CI for the log2 odds ratio using the Haldane-Anscombe (+0.5)
    correction, matching how log2_OR is computed upstream. Returns
    (log2_low, log2_high) as absolute log2-OR bounds."""
    a = row["pos_overlap"] + 0.5
    b = (row["pos_total"] - row["pos_overlap"]) + 0.5
    c = row["bg_overlap"] + 0.5
    d = (row["bg_total"] - row["bg_overlap"]) + 0.5
    ln_or = np.log((a * d) / (b * c))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)   # Woolf SE on ln(OR)
    lo, hi = ln_or - z * se, ln_or + z * se
    return lo / np.log(2), hi / np.log(2)          # convert to log2 scale


def _fdr_stars(fdr):
    """Significance stars from the BH-adjusted Fisher's-exact p-value (FDR)."""
    if fdr < 1e-3:
        return "***"
    if fdr < 1e-2:
        return "**"
    if fdr < 5e-2:
        return "*"
    if fdr < 1e-1:
        return "\u2020"      # dagger: 0.05 <= FDR < 0.1
    return ""


def _fmt_p(p):
    """Compact exact p-value label: 2 sig-figs, scientific below 1e-3."""
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "p=NA"
    if p < 1e-3:
        s = f"{p:.1e}"                     # e.g. 1.6e-06
        m, e = s.split("e")
        return f"p={m}\u00d710$^{{{int(e)}}}$"   # p=1.6×10^-6
    return f"p={p:.3f}"


def make_lollipop(data, top_n=18, sign="neg", background="powered",
                  outfile="enrichment/dotplot_neg_powered.png"):
    """Dot plot of per-target enrichment (FDR<0.1) with 95% CI error bars on
    the log2 odds ratio. A CI crossing 0 means the enrichment is not
    individually significant at the 95% level.

    top_n : int or None
        Show only the top_n targets by FDR (default 18). Pass None to show
        every FDR<0.1 target.
    sign : {'neg','pos'}
        Effect-size sign of the positives to test.
    background : {'powered','all'}
        Which background to compare against.
    outfile : str
        Output PNG path (.pdf written alongside).

    Row spacing is held CONSTANT across the two cell-type panels: both share a
    common y-extent set by whichever panel has more rows, so one row occupies
    the same vertical height in K562 and WTC11. If neither cell has an enriched
    target the figure is skipped and the function returns False.
    """
    sign_lab = {"neg": "neg-effect", "pos": "pos-effect"}[sign]
    bg_lab = {"powered": "well-powered background", "all": "all-non-sig background"}[background]

    # gather the rows to plot per cell, and the max row count (for shared spacing)
    tops = {}
    for cell in CELLS:
        enr = data[cell]
        r = enr[(enr["effect_sign"] == sign) & (enr["background"] == background)].copy()
        sig = r[r["FDR"] < 0.1]
        top = (sig if top_n is None else sig.nsmallest(top_n, "FDR")).sort_values("log2_OR")
        if len(top):
            ci = top.apply(_log2or_ci, axis=1, result_type="expand")
            top["ci_lo"], top["ci_hi"] = ci[0].values, ci[1].values
        tops[cell] = top
    nmax = max((len(t) for t in tops.values()), default=0)
    if nmax == 0:
        print(f"{os.path.basename(outfile)}: no targets at FDR<0.1 "
              f"({sign_lab} positives vs {bg_lab}) — figure skipped")
        return False

    fig, axes = plt.subplots(1, 2, figsize=(12, max(5.5, 0.28 * nmax + 1.5)))
    lab_fs = 7 if nmax <= 20 else (5 if nmax <= 45 else 4)
    for ax, cell in zip(axes, CELLS):
        top = tops[cell]
        ymax = len(top)
        if ymax:
            for i, (_, row) in enumerate(top.iterrows()):
                ec = "#e67e22" if row["assay"] == "Histone ChIP-seq" else "#34495e"
                ax.plot([row["ci_lo"], row["ci_hi"]], [i, i], color=ec, lw=1.2, alpha=0.55, zorder=2)
            colors = ["#e67e22" if a == "Histone ChIP-seq" else "#34495e" for a in top["assay"]]
            ax.scatter(top["log2_OR"], range(ymax), c=colors,
                       s=[30 + 180 * f for f in top["pos_frac"]], zorder=3,
                       edgecolors="k", linewidths=0.3)
            ax.set_yticks(range(ymax)); ax.set_yticklabels(top["target"], fontsize=lab_fs)
        else:
            ax.set_yticks([])
            ax.text(0.5, 0.5, "no targets at FDR<0.1", transform=ax.transAxes,
                    ha="center", va="center", fontsize=9, color="#999", style="italic")
        # SHARED y-extent -> identical row spacing across both panels
        ax.set_ylim(-0.7, nmax - 0.3)
        ax.set_xlabel("log2 odds ratio (95% CI)")
        scope = "top enriched targets" if top_n is not None else "all enriched targets"
        ax.set_title(f"{cell}: {scope}\n({sign_lab} positives, FDR<0.1; n={ymax})")
        ax.axvline(0, color="#c0392b", lw=0.8, ls="--", zorder=1)
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
        if ymax:
            xlo = min(0, top["ci_lo"].min())
            xhi = top["ci_hi"].max()
            xr = xhi - xlo if xhi > xlo else 1.0
            # extra right headroom for the overlap + exact-p + stars labels
            ax.set_xlim(xlo - 0.05 * xr, xhi + 0.75 * xr)
            for i, (_, row) in enumerate(top.iterrows()):
                stars = _fdr_stars(row["FDR"])
                ax.text(row["ci_hi"] + 0.03 * xr, i,
                        f"{row['pos_overlap']}/{row['pos_total']}  "
                        f"{_fmt_p(row['p_value'])} {stars}",
                        va="center", fontsize=5.5, color="#555")
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor="#34495e", label="TF", markersize=7),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="#e67e22", label="Histone", markersize=7)]
    fig.legend(handles=leg, loc="upper right", bbox_to_anchor=(0.995, 0.965),
               frameon=False, fontsize=8, ncol=2)
    fig.suptitle(f"Fisher's exact test, BH-FDR corrected ({sign_lab} positives vs {bg_lab}). "
                 "Point = log2 odds ratio "
                 "(size \u221d fraction of positives overlapping); bar = 95% CI (Woolf, "
                 "Haldane); label = overlap/total + exact Fisher p + FDR stars "
                 "(*** FDR<0.001, ** <0.01, * <0.05, \u2020 <0.1)", fontsize=7.5, y=0.995)
    fig.tight_layout()
    save_fig(fig, os.path.splitext(outfile)[0])
    plt.close(fig)
    print(f"{os.path.basename(outfile)} (dot plot with 95% CI + FDR stars, n={nmax} max rows)")
    return True


def main():
    man = pd.read_csv("elements/element_manifest.csv")
    per_cell = {}
    enr_data = {}
    for cell in CELLS:
        full, mat, ta, enr = load(cell)
        per_cell[cell] = (full, mat, ta, enr)
        enr_data[cell] = enr
        make_heatmap(cell, full, mat, ta, enr, man)
    make_volcano(enr_data)
    make_lollipop(enr_data, top_n=18, sign="neg", background="powered",
                  outfile="enrichment/dotplot_neg_powered.png")
    # full reference version: every FDR<0.1 target (not embedded in the report)
    make_lollipop(enr_data, top_n=None, sign="neg", background="powered",
                  outfile="enrichment/dotplot_neg_powered_all.png")
    # positive-effect positives (well-powered background): typically empty, but
    # produced for completeness; the function skips output if nothing is enriched.
    make_lollipop(enr_data, top_n=None, sign="pos", background="powered",
                  outfile="enrichment/dotplot_pos_powered_all.png")


if __name__ == "__main__":
    main()
