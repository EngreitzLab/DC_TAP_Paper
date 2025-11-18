import marimo

__generated_with = "0.16.5"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Select Guides Notebook V1

    This notebook (ver.1) explores the metadata from CRISPick outputs to appropriately select guides. 
    Version 1 relies on the CRISPick 'Pick Order' ranks as the main metric to select guides and uses 
    accompanying CRISPick output columns to filter out less promising guides.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Dependencies""")
    return


@app.cell
def _():
    from pathlib import Path

    import math
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

    from pybedtools import BedTool
    return (
        BedTool,
        BoundaryNorm,
        Line2D,
        ListedColormap,
        MultipleLocator,
        Path,
        np,
        pd,
        plt,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Main Notebook""")
    return


@app.cell
def _(Path):
    # Define Inputs
    out_dir = Path("results/2025-10-15/")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Define Inputs (WTC11)""")
    return


@app.cell
def _(BedTool, Path, pd):
    w_metadata = pd.read_csv("metadata/guide-metadata-wtc11-intermediate.tsv", sep="\t")
    w_gcutpos = BedTool("results/2025-10-05/dctapseq-validation-wtc11-sgrna-cutpositions.bed")

    out_wm_targetstart_outward_bed = Path("results/2025-10-15/wtc11-metadata-target-start-1bp-outward.bed")
    out_wm_targetend_outward_bed = Path("results/2025-10-15/wtc11-metadata-target-end-1bp-outward.bed")

    out_wm_targetstart_bed = Path("results/2025-10-15/wtc11-metadata-target-start-1bp.bed")
    out_wm_targetend_bed = Path("results/2025-10-15/wtc11-metadata-target-end-1bp.bed")

    out_jitterplot_50 = Path("results/2025-10-15/jitterplots-50/")
    out_jitterplot_100 = Path("results/2025-10-15/jitterplots-100/")
    out_jitterplot_250 = Path("results/2025-10-15/jitterplots-250/")
    return (
        out_jitterplot_100,
        out_jitterplot_250,
        out_jitterplot_50,
        out_wm_targetend_bed,
        out_wm_targetstart_bed,
        w_gcutpos,
        w_metadata,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Create Bed Files""")
    return


@app.cell
def _(BedTool, out_wm_targetend_bed, out_wm_targetstart_bed, pd, w_metadata):
    # Create 1bp reference bed files for each target ends
    def get_targetstart_1bp_bed(df, out):
        df = df.copy().drop_duplicates(subset=["intended_target_name"]).reset_index()
        bed_df = pd.DataFrame({
            "chr": df.loc[:, "intended_target_chr"],
            "start": df.loc[:, "intended_target_start"],
            "end": df.loc[:, "intended_target_start"] + 1, # 0-based
            "name": df.loc[:, "intended_target_name"] + "_start",
            "score": ".",
            "strand": df.loc[:, "strand"]
        })
        BedTool.from_dataframe(bed_df).saveas(out)
    
    def get_targetend_1bp_bed(df, out):
        df = df.copy().drop_duplicates(subset=["intended_target_name"]).reset_index()
        bed_df = pd.DataFrame({
            "chr": df.loc[:, "intended_target_chr"],
            "start": df.loc[:, "intended_target_end"] - 1, # 0-based
            "end": df.loc[:, "intended_target_end"],
            "name": df.loc[:, "intended_target_name"] + "_end",
            "score": ".",
            "strand": df.loc[:, "strand"]
        })
        BedTool.from_dataframe(bed_df).saveas(out)

    get_targetstart_1bp_bed(w_metadata, out_wm_targetstart_bed)
    get_targetend_1bp_bed(w_metadata, out_wm_targetend_bed)
    return


@app.cell
def _(
    BedTool,
    out_jitterplot_50,
    out_wm_targetend_bed,
    out_wm_targetstart_bed,
    pd,
    w_gcutpos,
):
    # Get list of guides (via cut position) that fall within 50bp (variable) of target ends
    def get_nearby_guides(guides_bed, target_start_bed, target_end_bed, offset, out=None):
        s = BedTool(target_start_bed)
        e = BedTool(target_end_bed)

        s_w = s.slop(genome="hg38", l=offset, r=offset)
        e_w = e.slop(genome="hg38", l=offset, r=offset)

        s_nearby_guides = guides_bed.intersect(s_w, wa=True, f=0.5).to_dataframe()
        e_nearby_guides = guides_bed.intersect(e_w, wa=True, f=0.5).to_dataframe()

        if out is not None:
            window = (
                BedTool
                    .from_dataframe(pd.concat([s_w.to_dataframe(), e_w.to_dataframe()]))
                    .saveas(out / f"wtc11-metadata-target-combined-1bp+{offset}.bed")
            )
            s_nearby_guides = guides_bed.intersect(s_w, wa=True, f=0.5).to_dataframe()
            e_nearby_guides = guides_bed.intersect(e_w, wa=True, f=0.5).to_dataframe()
            guides = (
                BedTool
                    .from_dataframe(pd.concat([s_nearby_guides, e_nearby_guides]))
                    .saveas(out / f"wtc11-metadata-target-combined-1bp+{offset}-guides.bed")
            )
    
        return s_nearby_guides.loc[:,"name"].to_list(), e_nearby_guides.loc[:,"name"].to_list()

    nearby_target_start, nearby_target_end = get_nearby_guides(
        w_gcutpos, 
        out_wm_targetstart_bed, 
        out_wm_targetend_bed, 
        offset=50,
        out=out_jitterplot_50
    )
    return get_nearby_guides, nearby_target_end, nearby_target_start


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Exploration""")
    return


@app.function
# Note, duplicate function. TODO: pull out to "src/python".
def get_per_target_subset(df, intended_target):
    subset = (
        df.copy()
            .loc[df.loc[:,"intended_target_name"] == intended_target]
            .reset_index()
    )
    return subset


@app.cell
def _(
    Line2D,
    MultipleLocator,
    nearby_target_end,
    nearby_target_start,
    np,
    plt,
    w_metadata,
):
    def plot_crispickrank_vs_reltargetend(axs, x, y, offset, ma, mi, color, x_title, y_title, round_threshold):
        c = np.where(color.to_numpy() >= round_threshold, "r", "b")
        axs.scatter(x, y, c=c)

        axs.set_ylabel(y_title)
        axs.set_ylim(65, -5)
        axs.yaxis.set_major_locator(MultipleLocator(10))
        axs.yaxis.set_minor_locator(MultipleLocator(5))
    
        axs.set_xlabel(x_title)
        axs.set_xlim(-offset - 5, offset + 5)
        axs.xaxis.set_major_locator(MultipleLocator(ma))
        axs.xaxis.set_minor_locator(MultipleLocator(mi))

        axs.axvline(x=0, color='#dddddd', linestyle='--', linewidth=1.5)
        return axs
    

    # For a one target, create a jitter plot for rank vs start & end
    def plot_jitter_for_target(subset, start, end, offset=50, ma=10, mi=5, round_threshold=5):
        df_start = subset.loc[subset.loc[:, "ranked_target_name"].isin(start)]
        df_end = subset.loc[subset.loc[:, "ranked_target_name"].isin(end)]
    
        # Plotting
        fig, axs = plt.subplots(1, 2, figsize=(9, 6), constrained_layout=True)
    
        # Y-Axis
        y_start = df_start.loc[:, "crispick_pickorder"]
        x_start = df_start.loc[:, "predicted_cutposition"] - df_start.loc[:, "intended_target_start"]
        r_start = df_start.loc[:, "crispick_pickround"]
    
        y_end = df_end.loc[:, "crispick_pickorder"]
        x_end = df_end.loc[:, "predicted_cutposition"] - df_end.loc[:, "intended_target_end"]
        r_end = df_end.loc[:, "crispick_pickround"]

        a_start = plot_crispickrank_vs_reltargetend(
            axs[0], 
            x_start,
            y_start,
            offset,
            ma,
            mi,
            r_start,
            x_title="Start",
            y_title="CRISPick Pick Order",
            round_threshold=round_threshold
        )
        a_end = plot_crispickrank_vs_reltargetend(
            axs[1], 
            x_end, 
            y_end,
            offset,
            ma,
            mi,
            r_end,
            x_title="End",
            y_title="",
            round_threshold=round_threshold
        )
        a_end.tick_params(axis="y", which="both", left=False, labelleft=False)

        legend_elements = [
            Line2D([0], [0], marker="o", color="none", markerfacecolor="blue", markersize=8, label=f"{round_threshold - 1} or less"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="red", markersize=8, label=f"{round_threshold} or more"),]
        a_end.legend(handles=legend_elements, title="CRISPick Rounds", bbox_to_anchor=(1, 1), loc="upper left")

        return fig


    plot_w_jitter_target1 = plot_jitter_for_target(
        get_per_target_subset(w_metadata, "chr17:7632480-7632781"), 
        nearby_target_start, 
        nearby_target_end
    )
    plot_w_jitter_target1.suptitle("chr17:7632480-7632781", x=0.3, ha="right")
    plot_w_jitter_target1.savefig("results/2025-10-15/wtc11-metadata-target-1-jitterplot.png")
    plot_w_jitter_target1
    return (plot_jitter_for_target,)


@app.cell
def _(
    nearby_target_end,
    nearby_target_start,
    out_jitterplot_50,
    plot_jitter_for_target,
    w_metadata,
):
    def plot_w_jitter_all(df, start, end, out, round_threshold=5, offset=50, ma=10, mi=5):
        df = df.copy()
        targets = list(df.loc[:,"intended_target_name"].unique())

        inc = 0
        for t in targets:
            p = plot_jitter_for_target(
                get_per_target_subset(w_metadata, t), 
                start, 
                end,
                offset=offset,
                ma=ma,
                mi=mi
            )
            p.suptitle(f"WTC11-Target-{inc}:{t}", x=0.45, ha="right")
            p.savefig(out / f"wm-jitterplot-{inc}-{t}.png")
            inc += 1

    # 50 bp window
    plot_w_jitter_all(w_metadata, nearby_target_start, nearby_target_end, out_jitterplot_50)
    return (plot_w_jitter_all,)


@app.cell
def _(
    get_nearby_guides,
    out_jitterplot_100,
    out_jitterplot_250,
    out_wm_targetend_bed,
    out_wm_targetstart_bed,
    plot_w_jitter_all,
    w_gcutpos,
    w_metadata,
):
    # 100 bp window
    nearby_target_start_100, nearby_target_end_100 = get_nearby_guides(
        w_gcutpos, 
        out_wm_targetstart_bed, 
        out_wm_targetend_bed, 
        offset=100,
        out=out_jitterplot_100
    )

    plot_w_jitter_all(
        w_metadata, 
        nearby_target_start_100, 
        nearby_target_end_100, 
        out_jitterplot_100,
        offset=100,
        ma=20,
        mi=10
    )

    # 250 bp window
    nearby_target_start_250, nearby_target_end_250 = get_nearby_guides(
        w_gcutpos, 
        out_wm_targetstart_bed, 
        out_wm_targetend_bed, 
        offset=250,
        out=out_jitterplot_250
    )

    plot_w_jitter_all(
        w_metadata, 
        nearby_target_start_250, 
        nearby_target_end_250, 
        out_jitterplot_250,
        offset=250,
        ma=50,
        mi=25
    )
    return nearby_target_end_100, nearby_target_start_100


@app.cell
def _(nearby_target_end, nearby_target_start, np, plt, w_metadata):
    def plot_scatter_onoff_for_target(subset, start, end, round_threshold=5):
        df_start = subset.loc[subset.loc[:, "ranked_target_name"].isin(start)]
        df_end = subset.loc[subset.loc[:, "ranked_target_name"].isin(end)]

        # Fit CFD100 hits to [0, 1]
        y_s = np.asarray(df_start.loc[:, "offtarget_cfd100_hits"], dtype=float)
        y_s_min, y_s_max = np.nanmin(y_s), np.nanmax(y_s)
        y_s_norm = (y_s - y_s_min) / (y_s_max - y_s_min) if y_s_max > y_s_min else np.zeros_like(y_s)
    
        y_e = np.asarray(df_end.loc[:, "offtarget_cfd100_hits"], dtype=float)
        y_e_min, y_e_max = np.nanmin(y_e), np.nanmax(y_e)
        y_e_norm = (y_e - y_e_min) / (y_e_max - y_e_min) if y_e_max > y_e_min else np.zeros_like(y_e)

        # Plotting   
        fig, axs = plt.subplots(figsize=(9, 6), constrained_layout=True)
    
        # Data
        x_s = df_start.loc[:, "ontarget_eff_score"]
        y_s = y_s_norm
    
        x_e = df_end.loc[:, "ontarget_eff_score"]
        y_e = y_e_norm
    
        axs.scatter(x_s, y_s)
        axs.scatter(x_e, y_e)

        # Labels
        axs.set_xlabel("On-Target Efficacy Score")
        axs.set_ylabel("Fitted Off-Target CFD100 Hits")
        return fig

    plot_w_scatter_rank_target1 = plot_scatter_onoff_for_target(
        get_per_target_subset(w_metadata, "chr17:7632480-7632781"), 
        nearby_target_start, 
        nearby_target_end
    )
    plot_w_scatter_rank_target1.suptitle(
        "chr17:7632480-7632781", 
        x=0.3, 
        ha="right"
    )
    plot_w_scatter_rank_target1.savefig("results/2025-10-15/wtc11-metadata-target-1-scatterrankplot.png")
    plot_w_scatter_rank_target1
    return


@app.cell
def _(w_metadata):
    w_metadata.columns
    return


@app.cell
def _(
    BoundaryNorm,
    Line2D,
    ListedColormap,
    MultipleLocator,
    np,
    pd,
    plt,
    w_metadata,
):
    def plot_rank_vs_round(df):
        df = df.copy()
        # Color category
        cats = pd.Categorical(df["intended_target_name"])
        codes = cats.codes
        labels = list(cats.categories)

        # Discrete colormap with K colors
        K = len(labels)
        cmap = ListedColormap(plt.cm.tab20.colors[:K])
        norm = BoundaryNorm(np.arange(-0.5, K + 0.5), K)

        # Plotting
        fig, axs = plt.subplots(figsize=(9, 6))
        x = df.loc[:,"crispick_pickround"]
        y = df.loc[:,"crispick_pickorder"]
        axs.scatter(x, y, c=codes, cmap=cmap, norm=norm)
        axs.axvline(x=5, color='#dddddd', linestyle='--', linewidth=1.5)

        legend_handles = [
            Line2D(
                [0], [0], 
                marker='o', 
                linestyle='', 
                markerfacecolor=cmap(i),
                markeredgecolor='none', 
                label=lab
            )
            for i, lab in enumerate(labels)
        ]
        axs.legend(
            handles=legend_handles,
            title="Target",
            bbox_to_anchor=(1, 1),
            loc="best"
        )

        axs.set_ylabel("CRISPick Pick Order")
        axs.set_ylim(65, -5)
        axs.yaxis.set_major_locator(MultipleLocator(10))
        axs.yaxis.set_minor_locator(MultipleLocator(5))
    
        axs.set_xlabel("CRISPick Picking Round")
        axs.xaxis.set_major_locator(MultipleLocator(2))
        axs.xaxis.set_minor_locator(MultipleLocator(1))

        fig.suptitle("CRISPick Pick Order vs Picking Rounds")

        return fig
        

    plot_w_rankvround = plot_rank_vs_round(w_metadata)
    plot_w_rankvround.savefig("results/2025-10-15/wtc11-metadata-rankvround-plot.png", bbox_inches='tight')
    plot_w_rankvround
    return


@app.cell
def _(
    nearby_target_end_100,
    nearby_target_start_100,
    plot_jitter_for_target,
    w_metadata,
):
    plot_w_jitter_target3 = plot_jitter_for_target(
        get_per_target_subset(w_metadata, "chr17:7560401-7560702"), 
        nearby_target_start_100, 
        nearby_target_end_100,    
        offset=100,
        ma=20,
        mi=10,
        round_threshold=4
    )
    plot_w_jitter_target3.suptitle("WTC11-Target-3:chr17:7560401-7560702", x=0.44, ha="right")
    plot_w_jitter_target3.savefig("results/2025-10-15/jitterplots-100/alt/wtc11-metadata-target-3-jitterplot.png")
    plot_w_jitter_target3
    return


@app.cell
def _(BedTool, Path, pd):
    # One-off
    def convert_cutsite_to_bed(design_file, out):
        df = pd.read_csv(design_file, sep="\t")
        bed_df = pd.DataFrame({
            "chr": "chr17",
            "start": df.loc[:, "sgRNA Cut Position (1-based)"] - 1, # Convert to 0-based
            "end": df.loc[:, "sgRNA Cut Position (1-based)"],
            "name": "one-off",
            "score": ".",
            "strand": "."
        })
        BedTool.from_dataframe(bed_df).saveas(out) 

    convert_cutsite_to_bed(
        design_file=Path("results/2025-10-15/followup-workspace/repeatmasker/repeatmasker-naive-check-sgrna-designs.txt"),
        out=Path("results/2025-10-15/followup-workspace/repeatmasker/repeatmasker-naive-check-sgrna-designs.bed")
    )
    return


if __name__ == "__main__":
    app.run()
