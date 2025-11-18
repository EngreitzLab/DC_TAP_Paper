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
    # Select Guides Notebook

    This notebook explores the metadata from CRISPick outputs to appropriately select guides.
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
    return Path, math, np, pd, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Main Notebook""")
    return


@app.cell
def _(Path):
    # Define Inputs
    out_dir = Path("results/2025-10-13/")
    return (out_dir,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Define Inputs (WTC11)""")
    return


@app.cell
def _(pd):
    w_metadata = pd.read_csv("metadata/guide-metadata-wtc11-intermediate.tsv", sep="\t")
    return (w_metadata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Define Inputs (K562)""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Exploration""")
    return


@app.function
def get_per_target_subset(df, intended_target):
    subset = (
        df.copy()
            .loc[df.loc[:,"intended_target_name"] == intended_target]
            .reset_index()
    )
    return subset


@app.cell
def _(plt):
    def plot_cp_offtarget_distribution(subset):
        subset = subset.sort_values(by="crispick_pickorder")
        fig, axs = plt.subplots()
        x = subset.loc[:,"crispick_pickorder"]
        y = subset.loc[:, "agg_cfd_score"]
        axs.plot(x, y)

        return fig, axs
    return (plot_cp_offtarget_distribution,)


@app.cell
def _(plot_cp_offtarget_distribution, w_metadata):
    w_cp_offtarget_dplot, w_cp_offtarget_dplot_a = plot_cp_offtarget_distribution(
        get_per_target_subset(w_metadata, "chr17:7632480-7632781")
    )
    w_cp_offtarget_dplot_a.set_xlabel("CRISPick Pick Order")
    w_cp_offtarget_dplot_a.set_ylabel("Agg CFD Score")
    # w_cp_offtarget_dplot.savefig(out_dir / "CPvsCFD.png")
    w_cp_offtarget_dplot
    return


@app.cell
def _(plt):
    def plot_cp_offtarget_hits(subset):
        subset = subset.sort_values(by="crispick_pickorder")
        fig, axs = plt.subplots()
        x = subset.loc[:,"crispick_pickorder"]
        y = subset.loc[:, "offtarget_cfd100_hits"]
        axs.plot(x, y)

        return fig, axs
    return (plot_cp_offtarget_hits,)


@app.cell
def _(plot_cp_offtarget_hits, w_metadata):
    w_cp_offtarget_hplot, w_cp_offtarget_hplot_a = plot_cp_offtarget_hits(
        get_per_target_subset(w_metadata, "chr17:7632480-7632781")
    )
    w_cp_offtarget_hplot_a.set_xlabel("CRISPick Pick Order")
    w_cp_offtarget_hplot_a.set_ylabel("CFD1 Off-Target Hits")
    # w_cp_offtarget_hplot.savefig(out_dir / "CPvsCFD1hits.png")
    w_cp_offtarget_hplot
    return


@app.cell
def _(plt):
    def plot_cp_ontarget_distribution(subset):
        subset = subset.sort_values(by="crispick_pickorder")
        fig, axs = plt.subplots()
        x = subset.loc[:,"crispick_pickorder"]
        y = subset.loc[:, "ontarget_eff_score"]
        axs.plot(x, y)

        return fig, axs
    return (plot_cp_ontarget_distribution,)


@app.cell
def _(plot_cp_ontarget_distribution, w_metadata):
    w_cp_ontarget_dplot, w_cp_ontarget_dplot_a = plot_cp_ontarget_distribution(
        get_per_target_subset(w_metadata, "chr17:7632480-7632781")
    )
    w_cp_ontarget_dplot_a.set_xlabel("CRISPick Pick Order")
    w_cp_ontarget_dplot_a.set_ylabel("On-Target Efficacy Score")
    # w_cp_ontarget_dplot.savefig(out_dir / "CPvsOnTarget.png")
    w_cp_ontarget_dplot
    return


@app.cell
def _(np, plt):
    def plot_offtarget_distribution(subset, crit_threshold=None, axs=None, title=None):
        query = "agg_cfd_score"
        subset = subset.sort_values(by=query, ascending=False)

        if axs is None:
            fig, axs = plt.subplots()

        x = subset.loc[:,"guide_id"].str.replace(r'^[^_]*_[^_]*_', '', regex=True)
        y = subset.loc[:, query]
        axs.plot(x, y)

        if crit_threshold is not None:
            axs.plot(x, np.full_like(x, crit_threshold))
            s = subset.loc[subset.loc[:, query] >= crit_threshold]
            ids = list(s.loc[:, "guide_id"])
            print(f"Total Guides above Threshold: {len(ids)}")
            print(f"Guide_IDs >= {crit_threshold} Agg CFD Score: {ids}")

        if title:
            axs.set_title(title)

        return axs.figure, axs
    return (plot_offtarget_distribution,)


@app.cell
def _(plot_offtarget_distribution, w_metadata):
    w_offtarget_dplot, w_offtarget_dplot_a = plot_offtarget_distribution(
        get_per_target_subset(w_metadata, "chr17:7632480-7632781"),
        crit_threshold=40
    )
    w_offtarget_dplot_a.set_xlabel("Guide ID")
    w_offtarget_dplot_a.set_ylabel("Agg. CFD Score")
    w_offtarget_dplot_a.tick_params(axis="x", rotation=45, labelsize=6)
    # w_offtarget_dplot.savefig(out_dir / "IDvsCFD.png")
    w_offtarget_dplot
    return


@app.cell
def _(np, plt):
    def plot_offtarget_hits(subset, crit_threshold=None, axs=None, title=None):
        query = "offtarget_cfd100_hits"
        subset = subset.sort_values(by=query, ascending=False)

        if axs is None:
            fig, axs = plt.subplots()

        x = subset.loc[:,"guide_id"].str.replace(r'^[^_]*_[^_]*_', '', regex=True)
        y = subset.loc[:, query]
        axs.plot(x, y)

        if crit_threshold is not None:
            axs.plot(x, np.full_like(x, crit_threshold))
            s = subset.loc[subset.loc[:, query] >= crit_threshold]
            ids = list(s.loc[:, "guide_id"])
            print(f"Total Guides above Threshold: {len(ids)}")
            print(f"Guide_IDs >= {crit_threshold} Off-Target Hits: {ids}")

        if title:
            axs.set_title(title)

        return axs.figure, axs
    return (plot_offtarget_hits,)


@app.cell
def _(plot_offtarget_hits, w_metadata):
    w_offtarget_hplot, w_offtarget_hplot_a = plot_offtarget_hits(
        get_per_target_subset(w_metadata, "chr17:7632480-7632781"),
        crit_threshold=1
    )

    w_offtarget_hplot_a.set_xlabel("Guide ID")
    w_offtarget_hplot_a.set_ylabel("CFD1 Off-Target Hits")
    w_offtarget_hplot_a.tick_params(axis="x", rotation=45, labelsize=6)
    # w_offtarget_hplot.savefig(out_dir / "IDvsCFD1hits.png")
    w_offtarget_hplot
    return


@app.cell
def _(np, plt):
    def plot_ontarget_distribution(subset, crit_threshold=None, axs=None, title=None):
        query = "ontarget_eff_score"
        subset = subset.sort_values(by=query, ascending=False)

        if axs is None:
            fig, axs = plt.subplots()

        x = subset.loc[:,"guide_id"].str.replace(r'^[^_]*_[^_]*_', '', regex=True)
        y = subset.loc[:, query]
        axs.plot(x, y)

        if crit_threshold is not None:
            axs.plot(x, np.full_like(x, crit_threshold))
            s = subset.loc[subset.loc[:, query] >= crit_threshold]
            ids = list(s.loc[:, "guide_id"])
            print(f"Total Guides above Threshold: {len(ids)}")
            print(f"Guide_IDs >= {crit_threshold} On-Target Efficacy Score: {ids}")

        if title:
            axs.set_title(title) 

        return axs.figure, axs
    return (plot_ontarget_distribution,)


@app.cell
def _(plot_ontarget_distribution, w_metadata):
    w_ontarget_dplot, w_ontarget_dplot_a = plot_ontarget_distribution(
        get_per_target_subset(w_metadata, "chr17:7632480-7632781"),
        crit_threshold=0.6
    )

    w_ontarget_dplot_a.set_xlabel("Guide ID")
    w_ontarget_dplot_a.set_ylabel("On-Target Efficacy Score")
    w_ontarget_dplot_a.tick_params(axis="x", rotation=45, labelsize=6)
    # w_ontarget_dplot.savefig(out_dir / "IDvsOnTarget.png")
    w_ontarget_dplot
    return


@app.cell
def _(math, np, plt):
    def plot_all_by_query(df, query_func, batches=5):
        df = df.copy()
        targets = list(df.loc[:,"intended_target_name"].unique())

        r = math.ceil(len(targets) / batches)
        c = batches
        fig , axs = plt.subplots(r, c, figsize=(4*c, 3*r), constrained_layout=True)

        increment = 0
        for i, j in np.ndindex(axs.shape):
            t = targets[increment]
            _, a = query_func(
                get_per_target_subset(df, t),
                axs = axs[i, j],
                title = t
            )
            increment += 1

        return fig
    return (plot_all_by_query,)


@app.cell
def _(
    plot_all_by_query,
    plot_offtarget_distribution,
    plot_offtarget_hits,
    plot_ontarget_distribution,
    w_metadata,
):
    w_plots_offtarget_dist = plot_all_by_query(w_metadata, plot_offtarget_distribution)
    w_plots_offtarget_hits = plot_all_by_query(w_metadata, plot_offtarget_hits)
    w_plots_ontarget_dist = plot_all_by_query(w_metadata, plot_ontarget_distribution)
    return (
        w_plots_offtarget_dist,
        w_plots_offtarget_hits,
        w_plots_ontarget_dist,
    )


@app.cell
def _(w_plots_offtarget_dist):
    w_plots_offtarget_dist
    return


@app.cell
def _(w_plots_offtarget_hits):
    w_plots_offtarget_hits
    # w_plots_offtarget_hits.savefig(out_dir / "IDvsCFD1hits_all.png")
    return


@app.cell
def _(w_plots_ontarget_dist):
    w_plots_ontarget_dist
    # w_plots_ontarget_dist.savefig(out_dir / "IDvsOnTarget_all.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Results Summary""")
    return


@app.cell
def _(w_metadata):
    print(w_metadata)
    return


@app.cell
def _(out_dir, plot_cp_offtarget_hits, w_metadata):
    w_cp_offtarget_hplot_chr19, w_cp_offtarget_hplot_a_chr19 = plot_cp_offtarget_hits(
        get_per_target_subset(w_metadata, "chr19:13163707-13164140")
    )
    w_cp_offtarget_hplot_a_chr19.set_xlabel("CRISPick Pick Order")
    w_cp_offtarget_hplot_a_chr19.set_ylabel("CFD1 Off-Target Hits")
    w_cp_offtarget_hplot_a_chr19.axvline(x=25, color='#dddddd', linestyle='--', linewidth=1.5)

    w_cp_offtarget_hplot_chr19.suptitle("WTC11-Target-4:chr19:13163707-13164140")
    w_cp_offtarget_hplot_chr19.savefig(out_dir / "CPvsCFD1hits-chr19.png")
    w_cp_offtarget_hplot_chr19
    return


@app.cell
def _(out_dir, plot_cp_offtarget_hits, w_metadata):
    w_cp_offtarget_hplot_ch17, w_cp_offtarget_hplot_a_ch17 = plot_cp_offtarget_hits(
        get_per_target_subset(w_metadata, "chr17:7560401-7560702")
    )
    w_cp_offtarget_hplot_a_ch17.set_xlabel("CRISPick Pick Order")
    w_cp_offtarget_hplot_a_ch17.set_ylabel("CFD1 Off-Target Hits")
    w_cp_offtarget_hplot_a_ch17.axvline(x=20, color='#dddddd', linestyle='--', linewidth=1.5)

    w_cp_offtarget_hplot_ch17.suptitle("WTC11-Target-3:chr17:7560401-7560702")
    w_cp_offtarget_hplot_ch17.savefig(out_dir / "CPvsCFD1hits-ch17.png")
    w_cp_offtarget_hplot_ch17
    return


if __name__ == "__main__":
    app.run()
