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
    # DC-TAP-seq Validation Experiment Power Analysis Notebook

    This notebook generates ... todo documentation

    Currently the notebook below is mostly scratch work.
    """
    )
    return


@app.cell
def _():
    import os
    import pandas as pd
    return (pd,)


@app.cell
def _(pd):
    raw_tbl_s3_wtc11_significant = pd.read_csv(
        "results/2025-11-14/lowMOIvshighMOI-Validation-Experiment/" + 
        "2025-11-20-dctapseq-validation-wtc11-sig-egpairs-table-s3.tsv", 
        sep="\t",
        dtype={
            "significant_wo_pos_controls_20fdr": str
        }
    )
    return (raw_tbl_s3_wtc11_significant,)


@app.cell
def _(raw_tbl_s3_wtc11_significant):
    pairs_btw_5_and_10 = raw_tbl_s3_wtc11_significant.copy()[
        (raw_tbl_s3_wtc11_significant["pct_change_effect_size"] > -10) &
        (raw_tbl_s3_wtc11_significant["pct_change_effect_size"] < -5)
    ]

    pairs_btw_5_and_10_ids = pairs_btw_5_and_10.loc[
        :, ["element_gene_pair_identifier_hg38", "gene_id", "pct_change_effect_size", "guide_list"]
    ]
    pairs_btw_5_and_10_ids.to_csv("data/2025-12-08/test-pairs-btw-5-and-10.csv", index=False)
    return pairs_btw_5_and_10, pairs_btw_5_and_10_ids


@app.cell
def _(mo, pairs_btw_5_and_10):
    mo.ui.dataframe(pairs_btw_5_and_10)
    return


@app.cell
def _(mo, pairs_btw_5_and_10_ids):
    mo.ui.dataframe(pairs_btw_5_and_10_ids)
    return


if __name__ == "__main__":
    app.run()
