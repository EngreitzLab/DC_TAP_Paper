import marimo

__generated_with = "0.16.5"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Select Targets Notebook

    This notebook explores the metadata from DC-TAP-seq manuscript Table S3 to appropriately select targets.
    """
    )
    return


@app.cell
def _():
    import subprocess
    from pathlib import Path

    import math
    import numpy as np
    import pandas as pd

    import itables
    import altair as alt
    import matplotlib.pyplot as plt

    from pybedtools import BedTool
    return BedTool, itables, pd


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Main Notebook""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Notebook Setup""")
    return


@app.cell
def _(itables):
    # Setup
    itables.init_notebook_mode(all_interactive=True)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Load Data""")
    return


@app.cell
def _(pd):
    # Get DC-TAP-seq Manuscript Table S3
    wtc_poshits = pd.read_csv("data/2025-11-10/dctapseq-tables3-wtc11-poshits.tsv", sep="\t")
    return (wtc_poshits,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Exploration""")
    return


@app.cell
def _(wtc_poshits):
    wtc_poshits
    return


@app.cell
def _(BedTool):
    wtc11_bed = BedTool("data/2025-11-10/dctapseq-tables3-wtc11-poshits.bed")
    rmsk_dict = {
        "DNA" : BedTool("data/2025-11-10/rmsk/DNA.bed.gz"),
        "LINE" : BedTool("data/2025-11-10/rmsk/LINE.bed.gz"),
        "LTR" : BedTool("data/2025-11-10/rmsk/LTR.bed.gz"),
        "RC" : BedTool("data/2025-11-10/rmsk/RC.bed.gz"),
        "Retroposon" : BedTool("data/2025-11-10/rmsk/Retroposon.bed.gz"),
        "RNA" : BedTool("data/2025-11-10/rmsk/RNA.bed.gz"),
        "rRNA" : BedTool("data/2025-11-10/rmsk/rRNA.bed.gz"),
        "Satellite" : BedTool("data/2025-11-10/rmsk/Satellite.bed.gz"),
        "scRNA" : BedTool("data/2025-11-10/rmsk/scRNA.bed.gz"),
        "Simplerepeats" : BedTool("data/2025-11-10/rmsk/Simple_repeat.bed.gz"),
        "SINE" : BedTool("data/2025-11-10/rmsk/SINE.bed.gz"),
        "snRNA" : BedTool("data/2025-11-10/rmsk/snRNA.bed.gz"),
        "srpRNA" : BedTool("data/2025-11-10/rmsk/srpRNA.bed.gz"),
        "tRNA" : BedTool("data/2025-11-10/rmsk/tRNA.bed.gz"),
        "UNKNOWN" : BedTool("data/2025-11-10/rmsk/UNKNOWN.bed.gz")
    }
    return rmsk_dict, wtc11_bed


@app.cell
def _():
    def overlap_rmsk(i, q):
        # 350 left and right because most elements 300bp in length
        nearby = i.slop(genome="hg38", l=350, r=350)
        intersect = nearby.intersect(q)
        return [interval.name for interval in intersect]

    def label_overlap(df, label, match_list):
        label = "rmsk_" + label
        df.loc[df['element_gene_pair_identifier_hg38'].isin(match_list), label] = True
        df.loc[~df['element_gene_pair_identifier_hg38'].isin(match_list), label] = False
        return df

    def find_nearby_rmsk(df, bed, rmsk):
        df_results = df.copy()

        # Find all nearby repetitive genomic regions within 1000bp of element center 
        for k, v in rmsk.items():
            nearby = overlap_rmsk(bed, v)
            df_results = label_overlap(df_results, k, nearby)

        # Label Pair if there exists any nearby repetitive elements
        rmsk_cols = df_results.filter(regex="^rmsk_").columns
        df_results["rmsk_nearby"] = df_results[rmsk_cols].any(axis=1)

        return df_results
    return (find_nearby_rmsk,)


@app.cell
def _(find_nearby_rmsk, rmsk_dict, wtc11_bed, wtc_poshits):
    df1 = find_nearby_rmsk(wtc_poshits, wtc11_bed, rmsk_dict)
    return (df1,)


@app.cell
def _(df1):
    df1
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""#### Adding Singleton Data""")
    return


@app.cell
def _(pd, wtc_poshits):
    wtc_sing = pd.read_csv("data/2025-11-10/wtc11_singleton.tsv", sep="\t", dtype={11: "boolean"}).loc[:, ["grna_id", "response_id", "pct_change_effect_size"]]
    wtc_poshits_sing = wtc_poshits.copy().loc[:,["element_gene_pair_identifier_hg38", "gene_id", "guide_ids", "pct_change_effect_size"]]
    return wtc_poshits_sing, wtc_sing


@app.cell
def _(wtc_sing):
    wtc_sing
    return


@app.cell
def _(wtc_poshits_sing):
    wtc_poshits_sing
    return


@app.function
def find_largest_singleton_effectsize(df1, df2):
    df1 = df1.copy()

    # Split comma-separated guide_ids into lists
    df1["guide_list"] = df1["guide_ids"].str.split(",")

    # Explode to one row per guid
    df1_expanded = df1.explode("guide_list")
    df1_expanded["guide_list"] = df1_expanded["guide_list"].str.strip()

    # Merge with df2 to get effect_size for each guide
    merged = df1_expanded.merge(
        df2[["grna_id", "response_id", "pct_change_effect_size"]],
        left_on=["guide_list", "gene_id"],
        right_on=["grna_id", "response_id"],
        how="left",
    )
    # For each original df1 row, take the large effect_size in negative direction
    largest_effect_per_pair = (
        merged
        .groupby("element_gene_pair_identifier_hg38")["pct_change_effect_size_y"]
        .min()
    )
    print(largest_effect_per_pair)

    # Attach back to df1
    df1["largest_effect_size"] = df1["element_gene_pair_identifier_hg38"].map(largest_effect_per_pair)
    return df1


@app.cell
def _(wtc_poshits_sing, wtc_sing):
    df_test = find_largest_singleton_effectsize(wtc_poshits_sing, wtc_sing)
    return (df_test,)


@app.cell
def _(df_test):
    df_test
    return


@app.cell
def _(df1, df_test):
    df1["largest_singleton_pct_change_effect_size"] = df_test.loc[:,"largest_effect_size"]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Exploration Part 2""")
    return


@app.cell
def _():
    df1_corelist = [
        "element_gene_pair_identifier_hg38",
        "pct_change_effect_size",
        "largest_singleton_pct_change_effect_size",
        "distance_to_gencode_gene_TSS",
        "direct_rate_negative",
        "element_category",
        "ubiq_category",
        "gencode_v44_coding_overlap",
        "overlapping_exon",
        "rmsk_nearby"
    ]

    df1_corelist_w_rmsk = [
        "element_gene_pair_identifier_hg38",
        "pct_change_effect_size",
        "largest_singleton_pct_change_effect_size",
        "distance_to_gencode_gene_TSS",
        "direct_rate_negative",
        "element_category",
        "ubiq_category",
        "gencode_v44_coding_overlap",
        "overlapping_exon",
        "rmsk_nearby",
        "rmsk_DNA",
        "rmsk_LINE",
        "rmsk_LTR",
        "rmsk_RC",
        "rmsk_Retroposon",
        "rmsk_RNA",
        "rmsk_rRNA",
        "rmsk_Satellite",
        "rmsk_scRNA",
        "rmsk_Simplerepeats",
        "rmsk_SINE",
        "rmsk_snRNA",
        "rmsk_srpRNA",
        "rmsk_tRNA",
        "rmsk_UNKNOWN"
    ]
    return (df1_corelist,)


@app.cell
def _(df1, df1_corelist):
    df1.loc[:, df1_corelist]
    return


@app.cell
def _(df1, df1_corelist):
    # From 66 Pairs
    # → 38 Pairs ⇒ less than 100kb from TSS. 
    # Ranges (1-59 kb). 
    # Pair that is 59 kb away from TSS corresponds to Direct Hit Rate of 2.38%

    # → 28 Pairs ⇒ more than 100kb from TSS. 
    # Ranges (150-197 kb).
    # Pair that is 150 kb away from TSS corresponds to Direct Hit Rate of 0.691%

    # Pairs where distance to TSS is less than 100kb away
    df1.loc[df1["distance_to_gencode_gene_TSS"] < 100000,:].loc[:, df1_corelist]
    return


@app.cell
def _(df1, df1_corelist):
    # Pairs where % Effect Size is less than -10% & avoiding overlaps 
    df1.loc[
        (df1["pct_change_effect_size"] < -10) & 
        (df1["distance_to_gencode_gene_TSS"] < 100000) &
        (df1["overlapping_exon"] == False)
    ].loc[:, df1_corelist]
    return


@app.cell
def _(df1, df1_corelist):
    # Pairs where % Effect Size is less than -10% (singleton) & avoiding overlaps 
    df1.loc[
        (df1["largest_singleton_pct_change_effect_size"] < -10) &
        (df1["distance_to_gencode_gene_TSS"] < 100000) &
        (df1["overlapping_exon"] == False)
    ].loc[:, df1_corelist]
    return


if __name__ == "__main__":
    app.run()
