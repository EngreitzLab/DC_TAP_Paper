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
    # DC-TAP-seq Validation Target Selection Notebook v2

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

    import matplotlib.pyplot as plt

    from pybedtools import BedTool
    return BedTool, pd, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Main Notebook""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Download data

    Required Files:
    1. DC-TAP-seq Paper Table S3
    2. WTC11 Singleton Analysis Result TSV

    Download Table S3 from [DC-TAP-seq Paper Preprint Supplementary](https://doi.org/10.1101/2025.09.16.676677). 
    Or for convenience, use the `CMD` snippet below at `validation` root directory to download Table S3 from 
    [github](https://github.com/EngreitzLab/DC_TAP_Paper/blob/v1.0.0).

    ```bash
    mkdir -p data/2025-11-14/
    curl -L "https://raw.githubusercontent.com/EngreitzLab/DC_TAP_Paper/refs/tags/v1.0.0/results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv" -o data/2025-11-14/dctapseq-tables3-wtc11.tsv
    ```

    To download the WTC11 sceptre singleton results.

    ```bash
    curl "https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/validations/WTC11/resources/2025-08-10-dctapseq-sceptre-wtc11-singleton.tsv" -o data/2025-11-14/wtc11-singleton.tsv
    ```
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Constants""")
    return


@app.cell
def _():
    CORELIST = [
        "element_gene_pair_identifier_hg38",
        "pct_change_effect_size",
    #     "largest_singleton_pct_change_effect_size",
        "distance_to_gencode_gene_TSS",
        "direct_rate_negative",
        "element_category",
        "ubiq_category",
    #    "gencode_v44_coding_overlap",
        "gencode_protein_coding_gene_body_overlap"
    #    "overlapping_exon",
    #    "rmsk_nearby"
    ]
    return (CORELIST,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Load Data""")
    return


@app.cell
def _(BedTool, pd):
    tbl_s3 = pd.read_csv(
        "data/2025-11-14/dctapseq-tables3-wtc11.tsv", 
        sep="\t",
        dtype={
            "significant_wo_pos_controls_20fdr": str
        }
    )

    wtc11_singleton = pd.read_csv(
        "data/2025-11-14/wtc11-singleton.tsv",
        sep="\t",
        dtype={
            "significant": "boolean"
        }
    )

    """
    For this analysis, repeatmasker bed files were obtained from
    Broad Institute's IGV public database. Find the referenced files at
    https://data.broadinstitute.org/igvdata/annotations/hg38/rmsk/{*.bed.gz}
    """
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
    return rmsk_dict, tbl_s3


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Pre-processing Data""")
    return


@app.cell
def _(tbl_s3):
    def preprocess():
        df = tbl_s3.copy()

        # Filter for significant WTC11 E-G pairs. Include negative and positive % change effect sizes.
        df = df[
            (df["cell_type"].eq("WTC11")) &
            (df["significant_wo_pos_controls_20fdr"].eq("TRUE")) &
            (df["Random_DistalElement_Gene"].eq(True))
        ]

        return df
    return (preprocess,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Functions & Methods""")
    return


@app.cell
def _(BedTool, pd):
    # Function to make a bed file from pre-defined columns
    def mkbed(df):
        """
        Function to make a BedTool Object from from pre-defined columns.
        """
        df = df.copy()
        bed_df = pd.DataFrame({
            "chr": df.loc[:, "targeting_chr_hg38"],
            "start": df.loc[:, "targeting_start_hg38"] - 1, # Convert to 0-based
            "end": df.loc[:, "targeting_end_hg38"],
            "name": df.loc[:, "element_gene_pair_identifier_hg38"],
            "score": ".",
            "strand": "."
        })
        bed = BedTool.from_dataframe(bed_df)
        return bed
    return (mkbed,)


@app.cell
def _(mkbed, rmsk_dict):
    # A set of functions to overlap RepeatMasker features

    def overlap_rmsk(i, q, w = 350):
        """
        Helper method for overlapping RepeatMasker (rmsk) Features.
        Finds nearby rmsk features given an set of genomic coordinate 
        interval inputs and a single rmsk reference (e.g. SINE).

        Parameters
        ----------
        i : pybedtools BedTool Object Input
        q : pybedtools BedTool Object RepeatMasker reference query
        w : defines a nearby window around an interval.
            Default is 350 bp which defines a 1000-bp interval
            for most DC-TAP-seq tested elements. A 1000-bp
            interval for defining nearby rmsk features is an
            arbitrary definition.

        Returns
        -------
        A list of genomic coordiantes names overlapping a nearby rmsk 
        feature within a specified window.
        """
        nearby = i.slop(genome="hg38", l=w, r=w)
        intersect = nearby.intersect(q)
        return [interval.name for interval in intersect]

    def label_overlap(df, label, match_list):
        """
        Helper method for labeling intervals for an existing dataframe
        given a rmsk feature label and a match list.
        """
        label = "rmsk_" + label
        df.loc[df['element_gene_pair_identifier_hg38'].isin(match_list), label] = True
        df.loc[~df['element_gene_pair_identifier_hg38'].isin(match_list), label] = False
        return df

    def find_nearby_rmsk(df, rmsk = rmsk_dict):
        """
        Function to find nearby RepeatMasker (rmsk) Features and append 
        results to an existing pre-defined dataframe.

        Parameters
        ----------
        df : an existing dataframe to append rmsk annotations onto.
        rmsk : refers to repeatmasker definitions and or references.
        """
        df_results = df.copy()
        bed = mkbed(df_results)

        # Find all nearby repetitive genomic regions within a 350-bp window around the element interval 
        # For most cases the total window frame is 1000-bp.
        for k, v in rmsk.items():
            nearby = overlap_rmsk(bed, v)
            df_results = label_overlap(df_results, k, nearby)

        # Label Pair if there exists any nearby repetitive elements
        rmsk_cols = df_results.filter(regex="^rmsk_").columns
        df_results["rmsk_nearby"] = df_results[rmsk_cols].any(axis=1)

        return df_results
    return


@app.cell
def _():
    # Function to add Singleton features
    return


@app.cell
def _():
    # Function to add relevant WTC11 Gene CPM to Pairs with overlapping exons
    # Referenced Dataset 2022-06-25 Elisa' Group TF-Perturb
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Run Scripts & Functions""")
    return


@app.cell
def _(preprocess):
    tbl_s3_preprocessed = preprocess()
    return (tbl_s3_preprocessed,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Exploration""")
    return


@app.cell
def _(mo, tbl_s3_preprocessed):
    # Show all 99 sig pairs.
    mo.ui.dataframe(tbl_s3_preprocessed)
    return


@app.cell
def _(CORELIST, mo, tbl_s3_preprocessed):
    # Filtering for CTCF element category and look at largest magnitude % pct change eff. size
    mo.ui.dataframe(
        tbl_s3_preprocessed[
            (tbl_s3_preprocessed["pct_change_effect_size"] > 10) |
            (tbl_s3_preprocessed["pct_change_effect_size"] < -10)
            ].loc[:, CORELIST]
    )
    return


@app.cell
def _(pd):
    # Check distribution of CPM across all gene from TF-Perturb
    df_cpm = pd.read_csv("data/2025-11-14/tpm_tab_use_for_D0_D2_design_TPMs_from_WTC11_TF_perturb_w_ubiq.csv")
    return (df_cpm,)


@app.cell
def _(df_cpm, mo):
    mo.ui.dataframe(
        df_cpm
    )
    return


@app.cell
def _(df_cpm, plt):
    df_cpm_for_distplot = (
        df_cpm[["gene", "tpm_D0"]]
            .sort_values(by="tpm_D0", ascending=False, na_position="last")[~((df_cpm["gene"].str.contains("MT")))]
            .reset_index(drop=True)
    )
    plt.hist(df_cpm_for_distplot.loc[:,"tpm_D0"], bins=100)
    plt.ylim(0,100)
    plt.show()
    return


@app.cell
def _(mkbed, tbl_s3_preprocessed):
    test_bed = mkbed(tbl_s3_preprocessed)
    test_bed
    return


if __name__ == "__main__":
    app.run()
