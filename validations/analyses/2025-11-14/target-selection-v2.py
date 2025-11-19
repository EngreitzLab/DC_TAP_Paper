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
    import re
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
def _(pd):
    # Function to add Singleton features

    def merge_singleton(df_s3, df_sing):
        """
        Helper method to merge singleton dataframe to table_s3 dataframe
        """
        # Split comma-separated guide_ids into lists
        df_s3["guide_list"] = df_s3["guide_ids"].str.split(",")

        # Explode to one row per guid
        df_s3_expanded = df_s3.explode("guide_list")
        df_s3_expanded["guide_list"] = df_s3_expanded["guide_list"].str.strip()

        # Merge with df_sing to get effect_size for each guide
        merged = df_s3_expanded.merge(
            df_sing[[
                "grna_id", 
                "response_id", 
                "pct_change_effect_size",
                "n_nonzero_trt",
                "n_nonzero_cntrl",
                "standard_error_pct_change"
            ]],
            left_on=["guide_list", "gene_id"],
            right_on=["grna_id", "response_id"],
            how="left",
        )

        # Re-label column names
        merged = merged.rename(columns={
            "response_id" : "singleton_response_id",
            "pct_change_effect_size_x" : "pct_change_effect_size", # fix original colname
            "pct_change_effect_size_y" : "singleton_pct_change_effect_size",
            "n_nonzero_trt" : "singleton_n_nonzero_trt",
            "n_nonzero_cntrl" : "singleton_n_nonzero_cntrl",
            "standard_error_pct_change_x" : "standard_error_pct_change", # fix original colname
            "standard_error_pct_change_y" : "singleton_standard_error_pct_change"
        })
    
        return merged

    def get_nth_performing_guide(df_merged, n, is_negative_direction):
        """
        Helper method to fetch n-th best performing guide with respect to
        direction of effect.
        """
        df_merged = df_merged.copy()
        is_ascending = True

        if is_negative_direction:
            df_merged = df_merged.loc[df_merged["pct_change_effect_size"] < 0 ]
        else:
            df_merged = df_merged.loc[df_merged["pct_change_effect_size"] > 0 ]
            is_ascending = False
    
        # Pre-process merged df. Sorts guides by singleton % change effect size
        df_nth_guide_features = (
            df_merged
                .sort_values(
                    [
                        "element_gene_pair_identifier_hg38",
                        "singleton_pct_change_effect_size"
                    ],
                    ascending = is_ascending
                ).groupby("element_gene_pair_identifier_hg38")[[
                    "element_gene_pair_identifier_hg38",
                    "grna_id",
                    "singleton_n_nonzero_trt",
                    "singleton_n_nonzero_cntrl",
                    "singleton_pct_change_effect_size",
                    "singleton_standard_error_pct_change"
                ]].head(n)
        )
    
        return df_nth_guide_features

    def format_merged_singleton_guides(df_results, df_ref, n):
        """
        Helper method to format and append statistics and feature columns of the n-th best
        performing guide.
        """
        l = f"singleton_g{n + 1}"               # label
        m = "element_gene_pair_identifier_hg38" # mapper
        df_ref = df_ref.copy().groupby(m).nth(n).set_index(m)

        df_results[f"{l}_guide_id"] = df_results[m].map(df_ref["grna_id"])
        df_results[f"{l}_counts_trt"] = df_results[m].map(df_ref["singleton_n_nonzero_trt"])
        df_results[f"{l}_counts_ctl"] = df_results[m].map(df_ref["singleton_n_nonzero_cntrl"])
        df_results[f"{l}_pctchange_effect_size"] = df_results[m].map(df_ref["singleton_pct_change_effect_size"].round(2))
        df_results[f"{l}_pctchange_stderror"] = df_results[m].map(df_ref["singleton_standard_error_pct_change"].round(2))
        df_results[f"{l}_pctchange_ci_interval"] = df_results[m].map(
            round(df_ref["singleton_standard_error_pct_change"] * 1.96, 2)
        )

    def format_merged_singleton_guides_summary_columns(df_results):
        """
        Helper method to format singleton guide stats in a single summary column.
        """
        sl = f"summary_singleton" # summary label

        # Select singleton columns
        cols = [c for c in df_results.columns if c.startswith("singleton_g")]

        # Group columns
        groups = {}
        for c in cols:
            stats = c.split("_", 2)[-1]
            groups.setdefault(stats, []).append(c)

        # Build summary columns
        for stats, cols in groups.items():
            summary_col = f"{sl}_{stats}"
            df_results[summary_col] = df_results[cols].astype(str).agg(", ".join, axis=1)

    def get_singletonfeatures(df, df_singleton, n = 3):
        """
        Function to get singleton features for each E-G pair.
        Features to get are the top 3 best performing guide 
        (by % change effect size and direction of effect) and 
        their associated statistics.

        Singleton results referenced was analyzed using SCEPTRE tool.
        See DC-TAP-seq Github for details.

        Statistics/Feature Column
        -------------------------
        singleton_*_guide_id : The guide id of the n-th best performing guide.
        singleton_*_guides : The total number of unique guide designed against a tested element. (NOT IMPLEMENTED)
        singleton_*_counts_trt : The total number of nonzero treatment. (N cells where guide was detected)
        singleton_*_counts_ctl : The total number of control treatment. (N cells where guide was not detected)
        singleton_*_pctchange_effect_size : The percentage change effect size of the indicated single guide.
        singleton_*_pctchange_stderror : The percentage change effect size standard error statistics of the 
                                         indicated single guide.
        singleton_*_pctchange_ci_interval : The percentage change effect size 95% confidence interval difference 
                                            of the indicated single guide. (ci_interval = +/- 1.96 * stderr)
                                        
        Note
        ----
        * denotes n-th best performing guide. For example, "g1" is the best performing and 
          "g3" is 3rd best performing guide.
    
        Parameters
        ----------
        df : an existing dataframe to append singleton features onto.
        df_singleton : refers to DC-TAP-seq singleton guide differential expression analysis dataframe via SCEPTRE tool.
        n : the n-th best performing singleton guide feature to fetch.

        Returns
        -------
        Append results to an existing pre-defined dataframe.
        """
        df = df.copy()

        # Merge Singleton annotations to Table_S3
        df_merged = merge_singleton(df, df_singleton)
    
        # Get the n-th best performing guide for when direction of effect is neg/pos
        df_nth_feat_neg = get_nth_performing_guide(df_merged, n, is_negative_direction=True)
        df_nth_feat_pos = get_nth_performing_guide(df_merged, n, is_negative_direction=False)

        # Concat both direction of effect guide features
        df_nth_performing_guides = pd.concat([df_nth_feat_neg, df_nth_feat_pos])

        # For each n-th guide, format and append singleton features to results dataframe
        for i in range(n):
            format_merged_singleton_guides(df, df_nth_performing_guides, i)
        format_merged_singleton_guides_summary_columns(df)
    
        return df
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


if __name__ == "__main__":
    app.run()
