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
    # Create CRISPick Input Targets Notebook
    This notebook returns a TSV file that can be uploaded to the CRISPick sgRNA Designer web tool given a bed file annotating intended targets.
    """
    )
    return


@app.cell
def _():
    import subprocess
    from pathlib import Path

    import pandas as pd
    from pybedtools import BedTool
    return BedTool, Path, pd


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Review Revelant DC-TAP-seq Paper Results

    View the IGV sessions for DC-TAP-seq results using the following below:

    |Cell Type| IGV Session Link           | 
    |---------|----------------------------|
    |`WTC11`  |https://tinyurl.com/5n8d9nnc|
    |`K562`   |N/A                         |
    """
    )
    return


@app.cell
def _():
    # Get input data
    # subprocess.run(["mkdir", "-p", "data/2025-10-05/"])
    # subprocess.run([
    #     "curl", 
    #     "-LO",
    #     "--output-dir", "data/2025-10-05/", 
    #     "https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/validations/WTC11/2025-10-09-target-selection/dctapseq-validation-wtc11-targets.bed"
    # ])
    # subprocess.run([
    #     "curl", 
    #     "-LO",
    #     "--output-dir", "data/2025-10-05/", 
    #     "https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/validations/K562/2025-10-09-target-selection/dctapseq-validation-k562-targets.bed"
    # ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Main Notebook""")
    return


@app.cell
def _(BedTool, Path):
    # Define inputs
    wtc11_targets = BedTool("data/2025-10-05/dctapseq-validation-wtc11-targets.bed")
    wtc11_out = Path("results/2025-10-05/dctapseq-validation-wtc11-crispick-input-targets.tsv")

    k562_targets = BedTool("data/2025-10-05/dctapseq-validation-k562-targets.bed")
    k562_out = Path("results/2025-10-05/dctapseq-validation-k562-crispick-input-targets.tsv")

    target_expansion_offset = 250 # DC-TAP-seq Targets are 300bp long with offset, final target region spans 800bp
    return (
        k562_out,
        k562_targets,
        target_expansion_offset,
        wtc11_out,
        wtc11_targets,
    )


@app.cell
def _(
    k562_out,
    k562_targets,
    target_expansion_offset,
    wtc11_out,
    wtc11_targets,
):
    # Set offset to expand guide designer target region equally upstream and downstream of target region.
    def create_crispick_inputs(targets, offset, out):
        guide_designer_window = targets.slop(
            genome="hg38", 
            l=offset - 1, # convert bed 0-based to CRISPick format
            r=offset - 1  # convert bed 0-based to CRISPick format
        )

        # Format to satisfy CRISPick input target specification.
        df = guide_designer_window.to_dataframe(header=None, usecols=[0, 1, 2, 5])
        df["targets"] = (df.chrom.astype(str) 
            + ":" + df.strand.astype(str) 
            + ":" + df.start.astype(str) + "-" + df.end.astype(str)
        )

        # Write To TSV. This file is uploaded to CRISPick to design guides for the selected target regions.
        out.parent.mkdir(parents=True, exist_ok=True)
        df.targets.to_csv(out, sep="\t", index=False, header=False)
        return guide_designer_window, df


    window_w, df_w = create_crispick_inputs(wtc11_targets, target_expansion_offset, wtc11_out)
    window_k, df_k = create_crispick_inputs(k562_targets, target_expansion_offset, k562_out)
    return df_k, df_w, window_k, window_w


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Results Summary from CRISPick Web Tool

    Find the web tool at [https://portals.broadinstitute.org/gppx/crispick/public](https://portals.broadinstitute.org/gppx/crispick/public).

    ### CRISPick Input Parameters

    |Parameters      |Selected Argument                        |
    |----------------|-----------------------------------------|
    |Reference Genome|NCBI RefSeq v.GCF_000001405.40-RS_2024_08|
    |Mechanism       |CRISPRko                                 |
    |Enzyme          |SpyoCas9 (NGG)                           |
    |On Target Scorer|RS3seq-Chen2013+RS3target                |
    |CRISPick Quota  |60                                       |

    ### Run Results from CRISPick

    All results are located at `results/2025-10-05/`.

    |Cell Type |Filename                                   |Description                           |
    |----------|-------------------------------------------|--------------------------------------|
    |`WTC11`   |dctapseq-validation-wtc11-sgrna-designs.txt|WTC11 picked sgRNA candidate sequences|
    |`K562`    |dctapseq-validation-k562-sgrna-designs.txt |K562 picked sgRNA candidate sequences | 

    For detailed descriptions of the design file's columns see the official 
    [CRISPick documentation](https://portals.broadinstitute.org/gppx/crispick/public/how-it-works).
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Supplementary Intermediate Files""")
    return


@app.cell
def _(Path, df_k, df_w, k562_targets, pd, window_k, window_w, wtc11_targets):
    # For matching intended target with CRISPick region to design against
    def create_intermediates(i_out, targets, df):
        i_out.parent.mkdir(parents=True, exist_ok=True)

        i_intended_df = targets.to_dataframe(dtype={"start": "Int64", "end": "Int64"}).drop(index=0).reset_index()
        i_match_df = pd.DataFrame({
            "intended_target_name": (
                i_intended_df.chrom.astype(str) 
                + ":" + i_intended_df.start.astype(str) + "-" + i_intended_df.end.astype(str)
            ),
            "crispick_targets": df["targets"]
        })

        i_match_df.to_csv(i_out, sep="\t", index=False)


    w_out = Path("results/2025-10-05/intermediates/dctapseq-validation-wtc11-sgrna-intermediate-match.tsv")
    k_out = Path("results/2025-10-05/intermediates/dctapseq-validation-k562-sgrna-intermediate-match.tsv")

    create_intermediates(w_out, wtc11_targets, df_w)
    create_intermediates(k_out, k562_targets, df_k)

    # For having target CRISPick regions as a bed file
    w_out = Path("results/2025-10-05/intermediates/dctapseq-validation-wtc11-crispick-input-targets.bed")
    window_w.saveas(w_out)
    k_out = Path("results/2025-10-05/intermediates/dctapseq-validation-k562-crispick-input-targets.bed")
    window_k.saveas(k_out)
    return


if __name__ == "__main__":
    app.run()
