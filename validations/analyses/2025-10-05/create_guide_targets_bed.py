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
    # Create Guide Targets Bed File Notebook
    This notebook returns the results of the initial CRISPRko sgRNA design for DC-TAP-seq paper validation in a friendly format for IGV. A paritial raw metadata file is also created following the definitions set by the IGVF E2G group.

    The goal of creating these files is to quality check and make the final selections of guides to purchase for the validation experiments.
    """
    )
    return


@app.cell
def _():
    import subprocess
    from pathlib import Path

    import pandas as pd
    from pybedtools import BedTool
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    return BedTool, Path, Seq, SeqIO, SeqRecord, pd, subprocess


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Download Data Dependencies""")
    return


@app.cell
def _():
    # Get Data Dependencies
    # subprocess.run(["mkdir", "-p", "data/refs"])
    # subprocess.run([
    #     "curl", 
    #     "-LO",
    #     "--output-dir", "data/refs/", 
    #     "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit"
    # ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Main Notebook""")
    return


@app.cell
def _():
    # Define inputs
    blat_stepSize = 5
    blat_tileSize = 11
    blat_minScore = 21 # Min number of matches. Thus, Num of allowed mismatches = seq_len of spacer_w_pam - minScore. 
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Define Inputs (WTC11)""")
    return


@app.cell
def _(Path, pd):
    wtc11_designs = pd.read_csv("results/2025-10-05/dctapseq-validation-wtc11-sgrna-designs.txt", sep="\t")
    wtc11_itarget_map = pd.read_csv(
        "results/2025-10-05/intermediates/dctapseq-validation-wtc11-sgrna-intermediate-match.tsv", 
        sep="\t"
    )

    out_wtc11_fasta = Path("results/2025-10-05/dctapseq-validation-wtc11-sgrna-designs.fasta")
    out_wtc11_psl = Path("results/2025-10-05/dctapseq-validation-wtc11-sgrna-designs.psl")
    out_wtc11_bed = Path("results/2025-10-05/dctapseq-validation-wtc11-sgrna-designs.bed")
    out_wtc11_cutpos_bed = Path("results/2025-10-05/dctapseq-validation-wtc11-sgrna-cutpositions.bed")
    out_wtc11_filtered_bed = Path("results/2025-10-05/dctapseq-validation-wtc11-sgrna-designs-filtered.bed")
    out_wtc11_metadata = Path("results/2025-10-05/intermediates/guide-metadata-wtc11-intermediate.tsv")
    return (
        out_wtc11_bed,
        out_wtc11_cutpos_bed,
        out_wtc11_fasta,
        out_wtc11_metadata,
        out_wtc11_psl,
        wtc11_designs,
        wtc11_itarget_map,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Define Inputs (K562)""")
    return


@app.cell
def _(Path, pd):
    k562_designs = pd.read_csv("results/2025-10-05/dctapseq-validation-k562-sgrna-designs.txt", sep="\t")
    k562_itarget_map = pd.read_csv(
        "results/2025-10-05/intermediates/dctapseq-validation-k562-sgrna-intermediate-match.tsv", 
        sep="\t"
    )

    out_k562_fasta = Path("results/2025-10-05/dctapseq-validation-k562-sgrna-designs.fasta")
    out_k562_psl = Path("results/2025-10-05/dctapseq-validation-k562-sgrna-designs.psl")
    out_k562_bed = Path("results/2025-10-05/dctapseq-validation-k562-sgrna-designs.bed")
    out_k562_cutpos_bed = Path("results/2025-10-05/dctapseq-validation-k562-sgrna-cutpositions.bed")
    out_k562_filtered_bed = Path("results/2025-10-05/dctapseq-validation-k562-sgrna-designs-filtered.bed")
    out_k562_metadata = Path("results/2025-10-05/intermediates/guide-metadata-k562-intermediate.tsv")
    return (
        k562_designs,
        k562_itarget_map,
        out_k562_bed,
        out_k562_cutpos_bed,
        out_k562_fasta,
        out_k562_metadata,
        out_k562_psl,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Format CRISPick Result Outputs""")
    return


@app.cell
def _(wtc11_designs, wtc11_itarget_map):
    # Get a unique guide_id
    wtc11_designs["guide_id"] = [f"251005_DCTAPseqValidationWTC11_{i}" for i in range(len(wtc11_designs))]

    # Get spacer with pam
    wtc11_designs["spacer_w_pam"] = (
        wtc11_designs.loc[:,"sgRNA Sequence"].astype(str) + 
        wtc11_designs.loc[:,"PAM Sequence"]
    )

    # Get intended target names
    lookup_w = wtc11_itarget_map.set_index("crispick_targets")["intended_target_name"]
    wtc11_designs["intended_target_name"] = wtc11_designs.loc[:,"Input"].map(lookup_w)
    wtc11_designs["ranked_target_name"] = (
        "R" + wtc11_designs.loc[:,"Pick Order"].astype(str) + 
        ":" + wtc11_designs.loc[:,"guide_id"]
    )
    wtc11_designs["chr"] = wtc11_designs.loc[:,"Input"].str.extract(r'^([^:]+)')
    return


@app.cell
def _(k562_designs, k562_itarget_map):
    # Get a unique guide_id
    k562_designs["guide_id"] = [f"251005_DCTAPseqValidationK562_{i}" for i in range(len(k562_designs))]

    # Get spacer with pam
    k562_designs["spacer_w_pam"] = (
        k562_designs.loc[:,"sgRNA Sequence"].astype(str) + 
        k562_designs.loc[:,"PAM Sequence"]
    )

    # Get intended target names
    lookup_k = k562_itarget_map.set_index("crispick_targets")["intended_target_name"]
    k562_designs["intended_target_name"] = k562_designs.loc[:,"Input"].map(lookup_k)
    k562_designs["ranked_target_name"] = (
        "R" + k562_designs.loc[:,"Pick Order"].astype(str) + 
        ":" + k562_designs.loc[:,"guide_id"]
    )
    k562_designs["chr"] = k562_designs.loc[:,"Input"].str.extract(r'^([^:]+)')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Convert to Bed for IGV""")
    return


@app.cell
def _(
    Seq,
    SeqIO,
    SeqRecord,
    k562_designs,
    out_k562_fasta,
    out_wtc11_fasta,
    wtc11_designs,
):
    def create_fasta(df, out):
        records = [
            SeqRecord(
                Seq(seq),
                id=seq_id,
                description=""
            ) for seq, seq_id in zip(df.loc[:,"spacer_w_pam"], df.loc[:,"ranked_target_name"])
        ]
        SeqIO.write(records, out, "fasta")

    create_fasta(wtc11_designs, out_wtc11_fasta)
    create_fasta(k562_designs, out_k562_fasta)
    return


@app.cell
def _(subprocess):
    # BLAT with 2 allowed mismatches (For WTC11, took ~ 2-mins to finish running)
    def blat(fasta, psl, stepSize, tileSize, minScore):
        subprocess.run([
            "blat",
            "data/refs/hg38.2bit",
            fasta,
            psl,
            f"-stepSize={stepSize}",
            f"-tileSize={tileSize}",
            f"-minScore={minScore}"
        ], check=True)

    # blat(out_wtc11_fasta, out_wtc11_psl, blat_stepSize, blat_tileSize, blat_minScore)
    # blat(out_k562_fasta, out_k562_psl, blat_stepSize, blat_tileSize, blat_minScore)
    return


@app.cell
def _(out_k562_bed, out_k562_psl, out_wtc11_bed, out_wtc11_psl, subprocess):
    def convert_psl_to_bed(celltype, psl, out):
        with open(out, "w") as o:
            subprocess.run([
                "echo",
                f'track name="{celltype} Validation Target Guides (all)" '
                f'description="All Ranked Target Guides for DC-TAP-seq knockout validations (hg38) - {celltype}"',
            ], check=True, stdout=o)
            subprocess.run([
                "awk",
                r'NR > 5 {print $14 "\t" $16 "\t" $17 "\t" $10 "\t" $1 "\t" ($9 == "+" ? "+" : "-")}',
                psl
            ], check=True, stdout=o)


    convert_psl_to_bed("WTC11", out_wtc11_psl, out_wtc11_bed)
    convert_psl_to_bed("K562", out_k562_psl, out_k562_bed)
    return


@app.cell
def _(
    BedTool,
    k562_designs,
    out_k562_cutpos_bed,
    out_wtc11_cutpos_bed,
    pd,
    subprocess,
    wtc11_designs,
):
    ## Prediction of where the guide cuts
    def get_cutposition_bed(celltype, df, out):
        bed_df = pd.DataFrame({
            "chr": df.loc[:, "chr"],
            "start": df.loc[:, "sgRNA Cut Position (1-based)"],
            "end": df.loc[:, "sgRNA Cut Position (1-based)"] + 1,
            "name": df.loc[:, "ranked_target_name"],
            "score": ".",
            "strand": df.loc[:,"Strand of sgRNA"]
        })
        BedTool.from_dataframe(bed_df).saveas(out)
        subprocess.run([
            "sed", 
            "-i", 
            "1i"
            f'track name="{celltype} Validation Target Guides Cut Positions (all)" '
            f'description="All Cut Positions of Target Guides for DC-TAP-seq knockout validations (hg38) - {celltype}"',
            out
        ])

    get_cutposition_bed("WTC11", wtc11_designs, out_wtc11_cutpos_bed)
    get_cutposition_bed("K562", k562_designs, out_k562_cutpos_bed)
    return


@app.cell
def _(
    BedTool,
    Path,
    k562_designs,
    out_k562_bed,
    out_k562_metadata,
    out_wtc11_bed,
    out_wtc11_metadata,
    pd,
    wtc11_designs,
):
    def create_intermediate_metadata(celltype, df, bed, out):
        targets = BedTool(
            Path(f"results/2025-10-05/intermediates/dctapseq-validation-{celltype.lower()}-crispick-input-targets.bed")
        )
        guides = BedTool(bed)
        ontargets = guides.intersect(targets, wa=True, f=0.5).to_dataframe()
        offtargets = guides.intersect(targets, v=True).to_dataframe()

        ontargets["guide_id"] = ontargets.loc[:,"name"].str.replace(r'^R\d+:', '', regex=True)
        offtargets["guide_id"] = offtargets.loc[:,"name"].str.replace(r'^R\d+:', '', regex=True)

        # Find number of off-targets
        offtargets_counts = (
            offtargets
                .copy()
                .loc[:,"guide_id"]
                .value_counts(dropna=False)
                .rename_axis("guide_id")
                .reset_index(name="off_targets")
        )

        print(f"{celltype}---------------------------------")
        print("Num on-targets  hits :", len(ontargets))
        print("Num off-targets hits :", len(offtargets))
        print("Total should be", len(guides), ":" ,len(ontargets) + len(offtargets))
        print(f"Total guides that failed to align: {len(df) - len(ontargets)}")
        print("Num of guides that off-targets   :", len(offtargets_counts))

        # Merge off-target metric to on-target bed
        ontargets = (
            ontargets
                .merge(offtargets_counts, on="guide_id", how="left")
                .assign(off_targets=lambda df: df["off_targets"].fillna(0).astype("Int64"))
        )

        # Merge on-target bed to designer dataframe
        df = df.merge(ontargets, on="guide_id", how="left")

        # Create metadata
        metadata_df = pd.DataFrame({
            "guide_id": df.loc[:,"guide_id"],
            "spacer": df.loc[:,"sgRNA Sequence"],
            "pam": df.loc[:,"PAM Sequence"],
            "targeting": df.loc[:,"off_targets"].apply(
                lambda e: True if (pd.notna(e) and str(e).strip()) else False
            ),
            "type": df.loc[:,"off_targets"].apply(
                lambda e: "validation" if (pd.notna(e) and str(e).strip()) else "failed to align"
            ),
            "guide_chr": df.loc[:,"chr"],
            "guide_start": df.loc[:,"start"].astype("Int64"),
            "guide_end": df.loc[:,"end"].astype("Int64"),
            "strand": df.loc[:,"strand"],
            "intended_target_name": df.loc[:,"intended_target_name"],
            "intended_target_chr": df.loc[:,"intended_target_name"].str.replace(r'(:).+', '', regex=True),
            "intended_target_start": df.loc[:,"intended_target_name"].str.replace(r'^[^:]+:(\d+)-\d+$', r'\1', regex=True),
            "intended_target_end": df.loc[:,"intended_target_name"].str.replace(r'^[^:]+:\d+-(\d+)$', r'\1', regex=True),
            "ranked_target_name": df.loc[:,"ranked_target_name"],
            "num_off_targets": df.loc[:,"off_targets"].apply(
                lambda e: e if (pd.notna(e) and str(e).strip()) else pd.NA
            )
        }).to_csv(out, header=True, index=False, sep="\t")

    create_intermediate_metadata("WTC11", wtc11_designs, out_wtc11_bed, out_wtc11_metadata)
    create_intermediate_metadata("K562", k562_designs, out_k562_bed, out_k562_metadata)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Results Summary

    ### IGV

    View the IGV sessions for DC-TAP-seq results using the following below:

    #### Intermediate Sessions

    |Cell Type| IGV Session Link (all tracks)| IGV Session Link (min tracks)| 
    |---------|------------------------------|------------------------------|
    |`WTC11`  |https://tinyurl.com/3hsrhtcc  |https://tinyurl.com/mnh7aj5x  |
    |`K562`   |https://tinyurl.com/49nmfu3p  |https://tinyurl.com/34ymhnvz  |

    #### Final Sessions 

    |Cell Type| IGV Session Link           | 
    |---------|----------------------------|
    |`WTC11`  |TBD                         |
    |`K562`   |TBD                         |

    ### Metadata

    TODO.
    """
    )
    return


if __name__ == "__main__":
    app.run()
