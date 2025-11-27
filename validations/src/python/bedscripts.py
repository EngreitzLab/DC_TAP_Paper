import subprocess

import pandas as pd
from pybedtools import BedTool
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import commonscripts as cs


def mkbed(df):
    """
    Function to make a BedTool Object from from pre-defined columns.
    """
    df = df.copy()
    bed_df = pd.DataFrame(
        {
            "chr": df.loc[:, "targeting_chr_hg38"],
            "start": df.loc[:, "targeting_start_hg38"] - 1,  # Convert to 0-based
            "end": df.loc[:, "targeting_end_hg38"],
            "name": df.loc[:, "element_gene_pair_identifier_hg38"],
            "score": ".",
            "strand": ".",
        }
    )
    bed = BedTool.from_dataframe(bed_df)
    return bed


def mkfastq_for_blat(df, seq_colname, seq_id_colname, out):
    """
    Function to make a fastq file to run blat as a downstream process.
    """
    records = [
        SeqRecord(Seq(seq), id=seq_id, description="")
        for seq, seq_id in zip(
            df.loc[:, seq_colname].astype(str), df.loc[:, seq_id_colname]
        )
    ]
    SeqIO.write(records, out, "fasta")


def blat(fasta, out, stepSize=5, tileSize=11, minScore=21):
    """
    Function to make psl file from blat.

    Note
    ----
    blat_minScore of 21 means min number of matches. Thus, Num of allowed
    mismatches = seq_len of spacer_w_pam - minScore.
    """

    ref = "data/refs/hg38.2bit"

    # Check data dependencies
    if not cs.file_exist(ref):
        print("INFO @ Running blat: No hg38.2bit reference file found.")
        print("INFO @ Running blat: Downloading hg38.2bit.")
        subprocess.run(["mkdir", "-p", "data/refs"])
        subprocess.run(
            [
                "curl",
                "-LO",
                "--output-dir",
                "data/refs/",
                "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit",
            ]
        )
        print("INFO @ Running blat: Downloaded.")
        print("INFO @ Running blat: Running blat...")

    subprocess.run(
        [
            "blat",
            "data/refs/hg38.2bit",
            fasta,
            out,
            f"-stepSize={stepSize}",
            f"-tileSize={tileSize}",
            f"-minScore={minScore}",
        ],
        check=True,
    )


def convert_psl_to_bed(psl, out):
    """
    Function to convert psl file to a bed file.
    """
    with open(out, "w") as o:
        subprocess.run(
            [
                "awk",
                r'NR > 5 {print $14 "\t" $16 "\t" $17 "\t" $10 "\t" $1 "\t" ($9 == "+" ? "+" : "-")}',
                psl,
            ],
            check=True,
            stdout=o,
        )
