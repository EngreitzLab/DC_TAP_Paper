"""Mean CTCF ChIP-seq signal over the FAM111A and CDC37 CTCF elements, per library.

Reads the 18 CTCF ChIP-seq bigWig files (CPM-normalized coverage, see
workflows/bam_to_bigwig.wdl) listed in metadata/ctcf_chip_libraries.tsv and writes
the mean signal over each 301-bp element.

Usage (from validations/), either with local copies of the bigWig files:
    bash scripts/bash/download_ctcf_bigwigs.sh
    python analyses/2026-06-29_CTCF_ChIP_validation/ctcf_signal_from_bigwigs.py
or reading the bigWig files directly from GEO:
    python analyses/2026-06-29_CTCF_ChIP_validation/ctcf_signal_from_bigwigs.py --from-geo
"""

import argparse
from pathlib import Path

import pandas as pd
import pyBigWig

VALIDATIONS_DIR = Path(__file__).resolve().parents[2]
GEO_SAMPLES_URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples"

ELEMENTS = {
    "FAM111A_peak": ("chr11", 59158502, 59158803),
    "CDC37_peak": ("chr19", 10390714, 10391015),
}

GRNA_IDS = {
    "G2": "WTC11_Random_Screen_Crop_3864",
    "G3": "WTC11_Random_Screen_Crop_3866",
    "G4": "WTC11_Random_Screen_Crop_7494",
    "G6": "WTC11_Random_Screen_Crop_7498",
}

INTENDED_TARGETS = {
    "FAM111A_CTCF_G2": "FAM111A",
    "FAM111A_CTCF_G3": "FAM111A",
    "CDC37_CTCF_G4": "CDC37",
    "CDC37_CTCF_G6": "CDC37",
    "NC2": "control",
    "WTC11_NT": "control",
}


def mean_signal(bigwig_path, chrom, start, end):
    bw = pyBigWig.open(bigwig_path)
    value = bw.stats(chrom, start, end, type="mean")[0] or 0.0
    bw.close()
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--libraries",
        default=VALIDATIONS_DIR / "metadata" / "ctcf_chip_libraries.tsv",
        type=Path,
    )
    parser.add_argument(
        "--bigwig-dir",
        default=VALIDATIONS_DIR / "data" / "ctcf_chip_bigwigs",
        type=Path,
    )
    parser.add_argument(
        "--from-geo",
        action="store_true",
        help="read the bigWig files from GEO instead of --bigwig-dir",
    )
    parser.add_argument(
        "--output",
        default=VALIDATIONS_DIR / "results" / "2026-06-29_CTCF_ChIP_validation" / "ctcf_signal_per_replicate.tsv",
        type=Path,
    )
    args = parser.parse_args()

    libraries = pd.read_csv(args.libraries, sep="\t")

    rows = []
    for library in libraries.itertuples():
        if args.from_geo:
            gsm = library.geo_gsm
            bigwig_path = f"{GEO_SAMPLES_URL}/{gsm[:7]}nnn/{gsm}/suppl/{library.geo_bigwig}"
        else:
            bigwig_path = str(args.bigwig_dir / library.geo_bigwig)
        for locus_label, (chrom, start, end) in ELEMENTS.items():
            rows.append(
                {
                    "sample_id": library.sample_id,
                    "group": library.group,
                    "locus": f"{chrom}:{start}-{end}",
                    "locus_label": locus_label,
                    "mean_signal": mean_signal(bigwig_path, chrom, start, end),
                }
            )

    signal = pd.DataFrame(rows)
    signal["replicate"] = signal["sample_id"].str.extract(r"_(rep\d+)$", expand=False)
    signal["guide"] = signal["group"].str.extract(r"_(G\d+)$", expand=False)
    signal["guide"] = signal["guide"].fillna(signal["group"])
    signal["grna_id"] = signal["guide"].map(GRNA_IDS)
    signal["intended_target"] = signal["group"].map(INTENDED_TARGETS)

    signal = (
        signal[["sample_id", "group", "guide", "grna_id", "replicate", "intended_target", "locus", "locus_label", "mean_signal"]]
        .sort_values(["locus", "group", "sample_id"])
        .reset_index(drop=True)
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    signal.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(signal)} rows to {args.output}")


if __name__ == "__main__":
    main()
