# 2026-06-29 CTCF ChIP-seq Validation

## Purpose

Measures CTCF ChIP-seq signal at the FAM111A and CDC37 CTCF elements in WTC11 CRISPRi lines carrying one sgRNA each: two sgRNAs per element, plus two non-targeting control lines. Each line has three ChIP replicates (18 libraries). The signal of a library is the mean CPM over the 301-bp tested element (hg38):

- `FAM111A_peak`: chr11:59158502-59158803
- `CDC37_peak`: chr19:10390714-10391015

## Inputs

- `metadata/ctcf_chip_libraries.tsv`: one row per library (sgRNA, ChIP replicate, sequencing library, GEO sample).
- CTCF ChIP-seq bigWig files from GEO [GSE303901](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303901) (one per GSM, about 830 MB each). They were made from deduplicated alignments with `workflows/bam_to_bigwig.wdl`: genome coverage of read pairs (`bedtools genomecov -pc`) scaled to counts per million mapped reads.
- `CTCF_jaspar_motifs.bed`: CTCF motif hits (JASPAR 2024, MA0139.2) in the two elements, from the UCSC JASPAR 2024 track (hg38).

## Code

Run from `validations/`, either on local copies of the bigWig files:

```bash
bash scripts/bash/download_ctcf_bigwigs.sh
python analyses/2026-06-29_CTCF_ChIP_validation/ctcf_signal_from_bigwigs.py
```

or reading the bigWig files directly from GEO:

```bash
python analyses/2026-06-29_CTCF_ChIP_validation/ctcf_signal_from_bigwigs.py --from-geo
```

Requires Python with `pandas` and `pyBigWig`.

## Outputs

- `results/2026-06-29_CTCF_ChIP_validation/ctcf_signal_per_replicate.tsv`

## Source data

`source_data/ctcf_signal_per_replicate.tsv` is the output of the script. `analyses/2026-06-30_CTCF_Plots/` and `analyses/2026-10-23_ctcf_and_validation_statistics/` read it.
