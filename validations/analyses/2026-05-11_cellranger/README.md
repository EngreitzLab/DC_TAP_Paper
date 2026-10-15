
# 2026-05-11 Cell Ranger Outputs

## Purpose

This analysis note documents the Cell Ranger count outputs used for downstream validation, the helper script used to download filtered feature-barcode matrices, and the R script used to extract per-sample RNA and CRISPR summary metrics.

## Inputs

Base GCS location:

- `gs://<validation-cellranger-bucket>/cellranger_outputs`

Example run directory:

- `gs://<validation-cellranger-bucket>/cellranger_outputs/validation_moi11_3_rep01_a7/`

Example files within a run:

- `gs://<validation-cellranger-bucket>/cellranger_outputs/validation_moi11_3_rep01_a7/outs/filtered_feature_bc_matrix.h5`
- `gs://<validation-cellranger-bucket>/cellranger_outputs/validation_moi11_3_rep01_a7/outs/metrics_summary.csv`
- `gs://<validation-cellranger-bucket>/cellranger_outputs/validation_moi11_3_rep01_a7/outs/possorted_genome_bam.bam`

The download script currently targets all listed validation runs under the bucket and copies the filtered H5 files locally.

## Code

Download script:

- `scripts/bash/download_filtered_h5.sh`

Metrics extraction script:

- `analyses/2026-05-11_cellranger/extract_metrics_summary.R`

Run from the repository root:

```bash
./scripts/bash/download_filtered_h5.sh
```

Generate the per-sample metrics table from the repository root with:

```bash
PATH='/usr/bin:/bin:/usr/sbin:/sbin:<path-to-google-cloud-sdk>/bin' \
/usr/local/bin/Rscript analyses/2026-05-11_cellranger/extract_metrics_summary.R
```

## Outputs

Local downloaded matrices:

- `data/filtered_h5/`

Extracted metrics table:

- `results/cellranger/cellranger_metrics_summary.tsv`

Each downloaded file is renamed to include the sample ID, for example:

- `data/filtered_h5/validation_moi11_3_rep01_a7_filtered_feature_bc_matrix.h5`

The extracted table standardizes the following fields from `metrics_summary.csv`:

- `cells_recovered`
- `mRNA_reads`
- `mean_mRNA_reads_per_cell`
- `Reads_mapped_to_targeted_transcriptome`
- `Median_targeted_mRNA_UMI_per_cell`
- `mRNA_seq_saturation`
- `guide_reads`
- `percent_guide_reads_with_protospacer`
- `median_guide_UMIs_per_cell`
- `guide_sequencing_saturation`

For this `cellranger count` workflow, `Reads_mapped_to_targeted_transcriptome` is populated from the standard Cell Ranger metric `Reads Mapped Confidently to Transcriptome`, and `Median_targeted_mRNA_UMI_per_cell` is populated from `Median UMI Counts per Cell`.

## Tool versions and settings

- Cell Ranger: `v9.0.1`
- Chemistry: `SC3Pv4`
- Transcriptome reference: `refdata-gex-GRCh38-2020-A`
