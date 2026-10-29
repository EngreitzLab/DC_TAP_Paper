
# 2026-05-11 Cell Ranger Outputs

## Purpose

This analysis note documents the Cell Ranger count outputs used for downstream validation, the helper script used to download filtered feature-barcode matrices, and the R script used to extract per-sample RNA and CRISPR summary metrics.

## Inputs

- Filtered feature-barcode matrices of the 16 10x channels, from GEO [GSE303901](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303901) (`GSE303901_validation_*_filtered_feature_bc_matrix.h5`). GEO names the two high-MOI samples by their measured MOI (12.6 and 14.4); the code keeps the Cell Ranger sample IDs (`moi11_3`, `moi13`). The mapping is in `metadata/validation_screen_samples.tsv`.
- Per-sample Cell Ranger `metrics_summary.csv` files. These are not public; their values are in `source_data/cellranger_metrics_summary.tsv`.

## Code

- `scripts/bash/download_filtered_h5.sh` downloads the matrices from GEO into `data/filtered_h5/`, named `<sample_id>_filtered_feature_bc_matrix.h5`.
- `analyses/2026-05-11_cellranger/extract_metrics_summary.R` plots the per-sample metrics. By default it reads `source_data/cellranger_metrics_summary.tsv`. If `VALIDATION_CELLRANGER_URI` points to a folder (local or `gs://`) with one Cell Ranger output folder per sample, it re-extracts the metrics from `<sample_id>/outs/metrics_summary.csv`.

Run from `validations/`:

```bash
bash scripts/bash/download_filtered_h5.sh
Rscript analyses/2026-05-11_cellranger/extract_metrics_summary.R
```

## Outputs

- `data/filtered_h5/<sample_id>_filtered_feature_bc_matrix.h5`
- `results/cellranger_june_2026/cellranger_metrics_summary.tsv`
- `results/cellranger_june_2026/pngs/<metric>.png`

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
