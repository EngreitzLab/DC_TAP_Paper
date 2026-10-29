# 2026-05-11 QC Experiment

## Purpose

Builds one merged Seurat object from the 16 filtered Cell Ranger matrices, attaches sample and feature metadata, keeps the guides that passed the plasmid-pool QC, and computes cell QC metrics. It also compares mean UMIs per targeted gene with a reference scRNA-seq dataset.

## Inputs

- `data/filtered_h5/`: from `scripts/bash/download_filtered_h5.sh` (GEO GSE303901)
- `data/GSE303901_wtc11_feature_reference.csv` and `data/gencode.v32.primary_assembly.annotation.gtf.gz`: from `scripts/bash/download_references.sh`
- `metadata/validation_screen_samples.tsv`, `metadata/Table_S4_primer_gene_list.tsv`
- Table S3: `../results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv`
- `analyses/2026-05-11_guide_pool_qc/source_data/guide_counts.tsv`
- Reference scRNA-seq: `data/D0_10x5prime_pilot_counts_summary_ribo15_mt15_umi2000.tsv`, an internal 10x 5' dataset that is not public. Without it, the comparison is plotted from `source_data/gene_umi_comparison_quadrants.tsv`.

## Code

- `combine_h5_counts.Rmd`: merged Seurat object with sample metadata, gene IDs/symbols and guide metadata.
- `qc_experiment.Rmd`: guide filtering, QC metrics and the reference comparison.
- Shared functions: `scripts/R/utils.R`.

Run from `validations/`:

```bash
Rscript -e "rmarkdown::render('analyses/2026-05-11_qc_experiment/combine_h5_counts.Rmd', output_dir = 'results/qc_experiment_June_2026')"
Rscript -e "rmarkdown::render('analyses/2026-05-11_qc_experiment/qc_experiment.Rmd', output_dir = 'results/qc_experiment_June_2026')"
```

## Outputs

In `results/qc_experiment_June_2026/`:

- `combined_cdna_seurat.rds`: merged object
- `combined_cdna_seurat_filtered.rds`: object with QC metrics and pool-QC guides only (input to SCEPTRE)
- `gencode_v32_gene_types.tsv`: gene-type cache
- `compare_to_reference_scRNA.{png,pdf}`, `gene_umi_comparison_quadrants.tsv`
