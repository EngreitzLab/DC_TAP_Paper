# 2026-05-11 QC Experiment

## Purpose

This analysis note documents the workflow used to build a merged Seurat object from filtered Cell Ranger H5 files, attach sample and feature metadata, and run QC comparisons including targeted-gene and scRNA-reference checks.

## Inputs

Core inputs:

- `data/filtered_h5/`
- `data/terra_table_fixed_pairing.tsv`
- `data/wtc11_encode_feature_reference_file_dups_removed_FIXED.csv`
- `data/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv`
- `data/Table_S4_primer_gene_list.tsv`
- `data/gencode.v32.primary_assembly.annotation.gtf.gz`
- `data/D0_10x5prime_pilot_counts_summary_ribo15_mt15_umi2000.tsv`
- `results/guide_pool_qc/guide_counts.tsv`

## Code

Main notebooks:

- `analyses/2026-05-11_qc_experiment/combine_h5_counts.Rmd`
- `analyses/2026-05-11_qc_experiment/qc_experiment.Rmd`

Shared helper code:

- `scripts/R/utils.R`

Render the combined object notebook from the repository root with:

```bash
RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_qc_experiment/combine_h5_counts.Rmd', output_format = 'html_document', output_dir = 'results/qc_experiment', output_file = 'combine_h5_counts.html')"
```

Render the QC notebook from the repository root with:

```bash
RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_qc_experiment/qc_experiment.Rmd', output_format = 'html_document', output_dir = 'results/qc_experiment', output_file = 'qc_experiment.html')"
```

## Outputs

Merged object and annotation cache:

- `results/qc_experiment/combined_cdna_seurat.rds`
- `results/qc_experiment/gencode_v32_gene_types.tsv`

Notebook outputs:

- `results/qc_experiment/combine_h5_counts.html`
- `results/qc_experiment/qc_experiment.html`

QC comparison outputs:

- `results/qc_experiment/compare_to_reference_scRNA.png`
- `results/qc_experiment/compare_to_reference_scRNA.pdf`
- `results/qc_experiment/gene_umi_comparison_quadrants.tsv`

## Tool versions and settings

- Seurat / SeuratObject workflow with RNA + gRNA assays
- GTF import via `rtracklayer`
- Data processing with `data.table` and sparse matrices from `Matrix`
- Rendering via `rmarkdown` with `RSTUDIO_PANDOC`
