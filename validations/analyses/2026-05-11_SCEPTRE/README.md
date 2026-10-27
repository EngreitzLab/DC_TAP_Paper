# 2026-05-11 SCEPTRE

## Purpose

Tests the element–gene pairs of the validation screen with SCEPTRE (v0.10.3). Six analyses are run:

| Prefix | Samples | MOI | Covariates (besides per-cell log counts) |
|---|---|---|---|
| `MOI1_1` | replicate 1 | low (1.1) | channel |
| `MOI1_5` | replicate 2 | low (1.5) | channel |
| `MOI11_3` | replicate 1 | high (12.6) | channel |
| `MOI13` | replicate 2 | high (14.4) | channel |
| `low_moi_merged` | both low-MOI replicates | low | channel, sample |
| `high_moi_merged` | both high-MOI replicates | high | channel, sample |

The prefixes are the Cell Ranger sample IDs. The MOIs in parentheses are the measured values used on GEO (see `metadata/validation_screen_samples.tsv`).

`01_create_sceptre_objects.Rmd`:
- reads the filtered Seurat object and keeps the readout genes;
- builds the discovery pairs from Table S3 (all WTC11 element–gene pairs of the 30 targeted elements);
- creates one SCEPTRE object per analysis.

The low-MOI analyses use low-MOI mode; the high-MOI analyses use high-MOI mode.

`02_run_sceptre_pipeline.Rmd` runs gRNA assignment, QC (`n_nonzero_trt_thresh = 15`, `n_nonzero_cntrl_thresh = 15`), the calibration check and the discovery analysis. Tests are two-sided, with gRNAs targeting the same element combined (`union`).

## Inputs

- `results/qc_experiment_June_2026/combined_cdna_seurat_filtered.rds` (from `analyses/2026-05-11_qc_experiment/`)
- Table S3: `../results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv`

## Code

Run from `validations/`; this takes several hours:

```bash
bash scripts/bash/run_validation_pipeline.sh 4
```

The script also builds the Seurat object. To run a single analysis:

```bash
Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/01_create_sceptre_objects.Rmd', params = list(moi_level = 'low', sample_list = c('MOI1_1', 'MOI1_5'), covariate_columns = c('channel', 'sample'), prefix = 'low_moi_merged'), output_dir = 'results/2026-06-10_SCEPTRE')"
Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', params = list(sceptre_object_path = 'results/2026-06-10_SCEPTRE/objects/low_moi_merged_sceptre_object_low_moi_merged.rds', prefix = 'low_moi_merged', integration_strategy = 'union', n_processors = 4), output_dir = 'results/2026-06-10_SCEPTRE')"
```

## Outputs

In `results/2026-06-10_SCEPTRE/`:

- `objects/<prefix>_sceptre_object_<moi>_moi_merged.rds`
- `runs/<prefix>_union/`: SCEPTRE output files, `sceptre_discovery_results.tsv`, `sceptre_calibration_results.tsv`, `sceptre_object_post_pipeline.rds`
- `tables/`: discovery pairs, `grna_target_data_frame.tsv`, per-object and per-run summaries

## Source data

`source_data/` holds the outputs that the downstream analyses use:

- `runs/<prefix>_union/sceptre_discovery_results.tsv`, `sceptre_calibration_results.tsv` and `sceptre_outputs/analysis_summary.txt` for the six analyses
- `tables/grna_target_data_frame.tsv`

`analyses/2026-05-13_get_statistics/`, `analyses/2026-05-14_compare_low_high_es/` and `analyses/2026-10-23_ctcf_and_validation_statistics/` read these tables by default. To use a fresh run instead, set their `sceptre_dir` parameter to `results/2026-06-10_SCEPTRE`.
