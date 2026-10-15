# 2026-05-11 SCEPTRE

## Purpose

This analysis note documents creation and execution of per-sample SCEPTRE objects for the DC-TAP validation dataset.

The workflow:

- reads the filtered Seurat object
- keeps only readout genes in the RNA assay
- uses existing gRNA-to-element metadata from the Seurat object
- builds discovery pairs from Table S3 for WTC11 using all element-to-gene pairs (including non-significant pairs)
- creates four independent SCEPTRE objects, one per sample
- runs the SCEPTRE pipeline separately for each sample

MOI grouping used in this analysis:

- high MOI: `MOI11_3`, `MOI13`
- low MOI: `MOI1_1`, `MOI1_5`

## Inputs

- `results/qc_experiment/combined_cdna_seurat.rds`
- `data/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv`

Required R package:

- `sceptre`

## Code

Object creation notebook (parameterized):

- `analyses/2026-05-11_SCEPTRE/01_create_sceptre_objects.Rmd`

Pipeline execution notebooks:

- `analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd` - Per-sample discovery
- `analyses/2026-05-11_SCEPTRE/03_run_sceptre_pipeline_singletons.Rmd` - Per-sample with singleton gRNA integration
- `analyses/2026-05-11_SCEPTRE/04_run_sceptre_pipeline_merge_replicates.Rmd` - Merged replicate analysis

### SCEPTRE Object Creation

All objects are created via `01_create_sceptre_objects.Rmd` with the following parameters:
- `moi_level`: "low" or "high" (determines gRNA assignment method: maximum vs. mixture)
- `sample_list`: c(...) vector of sample names to include
- `covariate_columns`: c(...) vector of covariates to include in the formula (always includes log-transformed expression metrics)
- `prefix`: output prefix for the .rds file

**Per-sample objects** (single sample, channel covariate only):

```bash
Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/01_create_sceptre_objects.Rmd', \
  params=list(moi_level='low', sample_list=c('MOI1_1'), covariate_columns=c('channel'), prefix='MOI1_1'), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='01_per_sample_moi1_1.html')"

Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/01_create_sceptre_objects.Rmd', \
  params=list(moi_level='low', sample_list=c('MOI1_5'), covariate_columns=c('channel'), prefix='MOI1_5'), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='01_per_sample_moi1_5.html')"

Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/01_create_sceptre_objects.Rmd', \
  params=list(moi_level='high', sample_list=c('MOI11_3'), covariate_columns=c('channel'), prefix='MOI11_3'), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='01_per_sample_moi11_3.html')"

Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/01_create_sceptre_objects.Rmd', \
  params=list(moi_level='high', sample_list=c('MOI13'), covariate_columns=c('channel'), prefix='MOI13'), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='01_per_sample_moi13.html')"
```

**Merged low-MOI object** (MOI1_1 + MOI1_5, with sample batch covariate):

```bash
Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/01_create_sceptre_objects.Rmd', \
  params=list(moi_level='low', sample_list=c('MOI1_1', 'MOI1_5'), covariate_columns=c('channel', 'sample'), prefix='low_moi_merged'), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='01_merged_low_moi.html')"
```

**Merged high-MOI object** (MOI11_3 + MOI13, with sample batch covariate):

```bash
Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/01_create_sceptre_objects.Rmd', \
  params=list(moi_level='high', sample_list=c('MOI11_3', 'MOI13'), covariate_columns=c('channel', 'sample'), prefix='high_moi_merged'), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='01_merged_high_moi.html')"
```

### Per-Sample Discovery Pipeline (Notebook 02)

Run the SCEPTRE discovery pipeline for each sample in **union** mode:

```bash
RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI1_1_sceptre_object_low_moi_merged.rds', \
              prefix='MOI1_1', integration_strategy='union', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI1_1_union.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI1_5_sceptre_object_low_moi_merged.rds', \
              prefix='MOI1_5', integration_strategy='union', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI1_5_union.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI11_3_sceptre_object_high_moi_merged.rds', \
              prefix='MOI11_3', integration_strategy='union', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI11_3_union.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI13_sceptre_object_high_moi_merged.rds', \
              prefix='MOI13', integration_strategy='union', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI13_union.html')" 2>&1 | tail -30
```

Run the SCEPTRE discovery pipeline for each sample in **singleton** mode:

```bash
RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI1_1_sceptre_object_low_moi_merged.rds', \
              prefix='MOI1_1', integration_strategy='singleton', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI1_1_singleton.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI1_5_sceptre_object_low_moi_merged.rds', \
              prefix='MOI1_5', integration_strategy='singleton', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI1_5_singleton.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI11_3_sceptre_object_high_moi_merged.rds', \
              prefix='MOI11_3', integration_strategy='singleton', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI11_3_singleton.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI13_sceptre_object_high_moi_merged.rds', \
              prefix='MOI13', integration_strategy='singleton', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI13_singleton.html')" 2>&1 | tail -30
```
Run the SCEPTRE discovery pipeline for **merged low-MOI** (MOI1_1 + MOI1_5) in both modes:

```bash
Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/low_moi_merged_sceptre_object_low_moi_merged.rds', \
              prefix='low_moi_merged', integration_strategy='union', n_processors=2), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_low_moi_merged_union.html')" 2>&1 | tail -30

Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/low_moi_merged_sceptre_object_low_moi_merged.rds', \
              prefix='low_moi_merged', integration_strategy='singleton', n_processors=2), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_low_moi_merged_singleton.html')" 2>&1 | tail -30
```

Run the SCEPTRE discovery pipeline for **merged high-MOI** (MOI11_3 + MOI13) in both modes:

```bash
Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/high_moi_merged_sceptre_object_high_moi_merged.rds', \
              prefix='high_moi_merged', integration_strategy='union', n_processors=2), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_high_moi_merged_union.html')" 2>&1 | tail -30

Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/high_moi_merged_sceptre_object_high_moi_merged.rds', \
              prefix='high_moi_merged', integration_strategy='singleton', n_processors=2), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_high_moi_merged_singleton.html')" 2>&1 | tail -30
```

## Outputs

Root output directory:

- `results/2026-05-11_SCEPTRE/`

Notebook outputs:

- `results/2026-05-11_SCEPTRE/01_per_sample_moi*.html` - Per-sample object creation reports
- `results/2026-05-11_SCEPTRE/01_merged_*.html` - Merged object creation reports
- `results/2026-05-11_SCEPTRE/02_*.html` - Per-sample discovery pipeline outputs
- `results/2026-05-11_SCEPTRE/03_*.html` - Per-sample singleton pipeline outputs
- `results/2026-05-11_SCEPTRE/04_*.html` - Merged replicate pipeline outputs

SCEPTRE objects:

- `results/2026-05-11_SCEPTRE/objects/sceptre_object_MOI1_1.rds`
- `results/2026-05-11_SCEPTRE/objects/sceptre_object_MOI1_5.rds`
- `results/2026-05-11_SCEPTRE/objects/sceptre_object_MOI11_3.rds`
- `results/2026-05-11_SCEPTRE/objects/sceptre_object_MOI13.rds`
- `results/2026-05-11_SCEPTRE/objects/sceptre_object_low_moi_merged.rds`
- `results/2026-05-11_SCEPTRE/objects/sceptre_object_high_moi_merged.rds`

Shared tables:

- `results/2026-05-11_SCEPTRE/tables/discovery_pairs_wtc11_all.tsv`
- `results/2026-05-11_SCEPTRE/tables/grna_target_data_frame.tsv`
- `results/2026-05-11_SCEPTRE/tables/sceptre_objects_summary.tsv`
- `results/2026-05-11_SCEPTRE/tables/sceptre_pipeline_*.tsv`

Per-sample run outputs:

- `results/2026-05-11_SCEPTRE/runs/MOI1_1/`
- `results/2026-05-11_SCEPTRE/runs/MOI1_5/`
- `results/2026-05-11_SCEPTRE/runs/MOI11_3/`
- `results/2026-05-11_SCEPTRE/runs/MOI13/`
- `results/2026-05-11_SCEPTRE/runs/MOI*_singleton/`

Merged run outputs:

- `results/2026-05-11_SCEPTRE/runs/merged_lowmoi_union/`
- `results/2026-05-11_SCEPTRE/runs/merged_lowmoi_singleton/`
- `results/2026-05-11_SCEPTRE/runs/merged_highmoi_union/`
- `results/2026-05-11_SCEPTRE/runs/merged_highmoi_singleton/`

Each run folder contains SCEPTRE output files plus:

- `sceptre_object_post_pipeline.rds`
- `sceptre_discovery_results.tsv`
- `sceptre_calibration_results.tsv`

## Notes

If script execution fails with `there is no package called 'sceptre'`, install it in the active R environment before running:

```r
renv::install('katsevich-lab/sceptre')
```
