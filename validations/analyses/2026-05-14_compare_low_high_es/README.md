# 2026-05-14 Compare Low and High MOI Effect Sizes

## Purpose

Compares effect sizes (% change in expression, with 95% confidence intervals) for the selected pairs that were significant in the screen. Three panels are drawn:

- the screen vs the joint high-MOI analysis;
- the joint high-MOI vs the joint low-MOI analysis;
- the screen vs the joint low-MOI analysis.

Each panel is annotated with the Pearson correlation.

## Inputs

- SCEPTRE discovery results (param `sceptre_dir`, default `analyses/2026-05-11_SCEPTRE/source_data`)
- `metadata/validation_30_pairs.tsv`
- Table S3: `../results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv`
- `data/gencode.v32.primary_assembly.annotation.gtf.gz` (from `scripts/bash/download_references.sh`)

## Code

Run from `validations/`:

```bash
Rscript -e "rmarkdown::render('analyses/2026-05-14_compare_low_high_es/compare_es.Rmd', output_dir = 'results/2026-05-14_compare_low_high_es')"
```

## Outputs

In `results/2026-05-14_compare_low_high_es/`:

- `scatter_effect_size_correlation_high_vs_paper.pdf`
- `scatter_effect_size_correlation_high_vs_low.pdf`
- `scatter_effect_size_correlation_low_vs_paper.pdf`
