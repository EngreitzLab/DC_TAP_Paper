# 2026-05-13 Get Statistics

## Purpose

Summarizes the SCEPTRE results of the validation screen for the 30 selected element–gene pairs. It draws the heatmap of significance in the joint low-MOI analysis, annotated with element category, the screen's effect size and whether the pair was significant in the screen. It also counts significant pairs per analysis and per element category.

## Inputs

- SCEPTRE discovery results and `tables/grna_target_data_frame.tsv` (param `sceptre_dir`, default `analyses/2026-05-11_SCEPTRE/source_data`)
- `metadata/validation_30_pairs.tsv`: the 30 selected pairs and whether each was significant in the screen (the TP53 self-promoter pair, which was outside the screen's FDR set, is marked significant)
- Table S3: `../results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv`
- `data/gencode.v32.primary_assembly.annotation.gtf.gz` (from `scripts/bash/download_references.sh`)

## Code

Run from `validations/`:

```bash
Rscript -e "rmarkdown::render('analyses/2026-05-13_get_statistics/results_statistics.Rmd', output_dir = 'results/2026-05-13_get_statistics')"
```

## Outputs

- `results/2026-05-13_get_statistics/manual_pairs_low_moi_heatmap.pdf`
