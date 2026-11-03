# 2026-10-23 CTCF and Validation Statistics

## Purpose

- `01_ctcf_chip_statistics.Rmd`: tests of the CTCF ChIP-seq signal. Each targeting sgRNA (3 libraries) and both sgRNAs of a locus (6 libraries) are compared with the pooled non-targeting controls (6 libraries) by exact two-sided Wilcoxon rank-sum tests, with Holm and BH adjustment. Sensitivity analyses: each control line separately, a Welch t-test on log2 CPM, and the lines targeting the other locus. Also reports the per-guide CRISPRi effects behind the per-guide plots and the overlap of the sgRNAs with CTCF motifs.
- `02_validation_30_pairs_table.Rmd`: original screen and validation estimates of the 30 selected element–gene pairs. It covers the pooled low- and high-MOI runs and the two low-MOI replicates analysed separately, with correlations with the screen and replicate concordance.

## Inputs

- `analyses/2026-06-29_CTCF_ChIP_validation/source_data/ctcf_signal_per_replicate.tsv`
- `metadata/ctcf_chip_libraries.tsv`
- `analyses/2026-06-30_CTCF_Plots/chip_sgrnas_hg38.bed`, `analyses/2026-06-30_CTCF_Plots/ctcf_element_motifs_hg38.bed`, `analyses/2026-06-29_CTCF_ChIP_validation/CTCF_jaspar_motifs.bed`
- Singleton SCEPTRE results and guide positions of the WTC11 screen, and `data/hg19ToHg38.over.chain` (as in `analyses/2026-06-30_CTCF_Plots/`)
- SCEPTRE results of the validation screen (param `sceptre_dir`, default `analyses/2026-05-11_SCEPTRE/source_data`)
- `metadata/validation_30_pairs.tsv`
- Table S3: `../results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv`

## Code

Run from `validations/`:

```bash
bash scripts/bash/download_references.sh
Rscript -e "rmarkdown::render('analyses/2026-10-23_ctcf_and_validation_statistics/01_ctcf_chip_statistics.Rmd', output_dir = 'results/2026-10-23_ctcf_and_validation_statistics')"
Rscript -e "rmarkdown::render('analyses/2026-10-23_ctcf_and_validation_statistics/02_validation_30_pairs_table.Rmd', output_dir = 'results/2026-10-23_ctcf_and_validation_statistics')"
```

## Outputs

Both notebooks write their tables to `source_data/`:

- `ctcf_chip_signal_per_library.tsv`: CTCF signal per library, with sgRNA, replicate, sequencing library and fragment count
- `ctcf_chip_wilcoxon_tests.tsv`: guide-level and locus-level tests
- `ctcf_elements_per_guide_crispri_effects.tsv`: per-guide effects at the two elements, with hg38 positions and distances to the CTCF motifs
- `validation_30_selected_pairs.tsv`: the 30 selected pairs with original and validation effect sizes, 95% CIs, p-values and significance calls
