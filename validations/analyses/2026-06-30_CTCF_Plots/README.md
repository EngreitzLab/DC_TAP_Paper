# 2026-06-30 CTCF Plots

## Purpose

Figure panels for the FAM111A and CDC37 CTCF elements:

- CTCF ChIP-seq signal per sgRNA line (boxplots, one point per ChIP replicate);
- per-guide CRISPRi effects of the WTC11 screen on FAM111A and CDC37, against the guide's distance from the CTCF motif.

Statistical tests of the ChIP-seq signal are in `analyses/2026-10-23_ctcf_and_validation_statistics/01_ctcf_chip_statistics.Rmd`.

## Inputs

- `analyses/2026-06-29_CTCF_ChIP_validation/source_data/ctcf_signal_per_replicate.tsv`
- Singleton (per-guide) SCEPTRE results and guide positions (hg19) of the WTC11 screen: `../results/process_validation_datasets/WTC11_DC_TAP_Seq/singleton_differential_expression/results_run_discovery_analysis.rds` and `../results/process_validation_datasets/WTC11_DC_TAP_Seq/guide_targets.tsv`
- `data/hg19ToHg38.over.chain` (from `scripts/bash/download_references.sh`)
- `ctcf_element_motifs_hg38.bed`: element windows and the CTCF motifs marked on the per-guide plots
- `chip_sgrnas_hg38.bed`: the four sgRNAs used for the ChIP-seq lines (genome-browser track)

## Code

Run from `validations/`:

```bash
bash scripts/bash/download_references.sh
Rscript -e "rmarkdown::render('analyses/2026-06-30_CTCF_Plots/plots_for_paper.Rmd', output_dir = 'results/2026-06-30_CTCF_Plots')"
```

## Outputs

In `results/2026-06-30_CTCF_Plots/`:

- `fam111a_boxplots.pdf`, `cdc37_boxplots.pdf`
- `fam111a_per_guide_effect_relative_hg38.pdf`, `cdc37_per_guide_effect_relative_hg38.pdf`
