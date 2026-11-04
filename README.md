# DC-TAP-seq paper: analysis code

Code to analyze the K562 and WTC11 DC-TAP-seq CRISPRi screens and to make the figures and tables of:

> Ray et al. *A systematic survey of distal element-gene regulatory interactions with Direct-Capture Targeted Perturb-seq.* [bioRxiv 10.1101/2025.09.16.676677](https://www.biorxiv.org/content/10.1101/2025.09.16.676677v2)

Release used for the paper: `v2.0.1` (Zenodo [10.5281/zenodo.23964909](https://doi.org/10.5281/zenodo.23964909)).

## Data

| Data | Where |
|---|---|
| Raw sequencing and Cell Ranger outputs (pilot, screens, validation screens, CTCF ChIP-seq) | GEO [GSE303901](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303901) |
| Large intermediate files used by the figures | Zenodo [10.5281/zenodo.23775358](https://doi.org/10.5281/zenodo.23775358) (list: `config/large_inputs.txt`) |
| gRNA designs | IGVF Data Portal `IGVFFI7272YLVI` (K562), `IGVFFI0580WJFK` (WTC11) |
| ENCODE annotations | ENCODE Portal; accessions in Supplementary Table 8 |
| Final element-gene results (Table S3) | `results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv`; columns described in `resources/formatting_dc_tap_results/columns_of_Final_DC_TAP_Seq_Results_table.txt` |

## Repository layout

| Path | Contents |
|---|---|
| `workflow/` | Snakemake pipeline: `Snakefile`, `rules/*.smk`, `scripts/`, per-rule conda envs in `envs/` |
| `resources/` | Pipeline inputs (design files, published screens, qPCR tables) |
| `results/` | Pipeline outputs and figure scripts, one folder per figure |
| `direct_vs_indirect_effects/` | Workflow estimating indirect effect rates and the probability of direct effects (Fig. 4, ED Figs. 7, 15, 16) |
| `validations/` | Validation screen and CTCF ChIP-seq analyses |
| `crispr_benchmarking/` | Inputs for benchmarking E2G models on DC-TAP-seq pairs (Fig. 6) |
| `seqspec/` | seqspec templates for DC-TAP-seq libraries |
| `config/` | Pipeline config and list of Zenodo inputs |

Pipeline steps (`workflow/rules/`):

1. `sceptre_setup.smk`: Cell Ranger outputs → SCEPTRE input objects.
2. `sceptre_power_analysis.smk`: SCEPTRE differential expression and power simulations.
3. `create_encode_output.smk`: element annotation (promoter, gene body), liftover, ENCODE CRISPR format.
4. `formatting_dc_tap_results.smk`: confidence intervals, element-gene categories, chromatin categories → Table S3.
5. `main_figure_*.smk`, `figures_wo_positive_controls.smk`, `supplementary_tables.smk`: figures and tables.

## System requirements

- OS: Linux x86_64 (tested: CentOS 7, Stanford Sherlock cluster) or macOS arm64 (tested: macOS 26.6, figure rules only).
- Software, installed by pixi from `pixi.toml` / `pixi.lock`: Snakemake 7.32.4, Python 3.11, R 4.4.3 with tidyverse, data.table, cowplot, ggExtra, rmarkdown.
- Per-rule software (SCEPTRE, Bioconductor, UCSC liftOver), installed by conda from `workflow/envs/*.yml`. SCEPTRE comes from the `jgalante` Anaconda channel, pinned by version, build and sha256:
  - `r-sceptre 0.10.1`: [Katsevich-Lab/sceptre](https://github.com/Katsevich-Lab/sceptre) commit `43a4534` (tag `v0.10.1`)
  - `r-sceptre-dev 0.10.1dev` (confidence intervals): commit `1517ebc`
  - Alternative install: `remotes::install_github("Katsevich-Lab/sceptre@<commit>")`
- Hardware: figure rules run on a laptop. The full pipeline needs a compute cluster: rules request up to 72 GB RAM and 12 h.

## Installation

```
git clone https://github.com/EngreitzLab/DC_TAP_Paper.git
cd DC_TAP_Paper
git checkout v2.0.1
pixi install        # https://pixi.sh; ~1 min
```

Per-rule conda environments are built on first use (tens of minutes).

## Demo

Make Fig. 4b–d and Extended Data Figs. 7, 15 and 16b,c from files in the repository:

```
pixi run snakemake --cores 4 \
  results/main_figure_4_and_extdata_figure_7/dc_tapseq_indirect_effects_main_analyses.html \
  results/extended_data_figure_15/dc_tapseq_indirect_effects_extdata_fig15.html \
  results/extended_data_figure_16/3_genome_wide_indirect_effects_extdata_fig16bc.html
```

Output: PDFs in `results/main_figure_4_and_extdata_figure_7/plots/`, `results/extended_data_figure_15/` and `results/extended_data_figure_16/`, plus an HTML report for each. Run time: ~20 s on a laptop (Apple M-series, 4 cores).

## Reproduce the figures

```
bash workflow/scripts/download_large_inputs_from_zenodo.sh   # ~0.5 GB from Zenodo, MD5-checked
pixi run main-figures                                         # ~90 min on 8 cores, incl. conda env builds
```

Figures 4–6 and Extended Data Figs. 7, 15, 16b,c only: `pixi run snakemake --cores 4 create_main_figures_4_5_6`.

| Panel | Rule | Script | Output |
|---|---|---|---|
| 1c, 1d | `main_figure_1cd_guide_capture_and_moi` | `workflow/scripts/main_figure_1_and_2/plot_guide_capture_and_moi_knockdown.R` | `results/main_figure_1_and_2/guide_capture_and_moi_knockdown/` |
| 1e | `main_figure_1e_umi_saturation` | `workflow/scripts/main_figure_1_and_2/plot_umi_saturation_dctap_vs_scrna.R` | `results/main_figure_1_and_2/umi_saturation_dctap_vs_scrna/` |
| 2b | `distance_by_effect_size_plots_FDR_wo_pos_controls` | `workflow/scripts/figures_wo_positive_controls/distance_by_effect_size_plots.R` | `results/figures_wo_positive_controls/{k562,wtc11}_distance_by_es.pdf` |
| 2c | `effect_size_plots_FDR_wo_pos_controls` | `workflow/scripts/figures_wo_positive_controls/effect_size_plots.R` | `results/figures_wo_positive_controls/effect_size_boxplot.pdf` |
| 2d | `wtc11_multi_moi`, `make_replicates_plots`, `main_figure_2d_validation_vs_screen` | `workflow/scripts/main_figure_1_and_2/`, `validations/analyses/2026-05-14_compare_low_high_es/compare_es.Rmd` | `results/main_figure_1_and_2/`, `validations/results/2026-05-14_compare_low_high_es/` |
| 2e | `main_figure_2e_dctap_vs_qpcr` | `workflow/scripts/main_figure_1_and_2/plot_dctap_vs_qpcr_effect_sizes.R` | `results/main_figure_1_and_2/qpcr_validation/` |
| 2f | `find_duplicate_pairs` | `workflow/scripts/main_figure_1_and_2/find_duplicate_pairs.R` | `results/main_figure_1_and_2/duplicate_pairs_analysis/duplicate_pairs.pdf` |
| 2h | `slc2a3_effect_size_plot` | `workflow/scripts/main_figure_1_and_2/slc2a3_effect_size_plot.R` | `results/main_figure_1_and_2/slc2a3_effect_size_barplot.pdf` |
| 2i | `main_figure_2i_slc2a3_qpcr` | `workflow/scripts/main_figure_1_and_2/plot_slc2a3_qpcr.R` | `results/main_figure_1_and_2/qpcr_validation/` |
| 3a, 3b | `compare_effect_sizes_to_other_screens_FDR_wo_pos_controls`, `power_to_detect_change_for_all_datasets_FDR_wo_pos_controls` | `workflow/scripts/figures_wo_positive_controls/` | `results/figures_wo_positive_controls/` |
| 3c, 5b–d | `main_figure_5_chromatin_categories` | `results/main_figure_5_and_6/chromatin_categories_and_power_estimation.R` | `results/main_figure_5_and_6/` (5c, 5d: `housekeeping_genes/fig5{c,d}_*.pdf`) |
| 4b–d, ED 7 | `main_figure_4_direct_indirect_effects` | `results/main_figure_4_and_extdata_figure_7/dc_tapseq_indirect_effects_main_analyses.Rmd` | `results/main_figure_4_and_extdata_figure_7/plots/` |
| 6a–d, ED 13 | `main_figure_6_model_scores` | `results/main_figure_5_and_6/model_scores.R` | `results/main_figure_5_and_6/e2g_model_scores/` |
| ED 15 | `extended_data_figure_15_direct_indirect_rates` | `results/extended_data_figure_15/dc_tapseq_indirect_effects_extdata_fig15.Rmd` | `results/extended_data_figure_15/` |
| ED 16 | `extended_data_figure_16bc_genome_wide_indirect_effects` (b,c) | `results/extended_data_figure_16/` (steps in its README) | `results/extended_data_figure_16/` |

Derived inputs to Figs. 4–6:

- `direct_vs_indirect_effects/compute_direct_vs_indirect_effects/results/`: cis/trans hit rates and direct-effect models (made by that folder's Snakemake workflow).
- `resources/main_figure_6/expt_pred_merged_annot.txt.gz`: output of [CRISPR_comparison](https://github.com/EngreitzLab/CRISPR_comparison) for DC-TAP-seq pairs.
- `results/main_figure_5_and_6/categorized_pairs/`: chromatin-categorized CRISPR pairs (`chromatin_categories_and_power_estimation.R categorize`).

## Run the pipeline from Cell Ranger outputs

Download the Cell Ranger outputs from GEO GSE303901 into `resources/process_validation_datasets/{K562,WTC11}_DC_TAP_Seq/cell_ranger_output/`, then:

```
snakemake --use-conda --cores <N>
```

This makes the final results table (`results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements.tsv`), the supplementary tables and Figs. 1–2.

To analyze a new DC-TAP-seq screen, build its SCEPTRE inputs (as in `sceptre_setup.smk`) in `results/process_validation_datasets/<sample>/`: `differential_expression/sceptre_diffex_input.rds`, `gene_gRNA_group_pairs.rds`, `gRNA_groups_table.rds`, `guide_targets.tsv`, `distances.tsv`. Then run differential expression and power simulations:

```
snakemake --use-conda --cores 4 results/process_validation_datasets/<sample>/differential_expression/results_run_discovery_analysis.rds
snakemake --use-conda --cores 4 results/process_validation_datasets/<sample>/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv
```

## License

MIT; see [`LICENSE`](LICENSE).
