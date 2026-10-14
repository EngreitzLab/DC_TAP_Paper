#!/usr/bin/env bash
# Re-create all data panels of main Figures 1-6 from the tables in the repository (plus the large
# inputs listed in config/large_inputs_on_sherlock.txt). Run inside the pixi environment:
#   pixi run bash workflow/scripts/run_main_figures.sh [cores]
# Only figure rules may run (--allowed-rules): a fresh clone gives every file a new timestamp, which
# would otherwise make Snakemake re-run the full upstream pipeline from the raw CellRanger outputs.
set -euo pipefail
cores="${1:-8}"
figure_rules="distance_by_effect_size_plots effect_size_plots slc2a3_effect_size_plot make_replicates_plots
  compare_k562_effect_sizes_to_qPCR downsampling_reads_for_UMIs find_duplicate_pairs wtc11_multi_moi
  compare_effect_sizes_to_other_screens power_to_detect_change_for_all_datasets
  main_figure_4_direct_indirect_effects main_figure_5_chromatin_categories main_figure_5_housekeeping_genes
  main_figure_6_model_scores"
snakemake --use-conda --conda-frontend conda --cores "$cores" --keep-going \
  --allowed-rules $figure_rules negative_targeting_confidence_interval \
    create_main_figure_1_and_2 create_main_figure_3 create_main_figures_4_5_6 make_main_figures \
  --forcerun $figure_rules \
  make_main_figures
