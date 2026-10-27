#!/usr/bin/env bash
# Figures and tables of the validation-screen and CTCF analyses from the committed source data
# (source_data/ folders); takes a few minutes. Run from validations/ after:
#   bash scripts/bash/download_references.sh
# Requires R with the packages in renv.lock and pandoc on PATH.
# Usage: bash scripts/bash/render_validation_figures.sh

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

render() {
  Rscript -e "rmarkdown::render('$1', output_dir = '$2')"
}

Rscript analyses/2026-05-11_cellranger/extract_metrics_summary.R
Rscript -e "rmarkdown::render('analyses/2026-05-11_guide_pool_qc/pool_QC.Rmd', output_dir = 'results/guide_pool_qc', output_file = 'guide_pool_qc.pdf')"
render analyses/2026-05-13_get_statistics/results_statistics.Rmd results/2026-05-13_get_statistics
render analyses/2026-05-14_compare_low_high_es/compare_es.Rmd results/2026-05-14_compare_low_high_es
render analyses/2026-06-30_CTCF_Plots/plots_for_paper.Rmd results/2026-06-30_CTCF_Plots
render analyses/2026-10-23_ctcf_and_validation_statistics/01_ctcf_chip_statistics.Rmd results/2026-10-23_ctcf_and_validation_statistics
render analyses/2026-10-23_ctcf_and_validation_statistics/02_validation_30_pairs_table.Rmd results/2026-10-23_ctcf_and_validation_statistics
