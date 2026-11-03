#!/usr/bin/env bash
# Full validation-screen pipeline: filtered Cell Ranger matrices -> merged Seurat object -> cell QC
# -> six SCEPTRE objects -> six SCEPTRE runs (union gRNA integration). Takes several hours and
# about 5 GB of disk. Run from validations/ after:
#   bash scripts/bash/download_filtered_h5.sh
#   bash scripts/bash/download_references.sh
# Requires R with the packages in renv.lock and pandoc on PATH.
# Usage: bash scripts/bash/run_validation_pipeline.sh [n_processors]

set -euo pipefail

n_processors="${1:-4}"
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

render() {
  Rscript -e "rmarkdown::render('$1', params = $2, output_dir = '$3', output_file = '$4')"
}

render analyses/2026-05-11_qc_experiment/combine_h5_counts.Rmd "list()" results/qc_experiment_June_2026 combine_h5_counts.html
render analyses/2026-05-11_qc_experiment/qc_experiment.Rmd "list()" results/qc_experiment_June_2026 qc_experiment.html

# SCEPTRE objects: one per sample, and one per MOI with both replicates (sample as a covariate)
objects=(
  "low|c('MOI1_1')|c('channel')|MOI1_1"
  "low|c('MOI1_5')|c('channel')|MOI1_5"
  "high|c('MOI11_3')|c('channel')|MOI11_3"
  "high|c('MOI13')|c('channel')|MOI13"
  "low|c('MOI1_1', 'MOI1_5')|c('channel', 'sample')|low_moi_merged"
  "high|c('MOI11_3', 'MOI13')|c('channel', 'sample')|high_moi_merged"
)

for object in "${objects[@]}"; do
  IFS='|' read -r moi_level sample_list covariates prefix <<< "${object}"
  render analyses/2026-05-11_SCEPTRE/01_create_sceptre_objects.Rmd \
    "list(moi_level = '${moi_level}', sample_list = ${sample_list}, covariate_columns = ${covariates}, prefix = '${prefix}')" \
    results/2026-06-10_SCEPTRE "01_${prefix}.html"
done

for object in "${objects[@]}"; do
  IFS='|' read -r moi_level _ _ prefix <<< "${object}"
  render analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd \
    "list(sceptre_object_path = 'results/2026-06-10_SCEPTRE/objects/${prefix}_sceptre_object_${moi_level}_moi_merged.rds', prefix = '${prefix}', integration_strategy = 'union', n_processors = ${n_processors})" \
    results/2026-06-10_SCEPTRE "02_${prefix}_union.html"
done
