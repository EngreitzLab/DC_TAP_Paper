#!/usr/bin/env bash

set -euo pipefail

if ! command -v gsutil >/dev/null 2>&1; then
  echo "Error: gsutil is required but was not found in PATH." >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${REPO_ROOT}/data/filtered_h5"
BASE_URI="${VALIDATION_CELLRANGER_URI:?set VALIDATION_CELLRANGER_URI to the gs:// folder holding the Cell Ranger outputs}"

mkdir -p "${OUTPUT_DIR}"

samples=(
  "validation_moi1_1_rep01_a1"
  "validation_moi1_1_rep02_a2"
  "validation_moi1_1_rep03_a3"
  "validation_moi1_1_rep04_a4"
  "validation_moi1_1_rep05_a5"
  "validation_moi1_1_rep06_a6"
  "validation_moi11_3_rep01_a7"
  "validation_moi11_3_rep02_a8"
  "validation_moi1_5_rep01_b1"
  "validation_moi1_5_rep02_b2"
  "validation_moi1_5_rep03_b3"
  "validation_moi1_5_rep04_b4"
  "validation_moi1_5_rep05_b5"
  "validation_moi1_5_rep06_b6"
  "validation_moi13_rep01_b7"
  "validation_moi13_rep02_b8"
)

for sample in "${samples[@]}"; do
  source_uri="${BASE_URI}/${sample}/outs/filtered_feature_bc_matrix.h5"
  destination="${OUTPUT_DIR}/${sample}_filtered_feature_bc_matrix.h5"

  echo "Downloading ${source_uri} -> ${destination}"
  gsutil cp -n "${source_uri}" "${destination}"
done

echo "Downloaded ${#samples[@]} H5 files to ${OUTPUT_DIR}"