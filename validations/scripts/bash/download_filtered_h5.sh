#!/usr/bin/env bash
# Download the 16 filtered Cell Ranger matrices of the validation screen from GEO (GSE303901).
# GEO names the high-MOI samples by their measured MOI (moi12_6, moi14_4); the files are saved
# under the sample IDs used in the code (moi11_3, moi13). See metadata/validation_screen_samples.tsv.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SAMPLE_SHEET="${REPO_ROOT}/metadata/validation_screen_samples.tsv"
OUTPUT_DIR="${REPO_ROOT}/data/filtered_h5"
GEO_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303901/suppl"

mkdir -p "${OUTPUT_DIR}"

sample_col=$(head -1 "${SAMPLE_SHEET}" | tr '\t' '\n' | grep -n -x "sample_id" | cut -d: -f1)
file_col=$(head -1 "${SAMPLE_SHEET}" | tr '\t' '\n' | grep -n -x "geo_h5_file" | cut -d: -f1)

n=0
while IFS=$'\t' read -r sample geo_file; do
  destination="${OUTPUT_DIR}/${sample}_filtered_feature_bc_matrix.h5"
  if [[ ! -s "${destination}" ]]; then
    echo "Downloading ${geo_file} -> ${destination}"
    curl -fsSL --retry 3 -o "${destination}.part" "${GEO_URL}/${geo_file}"
    mv "${destination}.part" "${destination}"
  fi
  n=$((n + 1))
done < <(tail -n +2 "${SAMPLE_SHEET}" | cut -f "${sample_col},${file_col}")

echo "Downloaded ${n} H5 files to ${OUTPUT_DIR}"
