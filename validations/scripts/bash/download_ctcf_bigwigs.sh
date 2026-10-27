#!/usr/bin/env bash
# Download the 18 CTCF ChIP-seq bigWig files (CPM-normalized coverage) from GEO (GSE303901).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
LIBRARY_SHEET="${REPO_ROOT}/metadata/ctcf_chip_libraries.tsv"
OUTPUT_DIR="${REPO_ROOT}/data/ctcf_chip_bigwigs"

mkdir -p "${OUTPUT_DIR}"

gsm_col=$(head -1 "${LIBRARY_SHEET}" | tr '\t' '\n' | grep -n -x "geo_gsm" | cut -d: -f1)
file_col=$(head -1 "${LIBRARY_SHEET}" | tr '\t' '\n' | grep -n -x "geo_bigwig" | cut -d: -f1)

while IFS=$'\t' read -r gsm bigwig; do
  destination="${OUTPUT_DIR}/${bigwig}"
  if [[ ! -s "${destination}" ]]; then
    url="https://ftp.ncbi.nlm.nih.gov/geo/samples/${gsm:0:7}nnn/${gsm}/suppl/${bigwig}"
    echo "Downloading ${url}"
    curl -fsSL --retry 3 -o "${destination}.part" "${url}"
    mv "${destination}.part" "${destination}"
  fi
done < <(tail -n +2 "${LIBRARY_SHEET}" | cut -f "${gsm_col},${file_col}")
