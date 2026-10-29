#!/usr/bin/env bash
# Download the public reference files used by the validation analyses into data/.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${REPO_ROOT}/data"

mkdir -p "${DATA_DIR}"

download() {
  local url="$1" destination="$2"
  if [[ ! -s "${destination}" ]]; then
    echo "Downloading ${url} -> ${destination}"
    curl -fsSL --retry 3 -o "${destination}.part" "${url}"
    mv "${destination}.part" "${destination}"
  fi
}

# GENCODE v32 gene annotation (gene symbols and biotypes)
download "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz" \
  "${DATA_DIR}/gencode.v32.primary_assembly.annotation.gtf.gz"

# Guide feature reference of the WTC11 screen (Cell Ranger feature reference)
download "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303901/suppl/GSE303901_wtc11_feature_reference.csv.gz" \
  "${DATA_DIR}/GSE303901_wtc11_feature_reference.csv.gz"
if [[ ! -s "${DATA_DIR}/GSE303901_wtc11_feature_reference.csv" ]]; then
  gunzip -c "${DATA_DIR}/GSE303901_wtc11_feature_reference.csv.gz" > "${DATA_DIR}/GSE303901_wtc11_feature_reference.csv"
fi

# UCSC hg19 -> hg38 liftOver chain (sgRNA coordinates of the screen are in hg19)
download "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" \
  "${DATA_DIR}/hg19ToHg38.over.chain.gz"
if [[ ! -s "${DATA_DIR}/hg19ToHg38.over.chain" ]]; then
  gunzip -c "${DATA_DIR}/hg19ToHg38.over.chain.gz" > "${DATA_DIR}/hg19ToHg38.over.chain"
fi
