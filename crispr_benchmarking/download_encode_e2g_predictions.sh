#!/usr/bin/env bash
# Download the ENCODE-rE2G and ABC predictions used in dc_tap_seq_paper.yml from the ENCODE portal.
# usage: bash crispr_benchmarking/download_encode_e2g_predictions.sh   (run from the repository root)
set -euo pipefail
out_dir="resources/e2g_predictions"
mkdir -p "$out_dir"
for accession in ENCFF970QAX ENCFF071ZZU ENCFF681DDZ ENCFF660VQW; do
  [ -s "$out_dir/$accession.bed.gz" ] || \
    curl -fL --retry 3 -o "$out_dir/$accession.bed.gz" "https://www.encodeproject.org/files/$accession/@@download/$accession.bed.gz"
done
