#!/usr/bin/env bash
set -euo pipefail

gtf_filepath="${1:-}"
if [[ -z "$gtf_filepath" || "$gtf_filepath" != *.gtf.gz ]]; then
  echo "Usage: $0 /path/to/file.gtf.gz" >&2
  exit 1
fi

out="${gtf_filepath%.gtf.gz}.exon.bed"

zcat "$gtf_filepath" |
  awk 'OFS="\t" {
  if ($1 ~ /^chr/ && $3=="exon") {
    if ($7 == "+") { print $1,$4-1,$4,$16,".",$7,".",".",".",".",".",".",$10,$12}
    else           { print $1,$5-1,$5,$16,".",$7,".",".",".",".",".",".",$10,$12}
    }
  }' |
  tr -d '";' |
  sort -k1,1V -k2,2n >"$out"

echo "$out"
