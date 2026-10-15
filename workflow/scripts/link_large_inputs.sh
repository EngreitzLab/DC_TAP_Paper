#!/usr/bin/env bash
# Symlink the large figure inputs listed in config/large_inputs_on_sherlock.txt into this clone.
# Directories are linked file by file, so only files missing from the repository are linked.
# usage: bash workflow/scripts/link_large_inputs.sh <source_dir> [fdr20_run_dir]
#   <source_dir>     a copy of the original analysis directory that contains these files
#   [fdr20_run_dir]  a copy of the final FDR 20% analysis run, for entries marked "fdr20_run"
set -euo pipefail
source_dir="${1:?usage: link_large_inputs.sh <source_dir> [fdr20_run_dir]}"
fdr20_run_dir="${2:-$source_dir}"
grep -vE '^\s*(#|$)' config/large_inputs_on_sherlock.txt | while IFS=$'\t' read -r entry source_key; do
  base="$source_dir"; [ "${source_key:-}" = "fdr20_run" ] && base="$fdr20_run_dir"
  [ -e "$base/$entry" ] || { echo "missing: $base/$entry" >&2; exit 1; }
  (cd "$base" && find "$entry" -type f) | while read -r path; do
    [ -e "$path" ] && continue
    mkdir -p "$(dirname "$path")"
    ln -s "$base/$path" "$path"
    echo "linked: $path"
  done
done
