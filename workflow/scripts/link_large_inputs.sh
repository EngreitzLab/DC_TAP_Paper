#!/usr/bin/env bash
# Symlink the large figure inputs listed in config/large_inputs_on_sherlock.txt into this clone.
# Directories are linked file by file, so only files missing from the repository are linked.
# usage: bash workflow/scripts/link_large_inputs.sh <source_dir>
#   <source_dir> = a copy of the original analysis directory that contains these files
set -euo pipefail
source_dir="${1:?usage: link_large_inputs.sh <source_dir>}"
grep -vE '^\s*(#|$)' config/large_inputs_on_sherlock.txt | while read -r entry; do
  [ -e "$source_dir/$entry" ] || { echo "missing: $source_dir/$entry" >&2; exit 1; }
  (cd "$source_dir" && find "$entry" -type f) | while read -r path; do
    [ -e "$path" ] && continue
    mkdir -p "$(dirname "$path")"
    ln -s "$source_dir/$path" "$path"
    echo "linked: $path"
  done
done
