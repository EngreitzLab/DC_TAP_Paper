#!/usr/bin/env bash
# Download the large figure inputs (config/large_inputs_on_sherlock.txt) from Zenodo into this clone.
# Files already in the repository are kept (tar -k); the per-file checksums are then verified.
# usage: bash workflow/scripts/download_large_inputs_from_zenodo.sh   (run from the repository root)
set -euo pipefail
zenodo_record="23775358"
url="https://zenodo.org/records/${zenodo_record}/files"
download_dir="resources/zenodo_download"
mkdir -p "$download_dir"
for f in dc_tap_paper_large_inputs.tar dc_tap_paper_large_inputs.md5; do
  [ -s "$download_dir/$f" ] || curl -fL --retry 3 -o "$download_dir/$f" "$url/$f?download=1"
done
tar -xkf "$download_dir/dc_tap_paper_large_inputs.tar" 2>/dev/null || true
md5_tool=$(command -v md5sum || echo "md5 -r")
grep -vE '^\s*$' "$download_dir/dc_tap_paper_large_inputs.md5" | while read -r sum path; do
  [ "$($md5_tool "$path" | cut -d' ' -f1)" = "$sum" ] || { echo "checksum mismatch: $path" >&2; exit 1; }
done
echo "large inputs in place"
