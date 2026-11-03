#!/usr/bin/env bash
# Download the large figure inputs (config/large_inputs.txt) from Zenodo into this clone.
# Each Zenodo file is unpacked to its repository path as listed in file_manifest.tsv, then the
# per-file MD5 sums in file_checksums.md5 are verified. Files already downloaded are not fetched again.
# usage: bash workflow/scripts/download_large_inputs_from_zenodo.sh   (run from the repository root)
set -euo pipefail
zenodo_record="23775358"
url="https://zenodo.org/records/${zenodo_record}/files"
download_dir="resources/zenodo_download"
mkdir -p "$download_dir"
fetch() { [ -s "$download_dir/$1" ] || curl -fL --retry 3 -o "$download_dir/$1" "$url/$1?download=1"; }
fetch file_manifest.tsv
fetch file_checksums.md5
tail -n +2 "$download_dir/file_manifest.tsv" | while IFS=$'\t' read -r zenodo_file repo_path unpack description; do
  fetch "$zenodo_file"
  case "$unpack" in
    untar)  tar -xf "$download_dir/$zenodo_file" ;;
    copy)   mkdir -p "$(dirname "$repo_path")" && cp "$download_dir/$zenodo_file" "$repo_path" ;;
    gunzip) mkdir -p "$(dirname "$repo_path")" && gzip -dc "$download_dir/$zenodo_file" > "$repo_path" ;;
  esac
  echo "placed: $repo_path"
done
md5_tool=$(command -v md5sum || echo "md5 -r")
while read -r sum path; do
  [ "$($md5_tool "$path" | cut -d' ' -f1)" = "$sum" ] || { echo "checksum mismatch: $path" >&2; exit 1; }
done < "$download_dir/file_checksums.md5"
echo "large inputs in place"
