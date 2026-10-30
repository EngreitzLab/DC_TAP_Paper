#!/usr/bin/env bash

BASE_URL="https://zenodo.org/api/records/23775358/draft/files"
OUT_DIR="data"

FILES=(
    "ENCODE_Gasperini2019_0.13gStd_MAST_perCRE_GRCh38.tsv.gz"
    "Gasperini2019_output_MAST_perCRE.tsv.gz"
    "Gasperini2019_output_trans_effects_MAST_perCRE.tsv.gz"
)

mkdir -p "$OUT_DIR"

for f in "${FILES[@]}"; do
    echo "Downloading $f"
    curl -fSL -H "Authorization: Bearer $ZENODO_TOKEN" \
        -o "$OUT_DIR/$f" \
        "$BASE_URL/$f/content"
done

echo "Done. Files saved in $OUT_DIR/"
