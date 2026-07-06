#!/usr/bin/env python3
"""
Download ENCODE narrowPeak files selected for the DC-TAP-seq ENCODE overlap analysis.

Input : metadata/encode_selected_files.csv  (produced by the curation step; one row per
        kept, unperturbed, released GRCh38 ChIP-seq experiment with a chosen narrowPeak file)
Output: data/encode_peaks/<cellset>/<file_accession>.bed.gz

Usage:
    python scripts/download_encode.py \
        --manifest metadata/encode_selected_files.csv \
        --outdir data/encode_peaks

Features:
    * Resumable: files already present with the correct md5 are skipped.
    * Verifies md5sum against ENCODE metadata; re-downloads on mismatch.
    * No credentials required - all files are public.

Dependencies: python>=3.8 (standard library only)
"""
import argparse, hashlib, os, sys, time
import csv
import urllib.request

def md5sum(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

def download(url, dest, expected_md5=None, retries=3):
    for attempt in range(1, retries + 1):
        try:
            urllib.request.urlretrieve(url, dest)
            if expected_md5 and md5sum(dest) != expected_md5:
                raise ValueError("md5 mismatch")
            return True
        except Exception as e:
            if attempt == retries:
                sys.stderr.write(f"  FAILED {url}: {e}\n")
                if os.path.exists(dest):
                    os.remove(dest)
                return False
            time.sleep(2 * attempt)
    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="metadata/encode_selected_files.csv")
    ap.add_argument("--outdir", default="data/encode_peaks")
    args = ap.parse_args()

    with open(args.manifest) as f:
        rows = list(csv.DictReader(f))

    n_ok = n_skip = n_fail = 0
    for i, r in enumerate(rows, 1):
        cellset = r["cellset"]
        facc = r["file"]
        url = r["href"]
        md5 = r.get("md5sum") or None
        if not facc or not url:
            continue
        d = os.path.join(args.outdir, cellset)
        os.makedirs(d, exist_ok=True)
        dest = os.path.join(d, f"{facc}.bed.gz")
        if os.path.exists(dest) and (md5 is None or md5sum(dest) == md5):
            n_skip += 1
            continue
        ok = download(url, dest, md5)
        if ok:
            n_ok += 1
        else:
            n_fail += 1
        if i % 50 == 0:
            print(f"[{i}/{len(rows)}] downloaded={n_ok} skipped={n_skip} failed={n_fail}")

    print(f"DONE: downloaded={n_ok} skipped={n_skip} failed={n_fail} total={len(rows)}")

if __name__ == "__main__":
    main()
