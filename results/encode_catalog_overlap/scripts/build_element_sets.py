#!/usr/bin/env python3
"""
build_element_sets.py

Step 1 of the DC-TAP-seq x ENCODE overlap pipeline.

Parses Supplementary Table S3 and defines, per cell type (K562, WTC11),
the element-locus sets used downstream:

  pos_neg     positive elements with >=1 significant NEGATIVE-effect pair
              (repression lowers target expression -> activating enhancer)
  pos_pos     positive elements with >=1 significant POSITIVE-effect pair
  bg_powered  background: never significant AND well powered
              (power_at_effect_size_15 > 0.8)   [primary background]
  bg_all      background: never significant (any power)  [sensitivity]

Only Random_DistalElement_Gene == True pairs are used. Elements are the
resized/merged hg38 intervals (resized_merged_targeting_*_hg38), which are
the coordinates the manuscript uses for chromatin overlap. Pairs are
collapsed to unique element loci (elem_key = chr:start-end).

Outputs (written under --outdir, default ./elements):
  {CELL}_{set}.bed        BED4 (chr, start, end, elem_key) per set
  element_manifest.csv     one row per (element, set) with gene / power info

Usage:
  python build_element_sets.py --table-s3 path/to/Table_S3.tsv
"""
import argparse
import os
import numpy as np
import pandas as pd

# resized/merged hg38 element coordinates (used for all chromatin overlap)
CHR = "resized_merged_targeting_chr_hg38"
ST = "resized_merged_targeting_start_hg38"
EN = "resized_merged_targeting_end_hg38"
# element-gene distance (bp) from element to the gene's GENCODE TSS
DIST = "distance_to_gencode_gene_TSS"

CELLS = ["K562", "WTC11"]


def build(df, cell):
    """Collapse Random distal pairs for one cell type to unique element loci."""
    s = df[(df["cell_type"] == cell) & (df["Random_DistalElement_Gene"] == True)].copy()
    s["elem_key"] = s[CHR].astype(str) + ":" + s[ST].astype(str) + "-" + s[EN].astype(str)

    def agg(g):
        sig = g[g["significant"] == True]
        return pd.Series({
            "chr": g[CHR].iloc[0], "start": g[ST].iloc[0], "end": g[EN].iloc[0],
            "n_pairs": len(g), "n_sig": len(sig),
            "any_sig": len(sig) > 0,
            "any_sig_neg": (sig["log_2_FC_effect_size"] < 0).any(),
            "any_sig_pos": (sig["log_2_FC_effect_size"] > 0).any(),
            "sig_genes": ",".join(sorted(set(sig["gene_symbol"]))) if len(sig) else "",
            "all_genes": ",".join(sorted(set(g["gene_symbol"]))),
            "element_category": g["element_category"].iloc[0],
            "max_power15": g["power_at_effect_size_15"].max(),
            "min_sig_log2FC": sig["log_2_FC_effect_size"].min() if len(sig) else np.nan,
            # minimum element-to-TSS distance (bp) across pairs; sig-only and all
            "min_dist_tss": g[DIST].abs().min(),
            "min_sig_dist_tss": sig[DIST].abs().min() if len(sig) else np.nan,
        })

    e = s.groupby("elem_key").apply(agg, include_groups=False).reset_index()
    for c in ["start", "end", "n_pairs", "n_sig"]:
        e[c] = e[c].astype(int)
    return e


def write_bed(e_sub, path):
    (e_sub.sort_values(["chr", "start"])[["chr", "start", "end", "elem_key"]]
        .to_csv(path, sep="\t", header=False, index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table-s3", required=True, help="Path to Table S3 TSV")
    ap.add_argument("--outdir", default="elements")
    args = ap.parse_args()

    df = pd.read_csv(args.table_s3, sep="\t")
    os.makedirs(args.outdir, exist_ok=True)

    manifest = []
    for cell in CELLS:
        e = build(df, cell)
        groups = {
            "pos_neg": e[e["any_sig_neg"]],
            "pos_pos": e[e["any_sig_pos"]],
            "bg_powered": e[(~e["any_sig"]) & (e["max_power15"] > 0.8)],
            "bg_all": e[~e["any_sig"]],
        }
        for name, sub in groups.items():
            write_bed(sub, os.path.join(args.outdir, f"{cell}_{name}.bed"))
            tmp = sub.copy()
            tmp["cell"] = cell
            tmp["set"] = name
            manifest.append(tmp)
            print(f"{cell:6s} {name:11s} n={len(sub)}")

    man = pd.concat(manifest, ignore_index=True)
    man.to_csv(os.path.join(args.outdir, "element_manifest.csv"), index=False)
    print(f"\nWrote {args.outdir}/element_manifest.csv  ({len(man)} element-set rows)")


if __name__ == "__main__":
    main()
