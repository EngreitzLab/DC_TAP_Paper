#!/usr/bin/env python3
"""
compute_overlap_enrichment.py

Step 4-6 of the DC-TAP-seq x ENCODE overlap pipeline.

Requires (produced by earlier steps):
  elements/element_manifest.csv          (build_element_sets.py)
  metadata/encode_selected_files.csv     (curate_encode_metadata.py)
  data/encode_peaks/{K562,iPSC_ESC}/<ENCFF...>.bed.gz  (download_encode.py)

Computes, per cell type:
  1. Binary element x ChIP-seq-target overlap matrices (union of peaks per
     target). Histone marks are extended +/-175 bp before overlap (matching
     the manuscript); TF / point-source peaks are used unextended.
  2. Per-target enrichment of overlap in positives vs background, via a
     two-sided Fisher's exact test with Benjamini-Hochberg FDR, computed
     separately for negative- and positive-effect positives and for both
     backgrounds (well-powered = primary, all-non-significant = sensitivity).
     Background always excludes any element that is itself ever a positive.

Outputs:
  overlap/overlap_matrix_{CELL}.csv    elements x targets (with flag columns)
  overlap/overlap_long_{CELL}.csv      element-experiment-target-assay long form
  overlap/target_assay_{CELL}.csv      target -> assay (TF/Histone) map
  enrichment/enrichment_{CELL}.csv     all sign x background enrichment tables

Usage:
  python compute_overlap_enrichment.py
"""
import os
import warnings
import numpy as np
import pandas as pd
import pyranges as pr
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

CELLS = ["K562", "WTC11"]
CELLSET = {"K562": "K562", "WTC11": "iPSC_ESC"}  # peak-directory / metadata key
HIST_EXT = 175  # bp; histone-mark peak extension each side (manuscript rule)
# ENCODE serves histone modifications under two assay titles: conventional
# "Histone ChIP-seq" and low-input multiplexed "Mint-ChIP-seq". Both are histone
# ChIP and are treated identically here (extended +/-HIST_EXT, grouped as
# histone). All native WTC11 histone marks are Mint-ChIP-seq.
HIST_ASSAYS = {"Histone ChIP-seq", "Mint-ChIP-seq"}
HIST_CLASS = "Histone ChIP-seq"  # normalized class label stored for grouping


def elements_for(man, cell):
    m = man[man["cell"] == cell]
    e = (m.groupby("elem_key")
           .agg(chr=("chr", "first"), start=("start", "first"), end=("end", "first"))
           .reset_index())
    for flag, setname in [("is_pos_neg", "pos_neg"), ("is_pos_pos", "pos_pos"),
                          ("is_bg_powered", "bg_powered"), ("is_bg_all", "bg_all")]:
        e[flag] = e["elem_key"].isin(m[m["set"] == setname]["elem_key"])
    return e


def build_overlap(elem, sel, cell):
    cs = CELLSET[cell]
    e = elem.copy()
    gr_e = pr.PyRanges(pd.DataFrame({
        "Chromosome": e["chr"], "Start": e["start"].astype(int),
        "End": e["end"].astype(int), "elem_key": e["elem_key"]}))
    files = sel[sel["cellset"] == cs].reset_index(drop=True)
    long_rows = []
    missing = 0
    for _, fr in files.iterrows():
        path = f"data/encode_peaks/{cs}/{fr['file']}.bed.gz"
        if not os.path.exists(path):
            missing += 1
            continue
        pk = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                         names=["Chromosome", "Start", "End"], dtype={"Chromosome": str})
        is_hist = fr["assay"] in HIST_ASSAYS
        if is_hist:
            pk["Start"] = (pk["Start"] - HIST_EXT).clip(lower=0)
            pk["End"] = pk["End"] + HIST_EXT
        ov = gr_e.join(pr.PyRanges(pk), how=None)
        if len(ov) == 0:
            continue
        assay_class = HIST_CLASS if is_hist else fr["assay"]
        for k in set(ov.df["elem_key"].unique()):
            long_rows.append((k, fr["experiment"], fr["target"], assay_class))
    if missing:
        print(f"  [{cell}] warning: {missing} peak files not found on disk")
    long = pd.DataFrame(long_rows, columns=["elem_key", "experiment", "target", "assay"])
    return e, long


def enrich(full, mat, tgt_assay, pos_flag, bg_flag):
    pos_keys = full.index[full[pos_flag]]
    # background excludes any element that is ever a positive
    bg_keys = full.index[full[bg_flag] & ~full["is_pos_neg"] & ~full["is_pos_pos"]]
    P, B = len(pos_keys), len(bg_keys)
    rows = []
    for tgt in mat.columns:
        a = int(mat.loc[pos_keys, tgt].sum()); b = P - a
        c = int(mat.loc[bg_keys, tgt].sum()); d = B - c
        if a == 0 and c == 0:
            continue
        OR, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        rows.append({"target": tgt, "assay": tgt_assay.get(tgt, "?"),
                     "pos_overlap": a, "pos_total": P, "pos_frac": a / P if P else np.nan,
                     "bg_overlap": c, "bg_total": B, "bg_frac": c / B if B else np.nan,
                     "odds_ratio": OR, "p_value": p})
    r = pd.DataFrame(rows)
    if len(r):
        r["FDR"] = multipletests(r["p_value"], method="fdr_bh")[1]
        # log2 odds ratio with Haldane-Anscombe 0.5 correction
        a = r["pos_overlap"]; b = r["pos_total"] - r["pos_overlap"]
        c = r["bg_overlap"]; d = r["bg_total"] - r["bg_overlap"]
        r["log2_OR"] = np.log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
        r = r.sort_values("FDR")
    return r


def main():
    man = pd.read_csv("elements/element_manifest.csv")
    sel = pd.read_csv("metadata/encode_selected_files.csv")
    os.makedirs("overlap", exist_ok=True)
    os.makedirs("enrichment", exist_ok=True)

    matrices = {}
    for cell in CELLS:
        elem = elements_for(man, cell)
        e, long = build_overlap(elem, sel, cell)
        tgt_long = long[["elem_key", "target", "assay"]].drop_duplicates()
        mat = pd.crosstab(tgt_long["elem_key"], tgt_long["target"])
        mat = (mat > 0).astype(int).reindex(e["elem_key"].values, fill_value=0)
        meta = e.set_index("elem_key")[["chr", "start", "end",
                                        "is_pos_neg", "is_pos_pos",
                                        "is_bg_powered", "is_bg_all"]]
        full = meta.join(mat)
        tgt_assay = tgt_long.drop_duplicates("target").set_index("target")["assay"]

        full.to_csv(f"overlap/overlap_matrix_{cell}.csv")
        long.to_csv(f"overlap/overlap_long_{cell}.csv", index=False)
        tgt_assay.to_csv(f"overlap/target_assay_{cell}.csv")
        matrices[cell] = (full, mat, tgt_assay)
        print(f"{cell}: {mat.shape[0]} elements x {mat.shape[1]} targets")

    enr_by_cell = {}
    for cell in CELLS:
        full, mat, tgt_assay = matrices[cell]
        parts = []
        for sign, pf in [("neg", "is_pos_neg"), ("pos", "is_pos_pos")]:
            for bgname, bf in [("powered", "is_bg_powered"), ("all", "is_bg_all")]:
                r = enrich(full, mat, tgt_assay, pf, bf).copy()
                r.insert(0, "effect_sign", sign)
                r.insert(1, "background", bgname)
                parts.append(r)
                n_sig = int((r["FDR"] < 0.1).sum()) if len(r) else 0
                print(f"  {cell} {sign:3s} vs {bgname:7s}: {n_sig} targets FDR<0.1")
        etab = pd.concat(parts, ignore_index=True)
        etab.to_csv(f"enrichment/enrichment_{cell}.csv", index=False)
        enr_by_cell[cell] = etab

    write_consolidated_tables(man, matrices, enr_by_cell)


# ---- two consolidated supplementary tables -----------------------------------
# annotation columns pulled from the element manifest (Table S3-derived)
ANN_COLS = ["element_category", "sig_genes", "all_genes", "n_pairs", "n_sig",
            "any_sig_neg", "any_sig_pos", "max_power15", "min_sig_log2FC",
            "min_dist_tss", "min_sig_dist_tss"]
META_COLS = ["elem_key", "chr", "start", "end",
             "is_pos_neg", "is_pos_pos", "is_bg_powered", "is_bg_all"]


def write_consolidated_tables(man, matrices, enr_by_cell):
    """Write two full-result tables spanning both cell types:
      tables/overlap_by_element.csv  — one row per element (positives+background),
        with Table S3 annotations, per-assay overlap counts, and a 0/1 column
        per ChIP-seq target.
      tables/enrichment_by_assay.csv — one row per (cell, target, effect_sign,
        background) with the full Fisher's-exact + BH-FDR statistics.
    """
    os.makedirs("tables", exist_ok=True)

    # ---- Table 1: per-element overlap, annotated ----
    per_elem = []
    for cell in CELLS:
        full, mat, tgt_assay = matrices[cell]
        tgt_cols = list(mat.columns)
        base = full.reset_index().rename(columns={"index": "elem_key"})
        ann = (man[man["cell"] == cell]
               .sort_values("n_sig", ascending=False)
               .drop_duplicates("elem_key")[["elem_key"] + ANN_COLS])
        df = base.merge(ann, on="elem_key", how="left")
        df.insert(0, "cell", cell)
        hist_t = [t for t in tgt_cols if tgt_assay.get(t, "") == "Histone ChIP-seq"]
        tf_t = [t for t in tgt_cols if tgt_assay.get(t, "") != "Histone ChIP-seq"]
        counts = pd.DataFrame({
            "n_targets_overlapped": mat[tgt_cols].sum(axis=1),
            "n_histone_overlapped": mat[hist_t].sum(axis=1) if hist_t else 0,
            "n_TF_overlapped": mat[tf_t].sum(axis=1) if tf_t else 0,
        }, index=mat.index).reindex(df["elem_key"].values).reset_index(drop=True)
        df = pd.concat([df, counts], axis=1)
        per_elem.append(df)
    all_cols = []
    for d in per_elem:
        all_cols += [c for c in d.columns if c not in all_cols]
    t1 = pd.concat([d.reindex(columns=all_cols) for d in per_elem], ignore_index=True)
    lead = (["cell"] + META_COLS + ANN_COLS +
            ["n_targets_overlapped", "n_histone_overlapped", "n_TF_overlapped"])
    lead = [c for c in lead if c in t1.columns]
    tgt_all = [c for c in t1.columns if c not in lead]
    t1[tgt_all] = t1[tgt_all].fillna(0).astype(int)
    t1 = t1[lead + tgt_all]
    t1.to_csv("tables/overlap_by_element.csv", index=False)
    print(f"tables/overlap_by_element.csv: {t1.shape[0]} elements x "
          f"{len(tgt_all)} target columns")

    # ---- Table 2: per-assay enrichment statistics ----
    front = ["cell", "target", "assay", "effect_sign", "background",
             "pos_overlap", "pos_total", "pos_frac", "bg_overlap", "bg_total",
             "bg_frac", "odds_ratio", "log2_OR", "p_value", "FDR"]
    parts = []
    for cell in CELLS:
        e = enr_by_cell[cell].copy()
        e.insert(0, "cell", cell)
        parts.append(e)
    t2 = pd.concat(parts, ignore_index=True)
    t2 = t2[[c for c in front if c in t2.columns] +
            [c for c in t2.columns if c not in front]]
    t2 = t2.sort_values(["cell", "effect_sign", "background", "FDR"]).reset_index(drop=True)
    t2.to_csv("tables/enrichment_by_assay.csv", index=False)
    print(f"tables/enrichment_by_assay.csv: {t2.shape[0]} rows "
          f"({(t2['FDR'] < 0.1).sum()} at FDR<0.1)")


if __name__ == "__main__":
    main()
