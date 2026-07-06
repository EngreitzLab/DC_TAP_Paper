#!/usr/bin/env python3
"""
Re-query the ENCODE portal and rebuild the curated, unperturbed ChIP-seq file
manifest used by the DC-TAP-seq overlap analysis (hg38 / GRCh38 only).

Outputs (into --outdir, default metadata/):
    encode_metadata_K562.csv     all K562 experiments with keep/exclude decision
    encode_metadata_iPSC.csv     all WTC11/H1/H9 experiments with keep/exclude decision
    encode_excluded_audit.csv    every excluded experiment + reason
    encode_selected_files.csv    kept experiments + chosen narrowPeak file (download manifest)

Usage:
    python scripts/curate_encode_metadata.py --outdir metadata

Dependencies: python>=3.8, requests, pandas
"""
import argparse, os, requests, pandas as pd
from collections import Counter

PORTAL = "https://www.encodeproject.org/search/"
S = requests.Session(); S.headers.update({"Accept": "application/json"})

FILE_FIELDS = ["accession","assay_title","target.label","target.investigated_as",
    "biosample_ontology.term_name","biosample_summary","status",
    "lab.title","audit",
    "replicates.library.biosample.genetic_modifications.purpose",
    "replicates.library.biosample.genetic_modifications.method",
    "replicates.library.biosample.genetic_modifications.perturbation",
    "replicates.library.biosample.genetic_modifications.modified_site_by_target_id.label",
    "replicates.library.biosample.genetic_modifications.introduced_tags",
    "replicates.library.biosample.treatments.treatment_term_name",
    "replicates.library.biosample.treatments.treatment_type",
    "files.accession","files.file_type","files.output_type","files.assembly",
    "files.status","files.preferred_default","files.href","files.file_size","files.md5sum"]

PRIORITY = ["optimal IDR thresholded peaks","conservative IDR thresholded peaks",
    "IDR thresholded peaks","pseudoreplicated IDR thresholded peaks",
    "replicated peaks","pseudoreplicated peaks","peaks"]

def fetch(term_names, assays):
    out = []
    for term in term_names:
        for assay in assays:
            params = [("type","Experiment"),("assay_title",assay),
                ("biosample_ontology.term_name",term),("status","released"),
                ("assembly","GRCh38"),("format","json"),("limit","all")]
            for f in FILE_FIELDS: params.append(("field", f))
            r = S.get(PORTAL, params=params, timeout=180)
            if r.status_code == 200: out.extend(r.json().get("@graph", []))
    return out

def biosamples(exp):
    for rep in exp.get("replicates", []):
        bs = rep.get("library", {}).get("biosample", {})
        if bs: yield bs

def classify(exp):
    reasons, tags = [], set()
    for bs in biosamples(exp):
        for gm in bs.get("genetic_modifications", []):
            purpose, pert = gm.get("purpose"), gm.get("perturbation")
            if purpose == "tagging" and pert is False:
                for t in (gm.get("introduced_tags") or []): tags.add(t.get("name",""))
            elif pert is True or purpose not in ("tagging", None):
                tgt = (gm.get("modified_site_by_target_id") or {}).get("label","?")
                reasons.append(f"perturbing GM: purpose={purpose},method={gm.get('method')},target={tgt}")
        for t in bs.get("treatments", []):
            reasons.append(f"treatment: {t.get('treatment_term_name')} ({t.get('treatment_type')})")
    return len(reasons) == 0, "; ".join(sorted(set(reasons))), ",".join(sorted(t for t in tags if t))

def pick_file(exp):
    cands = [f for f in exp.get("files", [])
             if f.get("file_type") == "bed narrowPeak" and f.get("assembly") == "GRCh38"
             and f.get("status") == "released"]
    if not cands: return None
    pref = [f for f in cands if f.get("preferred_default")]
    pool = pref if pref else cands
    pool.sort(key=lambda f: PRIORITY.index(f["output_type"]) if f.get("output_type") in PRIORITY else len(PRIORITY))
    return pool[0]

def curate(explist, cellset):
    rows = []
    for exp in explist:
        keep, reason, tags = classify(exp)
        f = pick_file(exp)
        if f is None and keep:
            keep, reason = False, "no GRCh38 narrowPeak (hg19 only)"
        audit = exp.get("audit", {}) or {}
        rows.append({
            "cellset": cellset,
            "term_name": (exp.get("biosample_ontology") or {}).get("term_name",""),
            "experiment": exp["accession"], "assay": exp.get("assay_title",""),
            "target": (exp.get("target") or {}).get("label",""),
            "investigated_as": ";".join((exp.get("target") or {}).get("investigated_as", [])),
            "biosample_summary": exp.get("biosample_summary",""),
            "lab": (exp.get("lab") or {}).get("title",""),
            "n_audit_error": len(audit.get("ERROR", [])),
            "n_audit_notcompliant": len(audit.get("NOT_COMPLIANT", [])),
            "n_audit_warning": len(audit.get("WARNING", [])),
            "endogenous_tag": tags, "keep": keep, "exclude_reason": reason,
            "file": f["accession"] if f else None,
            "output_type": f.get("output_type") if f else None,
            "preferred_default": bool(f.get("preferred_default")) if f else None,
            "file_size": f.get("file_size") if f else None,
            "md5sum": f.get("md5sum") if f else None,
            "href": ("https://www.encodeproject.org" + f.get("href","")) if f else None,
        })
    return pd.DataFrame(rows)

def prefer_native_wtc11(df):
    """For the WTC11 (iPSC_ESC) arm: if a target has any kept native-WTC11
    experiment, drop the H1/H9 experiments for that SAME target (they are only
    a fallback for targets with no WTC11 data). Targets with no WTC11 experiment
    keep their H1/H9 data. Mutates 'keep'/'exclude_reason' in place and returns
    the dataframe. Only touches rows where cellset == 'iPSC_ESC'.
    """
    arm = df["cellset"] == "iPSC_ESC"
    # targets that have >=1 kept native-WTC11 experiment
    wtc11_targets = set(
        df.loc[arm & (df["term_name"] == "WTC11") & df["keep"], "target"]
    )
    drop = arm & df["keep"] & (df["term_name"] != "WTC11") & df["target"].isin(wtc11_targets)
    df.loc[drop, "exclude_reason"] = "WTC11 data available for this target (H1/H9 fallback not used)"
    df.loc[drop, "keep"] = False
    return df

LAB_RANK = [("bernstein", 0), ("snyder", 1), ("stamatoyannopoulos", 2)]

def _lab_priority(lab):
    l = str(lab).lower()
    for name, rank in LAB_RANK:
        if name in l:
            return rank
    return 3  # any other lab -> fall back to data-quality flags

def _assay_priority(assay):
    """Within a (cellset, target), prefer a conventional histone-mark assay over
    the low-input multiplexed Mint-ChIP-seq assay when BOTH exist for the same
    mark in the same cell line. Mint-ChIP-seq is used only where it is the sole
    (or preferred-biosample) source — e.g. WTC11, which has no conventional
    Histone ChIP-seq on ENCODE. TF ChIP-seq and Histone ChIP-seq rank equal (0);
    Mint-ChIP-seq ranks below (1). This tier is applied BEFORE lab/audit
    tie-breakers so a real conventional experiment is never displaced by Mint on
    the basis of lab or accession alone."""
    return 1 if str(assay) == "Mint-ChIP-seq" else 0

def select_one_per_target(df):
    """Keep exactly ONE experiment per (cellset, target) so overlap calls are
    comparable across targets. Selection rule, applied only to currently-kept
    rows: (1) prefer a conventional assay (TF or Histone ChIP-seq) over
    Mint-ChIP-seq when both exist for the same mark in the same cell line, so a
    real ChIP-seq experiment is never displaced by the low-input multiplexed
    Mint assay on the basis of lab/accession — Mint is retained only where it is
    the sole source (e.g. all native WTC11 histone marks, which have no
    conventional Histone ChIP-seq on ENCODE); (2) prefer lab in priority order
    Bernstein > Snyder > Stamatoyannopoulos; (3) among the top-priority lab,
    prefer the experiment with the fewest ENCODE audit flags (ERROR, then
    NOT_COMPLIANT, then WARNING); (4) break any remaining tie by accession for
    determinism. All other kept experiments for that target are set keep=False.
    Composes after prefer_native_wtc11 (which has already removed H1/H9 fallbacks
    where native WTC11 exists), so single-experiment selection runs within the
    preferred biosample pool. Mutates 'keep'/'exclude_reason' in place; returns
    df.
    """
    d = df.copy()
    d["_lab_priority"] = d["lab"].map(_lab_priority)
    d["_assay_priority"] = d["assay"].map(_assay_priority)
    kept = d[d["keep"] & d["file"].notna()].copy()
    order = kept.sort_values(
        ["cellset", "target", "_assay_priority", "_lab_priority", "n_audit_error",
         "n_audit_notcompliant", "n_audit_warning", "experiment"],
        kind="mergesort")
    chosen = order.groupby(["cellset", "target"]).first().reset_index()
    chosen_exp = set(zip(chosen["cellset"], chosen["target"], chosen["experiment"]))
    for idx, row in df.iterrows():
        if not (row["keep"] and pd.notna(row["file"])):
            continue
        if (row["cellset"], row["target"], row["experiment"]) not in chosen_exp:
            df.at[idx, "keep"] = False
            df.at[idx, "exclude_reason"] = (
                "not selected: one experiment per target "
                "(conventional ChIP-seq over Mint-ChIP, then lab priority "
                "Bernstein>Snyder>Stamatoyannopoulos, then audit flags)")
    return df

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--outdir", default="metadata")
    args = ap.parse_args(); os.makedirs(args.outdir, exist_ok=True)

    # Mint-ChIP-seq is ENCODE's low-input multiplexed histone-ChIP assay. All
    # native WTC11 histone marks (H3K27ac/me3, H3K36me3, H3K4me1/3, H3K9me3) are
    # released under this assay title, so it MUST be queried alongside the
    # conventional "Histone ChIP-seq" title or WTC11 histone coverage is missed.
    HIST_ASSAYS = ["Histone ChIP-seq", "Mint-ChIP-seq"]
    k = curate(fetch(["K562"], ["TF ChIP-seq"] + HIST_ASSAYS), "K562")
    i = curate(fetch(["WTC11","H1","H9"], ["TF ChIP-seq"] + HIST_ASSAYS), "iPSC_ESC")
    i = prefer_native_wtc11(i)   # native WTC11 preferred; H1/H9 only fill gaps
    allm = pd.concat([k, i], ignore_index=True)
    allm = select_one_per_target(allm)  # exactly one experiment per (cellset,target)

    k.to_csv(os.path.join(args.outdir, "encode_metadata_K562.csv"), index=False)
    i.to_csv(os.path.join(args.outdir, "encode_metadata_iPSC.csv"), index=False)
    allm[~allm["keep"]].to_csv(os.path.join(args.outdir, "encode_excluded_audit.csv"), index=False)
    sel = allm[allm["keep"] & allm["file"].notna()].copy()
    sel.to_csv(os.path.join(args.outdir, "encode_selected_files.csv"), index=False)
    print(f"kept={len(sel)} excluded={(~allm['keep']).sum()} "
          f"K562_targets={sel[sel.cellset=='K562'].target.nunique()} "
          f"iPSC_targets={sel[sel.cellset=='iPSC_ESC'].target.nunique()}")

if __name__ == "__main__":
    main()
