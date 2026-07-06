# DC-TAP-seq positives × ENCODE ChIP-seq overlap and enrichment

Reproducible pipeline characterizing DC-TAP-seq **positive** distal elements by
their overlap with the full ENCODE catalog of transcription-factor and
histone-mark ChIP-seq peaks, in **K562** and **WTC11**, on **hg38 / GRCh38**.

This addresses the reviewer request to extend the manuscript's chromatin
characterization (which used H3K27ac, H3K27me3, CTCF) to the complete set of
ENCODE-profiled TFs and additional histone marks.

## What the pipeline does

1. **Element sets** — from Supplementary Table S3, define per cell type the
   positive elements (significant `Random_DistalElement_Gene` pairs, split by
   effect-size sign) and background elements (never significant), at the level
   of unique resized/merged hg38 element loci.
2. **ENCODE curation** — query the ENCODE portal for released GRCh38 ChIP-seq,
   keep only **unperturbed** cells (unmodified or endogenously epitope-tagged),
   and select one preferred narrowPeak per experiment.
3. **Download** — fetch the curated peak files (~1,016 files, ~350 MB), md5-verified.
4. **Overlap + enrichment** — build element × target binary overlap matrices
   (histone marks — both `Histone ChIP-seq` and `Mint-ChIP-seq` — extended
   ±175 bp, matching the manuscript), then test
   per-target enrichment in positives vs background (Fisher's exact + BH-FDR),
   separately for negative- and positive-effect positives and for two
   backgrounds.
5. **Figures + report** — overlap heatmaps, enrichment volcano/dot plots,
   and a self-contained HTML report.

## Key result

DC-TAP-seq **negative-effect** positives (repression lowers target expression —
canonical activating enhancers) are strongly and specifically enriched for TF
and histone-mark occupancy vs background, in both cell types. **Positive-effect**
positives show essentially no enrichment (0 targets at FDR<0.1 vs the
well-powered background in either cell type). Overlap calls reproduce the
manuscript's `element_category` labels exactly (High H3K27ac → H3K27ac 100%,
No H3K27ac → 0%, CTCF element → CTCF 100%, H3K27me3 element → H3K27me3 100%).

Enrichment significance is assessed by a two-sided **Fisher's exact test** per
target (positives vs background), corrected across all tested targets within
each cell/effect-sign/background group by **Benjamini–Hochberg FDR**. Testing
the same targets for depletion (FDR<0.1, odds ratio < 1) yields a single
significant result: **WTC11 H2A.Z (`H2AFZ`), OR≈0.40, FDR≈0.044** vs the
well-powered background; K562 has no significant depletions.

## Inputs

- **Table S3** (the only biological input):
  `Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_*.tsv`
  Pass its path via `--table-s3`. Required columns: `cell_type`,
  `Random_DistalElement_Gene`, `significant`, `log_2_FC_effect_size`,
  `gene_symbol`, `power_at_effect_size_15`, and the resized/merged hg38
  coordinates `resized_merged_targeting_{chr,start,end}_hg38`.
- **ENCODE peak data** — downloaded by the pipeline (public, no account needed).

## Environment

Conda (recommended — includes bedtools):

    conda env create -f environment.yml
    conda activate encode

Or pip for the Python parts (install bedtools separately):

    pip install -r requirements.txt

Dependencies: python 3.11, pandas, numpy, scipy, matplotlib, seaborn,
statsmodels, requests, pyranges, bedtools ≥2.31.

## Run order

All scripts use paths relative to the project root; run from there.

    # 1. Element sets from Table S3  ->  elements/*.bed, element_manifest.csv
    python scripts/build_element_sets.py --table-s3 /path/to/Table_S3.tsv --outdir elements

    # 2. (optional) Re-query ENCODE and rebuild the file manifest
    #    A pre-built metadata/ is shipped; run this only to refresh it.
    python scripts/curate_encode_metadata.py --outdir metadata

    # 3. Download curated ENCODE peaks (~350 MB)  ->  data/encode_peaks/{K562,iPSC_ESC}/
    python scripts/download_encode.py --manifest metadata/encode_selected_files.csv --outdir data/encode_peaks

    # 4. Overlap matrices + enrichment tables
    python scripts/compute_overlap_enrichment.py

    # 5. Figures (heatmaps, volcano, dot plots with 95% CI)
    #    Every figure is written as BOTH .png and a vector .pdf with editable
    #    (TrueType-embedded) text for Illustrator.
    python scripts/make_figures.py

    # 6. Self-contained HTML report (embeds all figures as base64).
    #    Reads the precomputed enrichment tables from step 4 and the figures
    #    from step 5 — it does NOT re-read the ENCODE peak files, so the report
    #    can be rebuilt without the ~350 MB download present.
    python scripts/build_report.py --table-s3 /path/to/Table_S3.tsv

See `scripts/README_download.md` for details on the ENCODE curation and download.

## Outputs

    elements/
      {K562,WTC11}_{pos_neg,pos_pos,bg_powered,bg_all}.bed   element sets (BED4)
      element_manifest.csv                                    element -> set + gene/power info
    metadata/
      encode_metadata_K562.csv, encode_metadata_iPSC.csv      all experiments + keep/exclude
      encode_excluded_audit.csv                               every exclusion + reason
      encode_selected_files.csv                               download manifest (URLs, md5)
    overlap/
      overlap_matrix_{K562,WTC11}.csv                         elements x targets binary matrix
      overlap_long_{K562,WTC11}.csv                           element-experiment-target long form
      target_assay_{K562,WTC11}.csv                           target -> TF/Histone map
    enrichment/
      enrichment_{K562,WTC11}.csv                             all sign x background enrichment
                                                             (odds_ratio, p_value, FDR, log2_OR, counts)
      volcano_neg_powered.{png,pdf}                          log2 OR vs -log10 FDR
      dotplot_neg_powered.{png,pdf}                          top-18 targets, log2 OR + 95% CI + FDR stars
      dotplot_neg_powered_all.{png,pdf}                      ALL FDR<0.1 targets (reference)
    tables/
      overlap_by_element.csv                                 ONE ROW PER ELEMENT (both cells): Table S3
                                                             annotations + set membership + per-assay
                                                             overlap counts + a 0/1 column per target
      enrichment_by_assay.csv                                ONE ROW PER (cell,target,effect_sign,
                                                             background): full Fisher's-exact + BH-FDR stats
    heatmaps/
      heatmap_{K562,WTC11}.{png,pdf}                         rows ordered by Table S3 element_category
                                                             (active->insulator->repressive), split by
                                                             effect sign; H3K27ac/H3K27me3/CTCF columns
                                                             pinned far left; two left strips =
                                                             element_category color + <100 kb-to-TSS flag
    report.html                                               self-contained report (incl. draft reviewer reply)

All figures are additionally written as vector **.pdf** with editable text
alongside the .png shown above.

## Element set sizes (unique element loci)

| Cell  | pos (neg effect) | pos (pos effect) | bg (well-powered) | bg (all non-sig) |
|-------|------------------|------------------|-------------------|------------------|
| K562  | 18               | 29               | 373               | 599              |
| WTC11 | 47               | 25               | 476               | 655              |

A few elements are hits with both signs (1 in K562, 4 in WTC11) and appear in
both positive groups. Background always excludes any element that is ever a
positive.

## Caveats

Positive sets are small (K562 neg n=18, WTC11 neg n=47), so per-TF odds ratios
have wide confidence intervals; enrichment partly reflects general element
activity/accessibility rather than TF-specific binding. Exactly one ENCODE
experiment is used per target (see "One experiment per target" below), so
overlap calls are comparable across targets rather than pooled over
heterogeneous pipelines.

## One experiment per target

To keep overlap calls comparable across targets, exactly one ChIP-seq
experiment is retained per (cell set, target). Experiments are ranked first by
**assay type — conventional TF/Histone ChIP-seq is preferred over
Mint-ChIP-seq** whenever both exist for the same mark in the same cell line, so
a genuine ChIP-seq experiment is never displaced by the low-input multiplexed
Mint assay on the basis of lab or accession; Mint-ChIP-seq is retained only
where it is the sole source of a mark (e.g. all native WTC11 histone marks,
which have no conventional Histone ChIP-seq on ENCODE). Within the same assay
tier, experiments are ranked by lab in priority order **Bradley Bernstein >
Michael Snyder > John Stamatoyannopoulos**; within the top-priority lab, by
fewest ENCODE audit flags (ERROR, then NOT_COMPLIANT, then WARNING); ties
broken by accession for determinism. This composes after the native-WTC11 preference (below), so
single-experiment selection runs within the preferred biosample pool. The
three manuscript K562 marks resolve to the Bernstein-lab experiments
`ENCSR000AKP` (H3K27ac), `ENCSR000AKO` (CTCF), `ENCSR000AKQ` (H3K27me3).
Implemented in `curate_encode_metadata.py::select_one_per_target()`; 298
duplicate-target experiments are dropped (reason recorded in
`encode_excluded_audit.csv`).

**Histone assays — Mint-ChIP-seq included.** On ENCODE the native WTC11
histone data are generated *exclusively* by **Mint-ChIP-seq**; WTC11 has no
conventional `Histone ChIP-seq`. The curation therefore queries three assay
categories — `TF ChIP-seq`, `Histone ChIP-seq` and `Mint-ChIP-seq` — and treats
both `Histone ChIP-seq` and `Mint-ChIP-seq` as histone marks (classified as
`Histone ChIP-seq` downstream, extended ±175 bp). This is what allows the two
manuscript WTC11 histone categories, **H3K27ac** and **H3K27me3**, to be scored
from native WTC11 data rather than from the H1/H9 fallback. (For example, the
IGFBP2 positive-effect element chr2:216,370,766-216,371,731, categorized
"H3K27me3 element" in the manuscript, overlaps the native WTC11 H3K27me3
Mint-ChIP-seq peak ENCSR418YEV.)

**WTC11 arm — biosample preference.** The WTC11 arm draws on ENCODE data from
WTC11 plus the embryonic stem-cell lines H1 and H9 (native WTC11 coverage is
sparse for many targets). Data are combined with a *native-WTC11-preferred*
rule: for any target that has at least one WTC11 experiment, only the WTC11
experiment(s) are used and the H1/H9 experiments for that same target are
dropped; H1/H9 are retained only as a fallback for targets with no WTC11
experiment at all. This is applied in
`curate_encode_metadata.py::prefer_native_wtc11()`; the dropped H1/H9 files
appear in `encode_excluded_audit.csv` with reason "WTC11 data available for this
target". Of the 168 WTC11-arm targets, 93 have native WTC11 data (87 TF ChIP-seq
+ 6 Mint-ChIP-seq histone marks: H3K27ac, H3K27me3, H3K4me1, H3K4me3, H3K9me3,
H3K36me3) and 75 are H1/H9-only (53 TF + 22 histone marks not natively assayed
in WTC11). 16 targets (ATF2, ATF3, CREB1, CTCF, JUN, MAX, SP1, TEAD4, USF1,
USF2, and the 6 histone marks above) previously mixed WTC11 with H1/H9 and now
use WTC11 only (66 H1/H9 files dropped). This keeps the pooled overlaps
consistent with the manuscript's WTC11-specific `element_category`, so a site is
not called bound on the basis of ESC-line data when a native WTC11 experiment
exists.
