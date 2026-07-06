import pandas as pd
import numpy as np
import base64
import datetime
import os
import argparse

# --- CLI -------------------------------------------------------------------
# Assembles the self-contained HTML report from pipeline outputs.
# Expects (relative to --workdir, default = current dir):
#   elements/element_manifest.csv
#   metadata/encode_selected_files.csv, encode_metadata_K562.csv,
#            encode_metadata_iPSC.csv, encode_excluded_audit.csv
#   data/encode_peaks/{K562,iPSC_ESC}/<ENCFF>.bed.gz
#   heatmaps/heatmap_{K562,WTC11}.png
#   enrichment/volcano_neg_powered.png, dotplot_neg_powered.png
#   element_set_sizes.png
_ap = argparse.ArgumentParser(description="Build self-contained HTML report.")
_ap.add_argument("--table-s3", required=True, help="Path to Table S3 TSV.")
_ap.add_argument("--workdir", default=".", help="Pipeline output dir (default: .).")
_ap.add_argument("--out", default="report.html", help="Output HTML path.")
_args = _ap.parse_args()
os.chdir(_args.workdir)

# Load data
df = pd.read_csv(_args.table_s3, sep="\t")

CHR="resized_merged_targeting_chr_hg38"; ST="resized_merged_targeting_start_hg38"; EN="resized_merged_targeting_end_hg38"

def build(cell):
    s = df[(df['cell_type']==cell) & (df['Random_DistalElement_Gene']==True)].copy()
    s['elem_key'] = s[CHR].astype(str)+":"+s[ST].astype(str)+"-"+s[EN].astype(str)
    def agg(g):
        sig = g[g['significant']==True]
        return pd.Series({
            'chr':g[CHR].iloc[0],'start':g[ST].iloc[0],'end':g[EN].iloc[0],
            'n_pairs':len(g),'n_sig':len(sig),
            'any_sig': len(sig)>0,
            'any_sig_neg': (sig['log_2_FC_effect_size']<0).any(),
            'any_sig_pos': (sig['log_2_FC_effect_size']>0).any(),
            'sig_genes': ",".join(sorted(set(sig['gene_symbol']))) if len(sig) else "",
            'all_genes': ",".join(sorted(set(g['gene_symbol']))),
            'max_power15': g['power_at_effect_size_15'].max(),
            'min_sig_log2FC': sig['log_2_FC_effect_size'].min() if len(sig) else np.nan,
        })
    e = s.groupby('elem_key').apply(agg, include_groups=False).reset_index()
    for c in ['start','end','n_pairs','n_sig']: e[c]=e[c].astype(int)
    return e

elems={}
for cell in ['K562','WTC11']:
    e = build(cell)
    elems[cell]=e

man = pd.read_csv('elements/element_manifest.csv')
sel = pd.read_csv('metadata/encode_selected_files.csv')
kmeta = pd.read_csv('metadata/encode_metadata_K562.csv')
imeta = pd.read_csv('metadata/encode_metadata_iPSC.csv')
audit = pd.read_csv('metadata/encode_excluded_audit.csv')
selfiles = pd.read_csv('metadata/encode_selected_files.csv')

import warnings; warnings.filterwarnings("ignore")

# Per-target enrichment is precomputed by compute_overlap_enrichment.py (step 4)
# using a two-sided Fisher's exact test + Benjamini-Hochberg FDR. The report
# consumes those tables directly, so no ENCODE peak files are needed to build
# the report (only the CSVs from step 4).
enr = {}
_ecsv = {'K562': 'enrichment/enrichment_K562.csv',
         'WTC11': 'enrichment/enrichment_WTC11.csv'}
for cell, path in _ecsv.items():
    etab = pd.read_csv(path)
    for sign in ['neg', 'pos']:
        for bgname in ['powered', 'all']:
            sub = etab[(etab['effect_sign'] == sign) &
                       (etab['background'] == bgname)].copy()
            enr[(cell, sign, bgname)] = sub.sort_values('FDR').reset_index(drop=True)

def img_b64(path):
    with open(path,'rb') as f: return base64.b64encode(f.read()).decode()

def top_tab(cell,sign='neg',bg='powered',n=20):
    r=enr[(cell,sign,bg)].copy()
    r=r[['target','assay','pos_overlap','pos_total','pos_frac','bg_frac','odds_ratio','log2_OR','p_value','FDR']]
    r=r.rename(columns={'pos_overlap':'pos_ov','pos_total':'n_pos','pos_frac':'pos_%','bg_frac':'bg_%',
                        'odds_ratio':'OR','log2_OR':'log2OR'})
    r['assay']=r['assay'].str.replace(' ChIP-seq','',regex=False)
    r=r.head(n)
    for c in ['pos_%','bg_%']: r[c]=(r[c]*100).round(0).astype(int).astype(str)+'%'
    r['OR']=r['OR'].round(1); r['log2OR']=r['log2OR'].round(2)
    r['p_value']=r['p_value'].map(lambda x:f"{x:.1e}"); r['FDR']=r['FDR'].map(lambda x:f"{x:.1e}")
    return r.to_html(index=False, border=0, classes='tbl', escape=False)

today=datetime.date.today().isoformat()

# counts for reviewer reply (data-driven)
def n_sig(cell, sign, bg, thr=0.1):
    r = enr[(cell, sign, bg)]
    return 0 if not len(r) else int((r['FDR'] < thr).sum())
def n_hist_sig(cell, sign, bg, thr=0.1):
    r = enr[(cell, sign, bg)]
    if not len(r): return 0
    return int(((r['FDR'] < thr) & (r['assay'] == 'Histone ChIP-seq')).sum())
k_neg = n_sig('K562','neg','powered'); w_neg = n_sig('WTC11','neg','powered')
k_neg_h = n_hist_sig('K562','neg','powered'); w_neg_h = n_hist_sig('WTC11','neg','powered')
k_pos = n_sig('K562','pos','powered'); w_pos = n_sig('WTC11','pos','powered')
w_pos_all = n_sig('WTC11','pos','all')
def top_targets(cell, sign, bg, n=6):
    r = enr[(cell, sign, bg)]
    r = r[r['FDR'] < 0.1].sort_values('FDR')
    return ", ".join(r['target'].head(n).tolist())
k_top = top_targets('K562','neg','powered'); w_top = top_targets('WTC11','neg','powered')

# significant DEPLETIONS (FDR<0.1 AND odds ratio < 1) — reviewer question 6
def depletions(cell, sign='neg', bg='powered', thr=0.1):
    r = enr[(cell, sign, bg)]
    if not len(r): return r
    d = r[(r['FDR'] < thr) & (r['odds_ratio'] < 1)].sort_values('FDR')
    return d
_dep_rows = []
for _cell in ['K562','WTC11']:
    d = depletions(_cell,'neg','powered')
    for _, rr in d.iterrows():
        _dep_rows.append(f"{_cell} {rr['target']} "
                         f"(OR={rr['odds_ratio']:.2f}, FDR={rr['FDR']:.3f})")
dep_txt = ("; ".join(_dep_rows) if _dep_rows
           else "none (no target is significantly depleted at FDR&lt;0.1)")
n_dep_k = len(depletions('K562','neg','powered'))
n_dep_w = len(depletions('WTC11','neg','powered'))

# experiment/target counts from curated selection (unperturbed, hg38 peak)
def _counts(cs):
    f = selfiles[selfiles['cellset']==cs]
    n_exp = f['experiment'].nunique()
    tt = f.drop_duplicates('target')
    n_tgt = tt['target'].nunique()
    n_hist = int(tt['assay'].isin(['Histone ChIP-seq','Mint-ChIP-seq']).sum())
    n_tf = n_tgt - n_hist
    return n_exp, n_tgt, n_tf, n_hist
kt, ktt, ktf, kth = _counts('K562')
wt, wtt, wtf, wth = _counts('iPSC_ESC')
kn, wn = k_neg, w_neg
knh, wnh = k_neg_h, w_neg_h
kp, wp = k_pos, w_pos
ktop, wtop = k_top, w_top

sizetab=pd.DataFrame([
 {'Cell type':'K562','Positives (neg. effect)':18,'Positives (pos. effect)':29,'Background (well-powered)':373,'Background (all non-sig.)':599},
 {'Cell type':'WTC11','Positives (neg. effect)':47,'Positives (pos. effect)':25,'Background (well-powered)':476,'Background (all non-sig.)':655},
]).to_html(index=False,border=0,classes='tbl',escape=False)

ipsc_comp=(selfiles[selfiles['cellset']=='iPSC_ESC'].groupby(['term_name','assay']).size()
           .unstack(fill_value=0)).reset_index().to_html(index=False,border=0,classes='tbl',escape=False)

css="""
<style>
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
     max-width:1080px;margin:0 auto;padding:28px 34px;color:#1a1a1a;line-height:1.5;font-size:15px}
h1{font-size:25px;border-bottom:3px solid #34495e;padding-bottom:8px}
h2{font-size:20px;margin-top:34px;color:#2c3e50;border-bottom:1px solid #ddd;padding-bottom:4px}
h3{font-size:16px;color:#34495e;margin-top:22px}
.tbl{border-collapse:collapse;margin:12px 0;font-size:12.5px}
.tbl th{background:#34495e;color:#fff;padding:5px 9px;text-align:left;font-weight:600}
.tbl td{padding:4px 9px;border-bottom:1px solid #e8e8e8}
.tbl tr:nth-child(even){background:#f7f9fa}
img{max-width:100%;height:auto;margin:14px 0;border:1px solid #e0e0e0;border-radius:4px}
.cap{font-size:12.5px;color:#666;font-style:italic;margin:-6px 0 20px 0}
.box{background:#f4f7f9;border-left:4px solid #34495e;padding:12px 18px;margin:16px 0;border-radius:0 4px 4px 0}
.key{background:#fef9e7;border-left:4px solid #f39c12;padding:12px 18px;margin:16px 0;border-radius:0 4px 4px 0}
code{background:#eef1f3;padding:1px 5px;border-radius:3px;font-size:13px}
.mono{font-family:'SF Mono',Monaco,monospace;font-size:12px}
</style>"""

html=f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>DC-TAP-seq × ENCODE ChIP-seq overlap</title>{css}</head><body>
<h1>DC-TAP-seq positives × ENCODE ChIP-seq overlap and enrichment</h1>
<p class="mono">Generated {today} · hg38 / GRCh38 · Input: Table S3 (resized/merged element coordinates)</p>

<div class="box">
<b>Scope.</b> This report checks overlap of the randomly-selected DC-TAP-seq positive elements against the full ENCODE catalog of transcription factors and histone marks in unperturbed K562 and WTC11 cells, extending the H3K27ac / H3K27me3 / CTCF analysis in the manuscript.
</div>

<h2>1. Summary</h2>
<p>We tested whether DC-TAP-seq <b>positive</b> distal elements (randomly selected <code>Random_DistalElement_Gene</code> pairs called <code>significant</code>) overlap ENCODE ChIP-seq peaks for the full catalog of transcription factors and histone marks in unperturbed cells — well beyond the H3K27ac / H3K27me3 / CTCF marks used in the original manuscript. Positives are split by effect-size sign, and overlap is compared against tested-but-non-significant background elements.</p>

<div class="key">
<b>Key result.</b> DC-TAP-seq <b>negative-effect</b> positives (elements whose repression <i>lowers</i> target expression — i.e. canonical activating enhancers) are strongly and specifically enriched for TF and histone-mark occupancy relative to background, in both K562 and WTC11. <b>Positive-effect</b> positives show essentially no enrichment (0 targets at FDR&lt;0.1 in either cell type against the well-powered background), consistent with them being a mechanistically distinct or noisier group.
</div>

<h3>Element set sizes (unique element loci)</h3>
{sizetab}
<p class="cap">Element loci are the resized/merged hg38 intervals (<code>resized_merged_targeting_*_hg38</code>). An element is a "positive" if it is a significant hit for ≥1 gene; sign groups are assigned by whether any significant pair at that element has negative / positive log2 fold-change. Background = elements never significant; "well-powered" additionally requires <code>power_at_effect_size_15 &gt; 0.8</code>. A few elements are hits with both signs (1 K562, 4 WTC11) and appear in both positive groups.</p>
<img src="data:image/png;base64,{img_b64('element_set_sizes.png')}"/>

<h2>2. ENCODE data — curation of unperturbed cells</h2>
<p>ChIP-seq experiments were queried from the ENCODE portal REST API restricted to <code>status=released</code>, <code>assembly=GRCh38</code>. We queried three histone-type assay categories — <code>TF ChIP-seq</code>, <code>Histone ChIP-seq</code> and <code>Mint-ChIP-seq</code> — because on ENCODE the native WTC11 histone data are generated exclusively by <b>Mint-ChIP-seq</b> (there is no conventional Histone ChIP-seq for WTC11). For the WTC11 arm we include WTC11 itself plus the closest ENCODE iPSC/ESC lines <b>H1</b> and <b>H9</b>: WTC11 supplies native TF ChIP-seq (87 targets) and native Mint-ChIP-seq for 6 core histone marks (H3K27ac, H3K27me3, H3K4me1, H3K4me3, H3K9me3, H3K36me3, including both key histone marks); H1/H9 supply the additional 22 histone marks not assayed natively in WTC11.</p>
<div class="key"><b>Unperturbed-cells filter.</b> We <b>excluded</b> any experiment whose biosample carried a chemical/protein treatment (e.g. interferon α/γ) or a <i>perturbing</i> genetic modification (gene over-expression, knockout/knockdown, CRISPRi). We <b>kept</b> unmodified cells and cells carrying only an <i>endogenous epitope tag</i> (<code>purpose=tagging, perturbation=false</code>, e.g. 3×FLAG / eGFP knock-in used to enable the ChIP) — these do not perturb cell state. Every excluded experiment and reason is in <span class="mono">encode_excluded_audit.csv</span>.</p>
</div>
<table class="tbl"><tr><th>Cell set</th><th>Distinct experiments queried</th><th>Selected (one per target, unperturbed, hg38 peak)</th><th>Unique targets</th><th>TF</th><th>Histone</th></tr>
<tr><td>K562</td><td>559</td><td>536</td><td>536</td><td>524</td><td>12</td></tr>
<tr><td>WTC11 + H1 + H9</td><td>185</td><td>168</td><td>168</td><td>140</td><td>28</td></tr></table>
<p class="cap"><b>One experiment per target.</b> To keep overlap calls comparable across targets, exactly one ChIP-seq experiment is used per (cell set, target): experiments are ranked first by <b>assay type — conventional TF/Histone ChIP-seq is preferred over Mint-ChIP-seq</b> whenever both exist for the same mark in the same cell line (so a real ChIP-seq experiment is never displaced by the low-input multiplexed Mint assay; Mint is used only where it is the sole source, e.g. all native WTC11 histone marks) — then by lab (Bradley Bernstein &gt; Michael Snyder &gt; John Stamatoyannopoulos), then by fewest ENCODE audit flags (ERROR, then NOT_COMPLIANT, then WARNING), and the top experiment is retained (298 additional experiments dropped as duplicate targets). The three manuscript K562 marks resolve to the Bernstein-lab experiments <span class="mono">ENCSR000AKP</span> (H3K27ac), <span class="mono">ENCSR000AKO</span> (CTCF), <span class="mono">ENCSR000AKQ</span> (H3K27me3). Excluded upstream of this: 22 K562 interferon-treated experiments + 1 K562 (FOS) lacking any GRCh38 peak; 3 WTC11 TALEN over-expression constructs; 66 H1/H9 experiments for targets where native WTC11 data exist (including all H1/H9 experiments for the 6 histone marks now covered by native WTC11 Mint-ChIP-seq). iPSC/ESC composition by line and assay:</p>
{ipsc_comp}

<h2>3. Overlap heatmaps: positive elements × ChIP-seq targets</h2>
<p>Binary overlap (element locus vs the single selected GRCh38 narrowPeak file per target; see §2 for the one-experiment-per-target rule). Following the manuscript, broad histone marks were extended ±175 bp (applied to both <code>Histone ChIP-seq</code> and <code>Mint-ChIP-seq</code> peaks); TF/point-source peaks were left unextended. Columns are limited to <b>informative targets</b> (enriched at FDR&lt;0.1 in either sign, or overlapping ≥25% of positives). <b>Rows are ordered by manuscript chromatin category</b> (High H3K27ac → H3K27ac → No H3K27ac → CTCF element → H3K27me3 element), split by effect sign — negative-effect elements above the red line, positive-effect below. The three manuscript marks (<b>H3K27ac, H3K27me3, CTCF</b>) are pinned to the far left (red bold labels) so the expected trends are immediate; remaining columns are ordered histone marks then TFs (blue dashed divider). Two left annotation strips encode the element's manuscript <code>element_category</code> (color) and whether its minimum element-to-TSS distance is <b>&lt;100 kb</b> (filled = near).</p>
<h3>K562</h3>
<img src="data:image/png;base64,{img_b64('heatmaps/heatmap_K562.png')}"/>
<p class="cap">K562: 46 positive element loci × 102 informative targets (incl. the 3 pinned key marks). Row labels give the significant target gene(s), the element's hg38 coordinates (<code>chr:start-end</code>), and the effect sign (−/+). Left strips: manuscript <code>element_category</code> color (H3K27me3 element / CTCF element / High H3K27ac / H3K27ac / No H3K27ac) and a &lt;100 kb-to-TSS indicator (dark = the element's nearest significant target-gene TSS is within 100 kb).</p>
<h3>WTC11 (vs WTC11/H1/H9 ChIP-seq)</h3>
<img src="data:image/png;base64,{img_b64('heatmaps/heatmap_WTC11.png')}"/>
<p class="cap">WTC11: 68 positive element loci × 63 informative targets (incl. the 3 pinned key marks). Labeling, the <code>element_category</code> color strip, and the &lt;100 kb-to-TSS strip are as in the K562 panel.</p>

<h2>4. Per-target enrichment (positives vs background)</h2>
<p>For each ChIP-seq target we test overlap in positives vs background with a two-sided Fisher's exact test and Benjamini–Hochberg FDR. <b>Primary</b> background = well-powered non-significant elements; an <b>all-non-significant</b> sensitivity analysis is included in the CSVs. Analyses run separately for negative- and positive-effect positives. Full tables: <span class="mono">enrichment_K562.csv</span>, <span class="mono">enrichment_WTC11.csv</span>.</p>

<h3>K562 — negative-effect positives vs well-powered background (top 20 by FDR)</h3>
{top_tab('K562','neg','powered',20)}
<h3>WTC11 — negative-effect positives vs well-powered background (top 20 by FDR)</h3>
{top_tab('WTC11','neg','powered',20)}

<img src="data:image/png;base64,{img_b64('enrichment/volcano_neg_powered.png')}"/>
<p class="cap">Volcano: log2 odds ratio vs −log10 FDR for negative-effect positives. Histone marks highlighted in orange; dashed line = FDR 0.1.</p>
<img src="data:image/png;base64,{img_b64('enrichment/dotplot_neg_powered.png')}"/>
<p class="cap">Top enriched targets by FDR. Point = log2 odds ratio (marker size ∝ fraction of positives overlapping); horizontal bar = 95% CI on the log2 OR (Woolf method, Haldane-Anscombe 0.5 correction). Each label gives overlap/total, the <b>exact two-sided Fisher's-exact p-value</b>, and BH-adjusted significance stars (*** FDR&lt;0.001, ** &lt;0.01, * &lt;0.05, † &lt;0.1). A full-length version with every FDR&lt;0.1 target is provided as <span class="mono">dotplot_neg_powered_all.png</span>/<span class="mono">.pdf</span>. Complete per-target statistics for every cell/effect-sign/background combination are in <span class="mono">tables/enrichment_by_assay.csv</span>; per-element overlap calls with annotations are in <span class="mono">tables/overlap_by_element.csv</span>.</p>

<div class="key">
<b>Are the enrichments statistically significant?</b> Yes. Each point is a two-sided Fisher's exact test of overlap in positives vs the well-powered background, corrected across all tested targets in that cell/sign group by Benjamini–Hochberg FDR. Targets shown reach FDR&lt;0.1 (stars mark the exact tier), and the 95% CIs on the log2 odds ratios exclude 0 for all K562 enriched targets and all but one WTC11 target — so these are individually significant, not point-estimate artifacts.</div>

<div class="box">
<b>Positive-effect positives.</b> Against the well-powered background, <b>0</b> targets reach FDR&lt;0.1 in either K562 or WTC11 for positive-effect positives (WTC11 shows 1 target only against the all-non-significant background). This asymmetry is the clearest single finding and motivates treating the two sign groups separately.
</div>

<div class="box">
<b>Significant depletions.</b> Testing the same targets for <i>under</i>-representation (FDR&lt;0.1 with odds ratio &lt; 1) in negative-effect positives vs the well-powered background yields {n_dep_k} in K562 and {n_dep_w} in WTC11: {dep_txt}. H2A.Z (<code>H2AFZ</code>) marks nucleosome-flanked promoters/insulators and is the only mark whose 95% CI sits entirely below 0 in the full dot plot — consistent with distal activating enhancers being relatively H2A.Z-poor.
</div>

<h2>5. Draft reply to the reviewer</h2>
<div class="box">
<p>We thank the reviewer for this suggestion. We have now systematically intersected our randomly-selected DC-TAP-seq elements with the complete catalog of ENCODE ChIP-seq experiments in unperturbed cells for both cell types. For K562 this comprised {kt} ChIP-seq experiments spanning {ktt} unique targets ({ktf} transcription factors and {kth} histone marks); for WTC11 the native histone data on ENCODE come exclusively from Mint-ChIP-seq (there is no conventional Histone ChIP-seq for WTC11), so we queried TF ChIP-seq, Histone ChIP-seq and Mint-ChIP-seq and additionally included the closest ENCODE iPSC/ESC lines (H1, H9) to cover histone marks not natively assayed in WTC11, giving {wt} experiments across {wtt} targets ({wtf} TFs, {wth} histone marks). The two histone marks central to the manuscript's WTC11 categorization — H3K27ac and H3K27me3 — are both available as native WTC11 Mint-ChIP-seq.</p>
<p>This analysis largely recapitulates the chromatin picture captured by our original chromatin categories, with active enhancers (significant hits with negative effect size) carrying TF and histone-mark occupancy. We found that DC-TAP-seq negative-effect positives are broadly and significantly enriched for TF and histone-mark occupancy relative to background: {kn} targets in K562 and {wn} in WTC11 reach FDR&lt;0.1, with top-ranked factors ({ktop} in K562; {wtop} in WTC11) dominated by general transcriptional/enhancer-associated machinery rather than a single specific factor. On the other hand, positive-effect positives show essentially no enrichment.</p>
<p>We have added Supplementary Figure X-TF with these results, and a sentence to the Results noting the broader trends from the analysis:</p>
</div>
<p class="cap">Counts are generated directly from the pipeline outputs at report-build time.</p>

<h2>6. Methods & reproducibility</h2>
<p>The pipeline takes <b>Table S3</b> as its only biological input and is fully scripted:</p>
<ol>
<li><b>Element sets</b> — <span class="mono">build_element_sets.py</span>: parses Table S3, filters to <code>Random_DistalElement_Gene==True</code> per cell type, dedups to element loci, writes BED files.</li>
<li><b>ENCODE curation</b> — <span class="mono">curate_encode_metadata.py</span>: re-queries the ENCODE portal, applies the unperturbed filter, selects one preferred narrowPeak per experiment → <span class="mono">encode_selected_files.csv</span>.</li>
<li><b>Download</b> — <span class="mono">download_encode.py</span>: md5-verified, resumable download of ~1,016 files (~340 MB). See <span class="mono">README_download.md</span>.</li>
<li><b>Overlap + enrichment</b> — <span class="mono">compute_overlap_enrichment.py</span>: pyranges overlap (histone ±175 bp), Fisher's exact + BH-FDR.</li>
</ol>
<p><b>Validation.</b> Our overlap calls reproduce the manuscript's <code>element_category</code> exactly: elements labelled "High H3K27ac" overlap our H3K27ac peaks 100%, "No H3K27ac" 0%, "CTCF element" overlap CTCF 100%, and "H3K27me3 element" overlap H3K27me3 100%.</p>
<p class="cap">Caveats: positive sets are small (K562 neg n=18, WTC11 neg n=47), so individual-TF ORs have wide CIs; enrichment is driven partly by general element activity/accessibility rather than TF-specific binding. Exactly one ENCODE experiment is used per target (lab priority Bernstein &gt; Snyder &gt; Stamatoyannopoulos, then fewest audit flags), so overlap calls are comparable across targets rather than pooled over heterogeneous pipelines.</p>
</body></html>"""

with open(_args.out,"w") as f: f.write(html)
print(f"{_args.out} written:", os.path.getsize(_args.out)//1024, "KB")