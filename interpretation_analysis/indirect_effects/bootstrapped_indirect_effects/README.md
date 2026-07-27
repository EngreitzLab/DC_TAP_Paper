# Bootstrapped confidence intervals for the probability of direct effect

This small workflow adds **per-pair confidence intervals** to the "probability
of direct effect" reported for the DC-TAP-seq results. It responds to the
reviewer comment that the framework "provides no quantification of per-pair
uncertainty."

## What it does

The published estimate (from the upstream workflow
`distal_regulation_paper/CRISPR_indirect_effects_FDR_check`) is:

1. **indirect rate** per dataset = fraction of significant *trans* (other-chromosome)
   perturbation–gene pairs;
2. **cis rate** per distance-to-TSS bin (50 kb) = fraction of significant cis pairs;
3. **direct rate** = `max(cis_rate − indirect_rate, 0)`, averaged across the 7
   datasets, fit as a power law `direct_rate = a · dist^b` (`lm(log ~ log)`);
4. **per-pair directness** = `direct_rate(dist) / (direct_rate(dist) + indirect_rate)`,
   reported for positive and negative effects.

This workflow re-implements that chain and wraps it in a **non-parametric cluster
(block) bootstrap that resamples perturbations** (it does **not** re-run
sceptre/MAST — it reuses the existing differential-expression results). Each of
`B = 1000` replicates resamples, within every dataset, its perturbations with
replacement (a perturbation's cis and trans pairs move together), recomputes the
rates, re-fits the power law, and recomputes each pair's directness. The 95%
confidence interval is the percentile interval across replicates.

Resampling whole perturbations (rather than individual pairs) preserves the
within-perturbation correlation of trans effects and probes sensitivity to the
composition of the targeted panel.

## Inputs (read-only, from upstream workflows)

All 7 datasets that feed the fitted model are bootstrapped uniformly:

- 6 "held-out" datasets (sceptre pipeline), per-dataset under
  `.../CRISPR_indirect_effects_FDR_check/results/<dataset>/`:
  `results_run_discovery_analysis.rds` (trans DE) and
  `annotated_cis_results.tsv.gz` (cis DE; comma-separated despite the extension).
- 1 "training" dataset, Gasperini2019 (MAST pipeline) under
  `.../ENCODE_CRISPR_data/results/Gasperini2019/`:
  `diff_expr/output_MAST_perCRE.tsv.gz` (cis),
  `trans_effects/output_trans_effects_MAST_perCRE.tsv.gz` (trans), and the
  ENCODE-format `results/ENCODE/ENCODE_Gasperini2019_0.13gStd_MAST_perCRE_GRCh38.tsv.gz`
  (for the ValidConnection enhancer filter).
- Shared: the gene-universe TSS BED (`CollapsedGeneBounds.hg38.TSS500bp.bed`) and
  the DC-TAP results file to annotate
  (`Final_DC_TAP_Seq_Results_..._fdr20_direct_effects.tsv`).
- For validation only: the published `direct_effects_model.rds` and
  `trans_positive_hit_rates.tsv`.

Paths are set in the config block at the top of `bootstrap_directness_CIs.R`.

## How to run

R (>=4.2.0) with the sceptre and tidyverse packages needs to be available. 

```bash
WD=$(pwd)   # this directory

# 1) bootstrap (≈ a few minutes on 8 cores; SLURM_CPUS_PER_TASK enables parallelism)
srun -p normal -n1 -c8 --mem=24G -t 00:45:00 --chdir="$WD" \
  env -u R_HOME BOOT_B=1000 SLURM_CPUS_PER_TASK=8 "Rscript" bootstrap_directness_CIs.R

# 2) render the report
srun -p dev -n1 -c2 --mem=16G -t 00:20:00 --chdir="$WD" \
  env -u R_HOME "Rscript" -e 'rmarkdown::render("analyze_bootstrapped_directness.Rmd")'
```

Overridable environment variables: `BOOT_B` (replicates, default 1000),
`BOOT_SEED` (default 20250716), `SLURM_CPUS_PER_TASK` (cores).

The script asserts that its re-implementation reproduces the published power-law
coefficients and per-dataset indirect rates before running the bootstrap; it
stops with an error if they do not match.

## Outputs (`results/`)

- `bootstrap_powerlaw_coefficients.tsv` — per replicate (`replicate` 0 = point
  estimate): `type, intercept, slope, a, b`.
- `bootstrap_indirect_rates.tsv` — per replicate: `dataset, type, indirect_rate`.
- `directness_probability_CI_by_distance.tsv` — smooth curves for the ribbon plot:
  `dist_to_tss, cell_type, type, point, median, lower, upper`.
- `Final_DC_TAP_Seq_Results_..._fdr20_direct_effects_with_CIs.tsv` — the DC-TAP
  results file plus `direct_vs_indirect_{positive,negative}_{median,lower,upper}`
  for every pair.
- `bootstrap_summary.rds` — compact object consumed by the report.
- `analyze_bootstrapped_directness.html` — report creating supplementary figure.

## Files

- `bootstrap_directness_CIs.R` — computation (loaders → normalized table →
  cluster bootstrap → CIs).
- `analyze_bootstrapped_directness.Rmd` — report/plots.
