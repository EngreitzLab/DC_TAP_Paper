# Reviewer code-availability response (2026-10-12)

Goal: a versioned release (DC_TAP_Paper v1.2.0 + CRISPR_indirect_effects v1.0.0), Zenodo DOIs,
with one documented, tested command per revised figure (Figs 3c, 4, 5, 6 + S3/S4).

## A. CRISPR_indirect_effects repo
- [ ] Replace /oak paths in config/config.yml with repo-relative paths (inputs pulled from DC_TAP_Paper release / Zenodo)
- [ ] Deposit the outputs Fig 4/5 need: cis_positive_hit_rates.tsv, trans_positive_hit_rates.tsv, direct_effect_models/direct_effects_model.rds, *.direct_effect.tsv.gz (Zenodo; results/ is gitignored)
- [ ] Add env file pinning sceptre (see E); README run command
- [ ] Tag v1.0.0, mint Zenodo DOI

## B. Fig 4 (results/main_figure_4_direct_indirect_effects)
- [ ] indirect_effects_dir points to `../../CRISPR_indirect_effects_FDR_check/results` (unpublished dir name) -> point at downloaded release outputs
- [ ] Table S3 path `./Table_S3...tsv` not in that dir -> use results/formatted_dc_tap_results/ copy

## C. Figs 5-6 (results/main_figure_5_and_6)
- [ ] Fix missing comma: chromatin_categories_and_power_estimation.R:448
- [ ] housekeeping_genes.R:623, model_scores.R:177/181 hard-code /oak + /home paths -> relative paths + config block at top
- [ ] chromatin script reads chrom_annotate_results/, direct_effect_results/, ENCODE_rE2G_results/, ENCODE_rE2G/resources/... none in repo -> deposit (Zenodo) or ship the minimal derived tables
- [ ] crispr_benchmarking/dc_tap_seq_paper.yml: 7 /oak paths -> Zenodo-downloaded inputs
- [ ] Output dirs like 2026_0626_dc_tap_revisions -> results/figures/...

## D. One proven figure-generation command
- [ ] scripts/download_derived_inputs.sh (Zenodo -> resources/derived/)
- [ ] Snakemake target per figure (e.g. `snakemake --use-conda make_fig4 make_fig5 make_fig6`) or a single make_revised_figures rule
- [ ] Run it from a FRESH clone on Sherlock with no /oak access in env; diff PDFs/tables vs submitted panels; record runtime

## E. Environment / SCEPTRE pinning
- [ ] Record exact SCEPTRE source commit (Katsevich-Lab/sceptre SHA) behind jgalante::r-sceptre and r-sceptre-dev builds; ask Jud
- [ ] Pin in envs (version+build string) AND document `remotes::install_github("Katsevich-Lab/sceptre@<SHA>")` fallback
- [ ] Mirror the two jgalante conda packages (tarballs) on Zenodo in case the channel disappears
- [ ] Remove `#prefix: /oak/...` from sceptre_power_simulations.yml

## F. README refresh
- [ ] Remove stale `interpretation_analysis/` references (dir doesn't exist); point to results/main_figure_* dirs
- [ ] Checkout tag v1.0.0 -> v1.2.0; add figure -> script -> command table; cross-link CIE release
- [ ] Reconcile duplicated LICENSE + LICENSE.txt

## G. Release
- [ ] Tag v1.2.0, GitHub release, Zenodo DOI; update manuscript Code availability with both DOIs
- [ ] Response-to-reviewer paragraph listing each fix
