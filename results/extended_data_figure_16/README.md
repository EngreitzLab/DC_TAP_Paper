### Extended Data Fig. 16: Modeling probability of direct effects

To produce all required input for these analyses, first run the main DC-TAP-seq snakemake workflow
followed by the `direct_vs_indirect_effects/compute_direct_vs_indirect_effects` snakemake workflow
to compute main direct versus indirect effects outputs. Other external input data can be downloaded
from Zenodo using the `1_download_zenodo.sh` script.

To reproduce Extended Data Figure 16, run scripts in this order:
- `1_download_zenodo.sh`: Download all required external input files from Zenodo to `./data/`
- `2_indirect_effect_sizes_extdata_fig16a.Rmd`: Effect size distributions of indirect effects shown in Extended Data Fig. 16a
- `3_genome_wide_indirect_effects_extdata_fig16bc.Rmd`: Comparison of indirect effects of DC-TAP-seq elements and target genes versus genome-wide elements and target genes shown in Extended Data Fig. 16b,c. Input data produced by code in `../../direct_vs_indirect_effects/genome_wide_indirect_effects`
- `4_bootstrap_direct_effect_prob_CIs.R`: Bootstrap analysis to compute 95% confidence intervals for the modeled probability of direct effects
- `5_plot_direct_effect_prob_CIs_extdata_fig16d.Rmd`: Plot bootstrapped 95 CIs for modeled direct effects probability shown in Extended Data Fig. 16d
