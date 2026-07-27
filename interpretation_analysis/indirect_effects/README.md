### Direct versus indirect effects analyses

This directory contains code used for direct vs. indirect effects analyses presented in Figure 4 and
Supplementary Figures 3, 15, and 16. Direct, indirect effect rates and probabilities of direct
effects are computed by the separate workflow https://github.com/EngreitzLab/CRISPR_indirect_effects.
Running this workflow is required to fully reproduce analysis in this directory.

Code for the following analyses are in this directory:
- `dc_tapseq_indirect_effects_main_analyses.Rmd`: Main analyses shown in Figures 4 and S3
- `dc_tapseq_indirect_effects_supplementary_analyses.Rmd`: Supplementary analyses and figures for modeling direct and indirect effect rates shown in Supplementary Figure 15
- `dc_tapseq_indirect_effect_sizes.Rmd`: Effect size distributions of indirect effects shown in Supplementary Figure 16a
- `genome_wide_indirect_effects`: Comparison of indirect effects of DC-TAP-seq elements and target genes versus genome-wide elements and target genes shown in Supplementary Figures 16b,c. Scripts need to be run sequentially
- `bootstrapped_indirect_effects`: Bootstrap analysis to compute 95% confidence intervals for the modeled probability of direct effects

