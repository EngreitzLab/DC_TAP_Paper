This directory contains code, intermediate files, and results used for the chromatin category and housekeeping gene analysis of DC-TAP-seq results and comparison to previous screens.

- `chromatin_categories_and_power_estimation.R`: code to produce the chromatin category plot in Fig. 5 (and associated supplementary figure), and the power estimation plot in Fig. 3c
- `categorized_pairs`: outputs from `chromatin_categories_and_power_estimation.R` used for Fig. 5; produced using results from ENCODE-rE2G, chrom-annotate, and DC-TAP-seq analysis pipelines
  - `dc_tap_gasperini.random_distal_element_pairs.categorized.tsv.gz`: DC-TAP random screen DE-G pairs and Gasperini DE-G pairs with chromatin category, promoter class, direct effect rate, and power annotations
  - `other_crispr_screens.categorized.tsv.gz`: same as above but for CRISPR perturbation data from 8 other studies
  - `all_gw_pairs.categorized.summary.tsv`
  - `thresholds.tsv`: quantiles for various chromatin features for distal elements in 5 cell types used for categorization

- `housekeeping_genes.R`: code to reproduce the housekeeping gene analysis in Fig. 5 and associated supplementary figures from Table S3

- `model_scores.R`: code to reproduce the Fig. 6 and associated supplementary figures from CRISPR_comparison benchmarking pipeline output 