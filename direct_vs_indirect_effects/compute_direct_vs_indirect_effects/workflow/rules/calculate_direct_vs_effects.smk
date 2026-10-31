## Calculate probability of direct vs. indirect CRISPR hits

# fit powerlaw models to predict direct hit rates as function of distance to TSS
rule model_direct_effects_rate:
  input:
    cis_rates = "results/cis_positive_hit_rates.tsv",
    trans_rates = "results/trans_positive_hit_rates.tsv"
  output:
    models = "results/direct_effect_models/direct_effects_model.rds",
    direct_rates_datasets = "results/direct_effect_models/direct_rates_per_datasets.tsv",
    direct_rates_average  = "results/direct_effect_models/direct_rates_average_across_datasets.tsv"
  params:
    seed = 20250508,
    max_distance = 1e6
  conda: "../envs/sceptre.yml"
  script:
    "../scripts/model_direct_effects_rate.R"
    
# impute trans positive hit rates for CRISPR datasets without empirically computed hit rates
rule impute_trans_effects_rate:
  input:
    trans_rates = "results/trans_positive_hit_rates.tsv",
    dataset_ids = "config/crispr_dataset_ids.tsv"
  output: "results/trans_positive_hit_rates_imputed.tsv"
  conda: "../envs/sceptre.yml"
  script:
    "../scripts/impute_trans_effects_rate.R"

# plot rates and modelled direct rates across all datasets
rule analyze_modeled_direct_rates:
  input:
    cis_rates = "results/cis_positive_hit_rates.tsv",
    trans_rates = "results/trans_positive_hit_rates.tsv",
    models = "results/direct_effect_models/direct_effects_model.rds",
    direct_rates_datasets = "results/direct_effect_models/direct_rates_per_datasets.tsv",
    direct_rates_average  = "results/direct_effect_models/direct_rates_average_across_datasets.tsv",
    dataset_ids = "config/crispr_dataset_ids.tsv"
  output: "results/model_indirect_effects.html"
  params:
    max_distance = 1e6
  conda: "../envs/sceptre.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/model_indirect_effects.Rmd"
    
# predict direct vs indirect effects probability rate for DC-TAP-seq data
rule predict_direct_vs_indirect_effect_dctap:
  input:
    crispr = "resources/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_fdr20.tsv.gz",
    models = "results/direct_effect_models/direct_effects_model.rds",
    indirect_rates = "results/trans_positive_hit_rates_imputed.tsv",
    dataset_ids = "config/crispr_dataset_ids.tsv",
    gene_tss = "resources/CollapsedGeneBounds.hg38.TSS500bp.bed"
  output: "results/annotated_crispr_data/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_fdr20_direct_effects.tsv"
  params:
    dist_col = "distance_to_gencode_gene_TSS"
  conda: "../envs/sceptre.yml"
  script:
    "../scripts/predict_direct_vs_indirect_effect_dctap.R"        
