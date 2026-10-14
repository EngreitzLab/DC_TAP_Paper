## Rules to create main Figures 4-6 (and Fig. 3c) from tables in the repository.
## These rules have no conda directive: they run in the top-level pixi environment (pixi.toml).

# Figure 4 (and Fig. S3): direct vs indirect effects, from CRISPR_indirect_effects workflow outputs
rule main_figure_4_direct_indirect_effects:
  input:
    rmd = "results/main_figure_4_direct_indirect_effects/dc_tapseq_indirect_effects.Rmd",
    cis_rates = "resources/main_figure_4/indirect_effects/cis_positive_hit_rates.tsv",
    trans_rates = "resources/main_figure_4/indirect_effects/trans_positive_hit_rates.tsv",
    direct_effect_model = "resources/main_figure_4/indirect_effects/direct_effect_models/direct_effects_model.rds",
    table_s3 = "results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv"
  output:
    html = "results/main_figure_4_direct_indirect_effects/dc_tapseq_indirect_effects.html",
    rates = "results/main_figure_4_direct_indirect_effects/plots/cis_and_trans_effect_rates.pdf",
    p_direct = "results/main_figure_4_direct_indirect_effects/plots/probability_direct_vs_indirect.pdf",
    effect_direction = "results/main_figure_4_direct_indirect_effects/plots/negative_vs_positive_effects_combined.pdf",
    effect_size = "results/main_figure_4_direct_indirect_effects/plots/probability_direct_vs_effect_size.pdf"
  log: "results/main_figure_4_direct_indirect_effects/logs/main_figure_4.log"
  params:
    out_dir = "results/main_figure_4_direct_indirect_effects"
  shell:
    """
    Rscript -e 'wd <- getwd(); rmarkdown::render("{input.rmd}", output_dir = file.path(wd, "{params.out_dir}"),
      params = list(indirect_effects_dir = file.path(wd, dirname("{input.cis_rates}")),
                    table_s3_file = file.path(wd, "{input.table_s3}"),
                    plots_dir = file.path(wd, "{params.out_dir}/plots")))' > {log} 2>&1
    """

# Figure 5b, 5c,d (with Gasperini et al.) and Figure 3c, from categorized pairs in the repository
rule main_figure_5_chromatin_categories:
  input:
    script = "results/main_figure_5_and_6/chromatin_categories_and_power_estimation.R",
    functions = "results/main_figure_5_and_6/housekeeping_gene_plot_functions.R",
    dc_tap_gasperini = "results/main_figure_5_and_6/categorized_pairs/dc_tap_gasperini.random_distal_element_pairs.categorized.tsv.gz",
    genome_wide_pairs = "results/main_figure_5_and_6/categorized_pairs/all_gw_pairs.categorized.summary.tsv"
  output:
    fig5b = "results/main_figure_5_and_6/categorized_pairs/proportion_by_cell_type.element_category.de_weighted.combine_dc_tap.pdf",
    fig5b_by_cell_type = "results/main_figure_5_and_6/categorized_pairs/proportion_by_cell_type.element_category.de_weighted.for_dc_tap.pdf",
    fig5c = "results/main_figure_5_and_6/housekeeping_genes/v3_combinedDC_TAP_all_metrics_comparison_allPower/weighted_ubiq_category_crispr_results_combined_percent_positives.tsv",
    fig3c = "results/main_figure_5_and_6/housekeeping_genes/detected_positives/5_10_positives_by_effect_size_and_rate.Regulated.pdf"
  log: "results/main_figure_5_and_6/logs/main_figure_5_chromatin_categories.log"
  shell:
    "Rscript {input.script} plot results/main_figure_5_and_6 > {log} 2>&1"

# Figure 5c,d (DC-TAP-seq split and combined by cell type), from Table S3
rule main_figure_5_housekeeping_genes:
  input:
    script = "results/main_figure_5_and_6/housekeeping_genes.R",
    functions = "results/main_figure_5_and_6/housekeeping_gene_plot_functions.R",
    table_s3 = "results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv"
  output:
    significance = "results/main_figure_5_and_6/hkg_comparison/combinedDC_TAP_all_metrics_significance_allPower/weighted_ubiq_category.signficance_table.by_Category.tsv"
  log: "results/main_figure_5_and_6/logs/main_figure_5_housekeeping_genes.log"
  shell:
    "Rscript {input.script} {input.table_s3} results/main_figure_5_and_6 > {log} 2>&1"

# Figure 6 (ENCODE-rE2G) and Extended Data Fig. 13 (ABC, AlphaGenome), from CRISPR_comparison output
rule main_figure_6_model_scores:
  input:
    script = "results/main_figure_5_and_6/model_scores.R",
    crispr_benchmark = "resources/main_figure_6/expt_pred_merged_annot.txt.gz"
  output:
    fig6a = "results/main_figure_5_and_6/e2g_model_scores/ENCODE_rE2G/pDirect_category_by_score_Regulated.pdf",
    fig6b = "results/main_figure_5_and_6/e2g_model_scores/ENCODE_rE2G/filter50pct_EffectSize_5pct_by_score_Regulated.pdf",
    fig6c = "results/main_figure_5_and_6/e2g_model_scores/ENCODE_rE2G/filter50pct_element_category_by_score_Regulated.pdf",
    fig6d = "results/main_figure_5_and_6/e2g_model_scores/ENCODE_rE2G/filter50pct_ubiq_category_by_score_Regulated.pdf"
  log: "results/main_figure_5_and_6/logs/main_figure_6_model_scores.log"
  shell:
    "Rscript {input.script} {input.crispr_benchmark} results/main_figure_5_and_6 > {log} 2>&1"
