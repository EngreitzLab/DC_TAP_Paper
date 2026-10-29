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
    genome_wide_pairs = "results/main_figure_5_and_6/categorized_pairs/all_gw_pairs.categorized.summary.tsv",
    table_s3 = "results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv",
    gasperini_results = "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_element_gene_pair_categories_wo_pos_controls.tsv"
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

# Figure 2e: DC-TAP-seq vs qPCR effect sizes for selected gRNAs (K562 and WTC11)
rule main_figure_2e_dctap_vs_qpcr:
  input:
    script = "workflow/scripts/main_figure_1_and_2/plot_dctap_vs_qpcr_effect_sizes.R",
    effect_sizes = "resources/main_figure_1_and_2/qpcr_validation/k562_wtc11_qpcr_vs_dctap_effect_sizes.tsv"
  output:
    "results/main_figure_1_and_2/qpcr_validation/dctap_vs_qpcr_effect_sizes.pdf"
  log: "results/main_figure_1_and_2/logs/main_figure_2e_dctap_vs_qpcr.log"
  shell:
    "Rscript {input.script} {input.effect_sizes} {output} > {log} 2>&1"

# Figure 1c,d: guide capture efficiency of CROP-seq vs direct-capture vectors, and K562 knockdown vs MOI
rule main_figure_1cd_guide_capture_and_moi:
  input:
    script = "workflow/scripts/main_figure_1_and_2/plot_guide_capture_and_moi_knockdown.R",
    crop_guides = "resources/main_figure_1_and_2/guide_capture_and_moi_knockdown/P10_tap_all_guide_info.txt",
    hu6_calls = "resources/main_figure_1_and_2/guide_capture_and_moi_knockdown/protospacer_calls_per_cell_052322_analysis_L-SG-TAP-Whole.csv",
    mu6_calls = "resources/main_figure_1_and_2/guide_capture_and_moi_knockdown/protospacer_calls_per_cell_052322_analysis_H-DC-TAP-Whole.csv",
    moi_knockdown = "resources/main_figure_1_and_2/guide_capture_and_moi_knockdown/input_file_for_k562_multi_moi1_panel.txt"
  output:
    fig1c = "results/main_figure_1_and_2/guide_capture_and_moi_knockdown/fig1c_guide_umi_capture.pdf",
    fig1d = "results/main_figure_1_and_2/guide_capture_and_moi_knockdown/fig1d_moi_knockdown.pdf"
  log: "results/main_figure_1_and_2/logs/main_figure_1cd_guide_capture_and_moi.log"
  params:
    input_dir = "resources/main_figure_1_and_2/guide_capture_and_moi_knockdown"
  shell:
    "Rscript {input.script} {params.input_dir} {output.fig1c} {output.fig1d} > {log} 2>&1"

# Figure 2i: SLC2A3 qPCR for guide #5742 vs non-targeting (values are in the script)
rule main_figure_2i_slc2a3_qpcr:
  input:
    script = "workflow/scripts/main_figure_1_and_2/plot_slc2a3_qpcr.R"
  output:
    "results/main_figure_1_and_2/qpcr_validation/slc2a3_qpcr_guide5742.pdf"
  log: "results/main_figure_1_and_2/logs/main_figure_2i_slc2a3_qpcr.log"
  shell:
    "Rscript {input.script} {output} > {log} 2>&1"

# GENCODE v32 gene annotation, used by the validation-screen analyses to map gene IDs to symbols
rule download_gencode_v32_for_validations:
  output:
    "validations/data/gencode.v32.primary_assembly.annotation.gtf.gz"
  shell:
    "curl -fsSL --retry 3 -o {output} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"

# Figure 2d (top): low-MOI validation screen vs full screen effect sizes (E. Mattei's compare_es.Rmd).
# Also writes the high-MOI comparisons. Runs with --vanilla to skip the renv set-up in validations/.Rprofile.
rule main_figure_2d_validation_vs_screen:
  input:
    rmd = "validations/analyses/2026-05-14_compare_low_high_es/compare_es.Rmd",
    low_moi = "validations/analyses/2026-05-11_SCEPTRE/source_data/runs/low_moi_merged_union/sceptre_discovery_results.tsv",
    high_moi = "validations/analyses/2026-05-11_SCEPTRE/source_data/runs/high_moi_merged_union/sceptre_discovery_results.tsv",
    validation_pairs = "validations/metadata/validation_30_pairs.tsv",
    gencode = "validations/data/gencode.v32.primary_assembly.annotation.gtf.gz",
    table_s3 = "results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv"
  output:
    fig2d_top = "validations/results/2026-05-14_compare_low_high_es/scatter_effect_size_correlation_low_vs_paper.pdf",
    high_vs_paper = "validations/results/2026-05-14_compare_low_high_es/scatter_effect_size_correlation_high_vs_paper.pdf",
    high_vs_low = "validations/results/2026-05-14_compare_low_high_es/scatter_effect_size_correlation_high_vs_low.pdf"
  log: "results/main_figure_1_and_2/logs/main_figure_2d_validation_vs_screen.log"
  shell:
    """
    log=$(pwd)/{log}
    cd validations && Rscript --vanilla -e "rmarkdown::render('analyses/2026-05-14_compare_low_high_es/compare_es.Rmd', output_dir = 'results/2026-05-14_compare_low_high_es', knit_root_dir = getwd())" > $log 2>&1
    """

# Figure 1e: UMI saturation of DC-TAP-seq vs scRNA-seq with guide capture (J. Ray's read-downsampling tables)
rule main_figure_1e_umi_saturation:
  input:
    script = "workflow/scripts/main_figure_1_and_2/plot_umi_saturation_dctap_vs_scrna.R",
    dc_tap = "resources/main_figure_1_and_2/umi_saturation_dctap_vs_scrna/complexity_df_dctap_moi3.csv",
    scrna = "resources/main_figure_1_and_2/umi_saturation_dctap_vs_scrna/complexity_df_smaller_scrna_w_guide.csv"
  output:
    "results/main_figure_1_and_2/umi_saturation_dctap_vs_scrna/fig1e_umi_saturation.pdf"
  log: "results/main_figure_1_and_2/logs/main_figure_1e_umi_saturation.log"
  params:
    input_dir = "resources/main_figure_1_and_2/umi_saturation_dctap_vs_scrna"
  shell:
    "Rscript {input.script} {params.input_dir} {output} > {log} 2>&1"
