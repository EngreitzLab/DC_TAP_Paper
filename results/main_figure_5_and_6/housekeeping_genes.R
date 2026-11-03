# Fig. 5c,d (DC-TAP-seq, split and combined by cell type) from Table S3
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "housekeeping_gene_plot_functions.R"))

### MAIN
# usage: Rscript housekeeping_genes.R [table_s3_file] [results_dir]  (run from the repository root)
args <- commandArgs(trailingOnly = TRUE)
table_s3_file <- if (length(args) >= 1) args[1] else "results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv"
results_dir <- if (length(args) >= 2) args[2] else "results/main_figure_5_and_6"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

cell_types_cr <- c("GM12878", "HCT116", "Jurkat", "K562_training", "K562_validation", "WTC11")#, "H1_for_WTC11", "H9_for_WTC11")
ct_cp_cr <- c("#e9c54e", "#f29742", "#5eb342", "#792374", "#b778b3", "#5496ce")#, "#00488d", "#00488d")
ds_cp <- c(Gasperini2019 = "#d3a9ce", Nasser2021 = "#b778b3", Schraivogel2020 = "#a64791", K562_DC_TAP = "#006eae", WTC11_DC_TAP = "#00488d")

names(ct_cp_cr) <- cell_types_cr

# --- plot housekeeping gene comparisons, split by cell type --- #
if (TRUE) {
    crispr_merged <- read_dc_tap_table_s3(table_s3_file, filter_to_random_DEG = TRUE) %>% 
        #mutate(Dataset = ifelse(Dataset %in% c("K562_DC_TAP", "WTC11_DC_TAP"), "DC_TAP", Dataset)) %>% 
        mutate(enhancerness = ifelse(element_category %in% c("High H3K27ac", "H3K27ac"), "H3K27ac+ element", "Other element")) %>% 
        mutate(EffectSize = EffectSize * 100) %>% 
        mutate(direct_vs_indirect = ifelse(EffectSize < 0, direct_vs_indirect_negative, direct_vs_indirect_positive))

    group_vars <- c("ubiq_category")
    this_results_dir <- file.path(results_dir, "hkg_comparison"); dir.create(this_results_dir, showWarnings = FALSE)
    
    for (g in group_vars) {
        for (use_allpower in c(TRUE)){
            if (use_allpower) {
                out_dir <- file.path(this_results_dir, "splitDC_TAP_all_metrics_comparison_allPower"); dir.create(out_dir, showWarnings = FALSE)
                out_dir2 <- file.path(this_results_dir, "splitDC_TAP_all_metrics_significance_allPower"); dir.create(out_dir2, showWarnings = FALSE)
            } else {
                out_dir <- file.path(this_results_dir, "splitDC_TAP_all_metrics_comparison_wellPowered"); dir.create(out_dir, showWarnings = FALSE)
                out_dir2 <- file.path(this_results_dir, "splitDC_TAP_all_metrics_significance_wellPowered"); dir.create(out_dir2, showWarnings = FALSE)
            }

            # direct effect weighted
            out_prefix <- paste0(out_dir, "/weighted_")
            out_prefix2 <- paste0(out_dir2, "/weighted_")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = TRUE, direct_effect_threshold = NULL)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = TRUE, direct_effect_threshold = NULL)

            # direct effect filtered
            out_prefix <- paste0(out_dir, "/filter50_")
            out_prefix2 <- paste0(out_dir2, "/filter50_")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)

            # no adjustments for p(direct)
            out_prefix <- paste0(out_dir, "/")
            out_prefix2 <- paste0(out_dir2, "/")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0)
        }
    }
}

# --- plot housekeeping gene comparisons, combined cell type --- #
if (TRUE) {
    crispr_merged <- read_dc_tap_table_s3(table_s3_file, filter_to_random_DEG = TRUE) %>% 
        mutate(Dataset = ifelse(Dataset %in% c("K562_DC_TAP", "WTC11_DC_TAP"), "DC_TAP", Dataset)) %>% 
        mutate(enhancerness = ifelse(element_category %in% c("High H3K27ac", "H3K27ac"), "H3K27ac+ element", "Other element")) %>% 
        mutate(EffectSize = EffectSize * 100) %>% 
        mutate(direct_vs_indirect = ifelse(EffectSize < 0, direct_vs_indirect_negative, direct_vs_indirect_positive))

    group_vars <- c("ubiq_category")
    this_results_dir <- file.path(results_dir, "hkg_comparison"); dir.create(this_results_dir, showWarnings = FALSE)
    
    for (g in group_vars) {
        for (use_allpower in c(TRUE)){
            if (use_allpower) {
                out_dir <- file.path(this_results_dir, "combinedDC_TAP_all_metrics_comparison_allPower"); dir.create(out_dir, showWarnings = FALSE)
                out_dir2 <- file.path(this_results_dir, "combinedDC_TAP_all_metrics_significance_allPower"); dir.create(out_dir2, showWarnings = FALSE)
            } else {
                out_dir <- file.path(this_results_dir, "combinedDC_TAP_all_metrics_comparison_wellPowered"); dir.create(out_dir, showWarnings = FALSE)
                out_dir2 <- file.path(this_results_dir, "combinedDC_TAP_all_metrics_significance_wellPowered"); dir.create(out_dir2, showWarnings = FALSE)
            }

            # direct effect weighted
            out_prefix <- paste0(out_dir, "/weighted_")
            out_prefix2 <- paste0(out_dir2, "/weighted_")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = TRUE, direct_effect_threshold = NULL)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = TRUE, direct_effect_threshold = NULL)

            # direct effect filtered
            out_prefix <- paste0(out_dir, "/filter50_")
            out_prefix2 <- paste0(out_dir2, "/filter50_")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)

            # no adjustments for p(direct)
            out_prefix <- paste0(out_dir, "/")
            out_prefix2 <- paste0(out_dir2, "/")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0)
        }
    }
}