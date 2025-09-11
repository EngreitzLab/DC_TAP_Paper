# Script: add_results_wo_positive_controls_gasperini.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/add_results_wo_positive_controls_gasperini.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

message("Loading input files")
element_gene_pairs <- read_tsv(snakemake@input$gasperini_results)

# Read in power simulation data for ALL effect sizes
effect_sizes <- c(2, 3, 5, 10, 15, 20, 25, 50)

# Load Gasperini power simulations for all effect sizes using map_dfr
gasperini_power_simulation <- map_dfr(seq_along(effect_sizes), function(i) {
  fread(snakemake@input$combined_power_analysis_output_gasperini[[i]]) %>%
    as_tibble() %>%
    mutate(cell_type = "K562", effect_size = effect_sizes[i])
})

### CUSTOM FDR CORRECTION WITHOUT POSITIVE CONTROLS ==========================

message("Applying custom FDR correction excluding positive controls")

# Apply FDR correction by cell type, excluding positive controls
final_pairs <- element_gene_pairs %>%
  group_by(cell_type) %>%
  mutate(
    # Create a flag for pairs to include in FDR correction
    include_in_fdr = (selfPromoter == FALSE),
    
    # Apply Benjamini-Hochberg correction only to the subset
    sceptre_adj_p_value_wo_pos_controls = ifelse(
      include_in_fdr & !is.na(sceptre_p_value),
      p.adjust(ifelse(include_in_fdr & !is.na(sceptre_p_value), sceptre_p_value, NA), method = "BH"),
      NA_real_
    ),
    
    # Call significance at FDR specified padj_threshold (current default at 0.2)
    significant_wo_pos_controls_20fdr = case_when(
      include_in_fdr & !is.na(sceptre_adj_p_value_wo_pos_controls) ~ sceptre_adj_p_value_wo_pos_controls <= snakemake@params$padj_threshold,
      TRUE ~ NA
    )
  ) %>%
  ungroup()


### RECALCULATE POWER FOR WO_POS_CONTROLS FOR ALL EFFECT SIZES ===============

message("Recalculating power for wo_pos_controls significance threshold for all effect sizes")

# Get max nominal p-value from wo_pos_controls significant pairs for each cell type
max_nom_p_val_gasperini <- final_pairs %>% 
  filter(significant_wo_pos_controls_20fdr == TRUE, cell_type == "K562") %>% 
  pull(sceptre_p_value) %>% 
  max(na.rm = TRUE)

message(paste("Max nominal p-value for Gasperini_K562 wo_pos_controls:", max_nom_p_val_gasperini))

# Function to calculate power for a given effect size
calculate_power_for_effect_size <- function(df, effect_size_val, max_p_val, cell_type_val) {
  df %>%
    filter(cell_type == cell_type_val, effect_size == effect_size_val, !is.na(log_2_fold_change)) %>%
    group_by(grna_target, response_id) %>%
    summarize(
      !!paste0("power_at_effect_size_", effect_size_val, "_wo_pos_controls_20fdr") := 
        mean(p_value < max_p_val & log_2_fold_change < 0),
      .groups = "drop"
    )
}

# Calculate power for all effect sizes for Gasperini K562
power_wo_pos_controls_gasperini <- map(effect_sizes, function(es) {
  calculate_power_for_effect_size(gasperini_power_simulation, es, max_nom_p_val_gasperini, "K562")
}) %>%
  reduce(left_join, by = c("grna_target", "response_id")) %>%
  mutate(cell_type = "K562")

# Calculate mean perturbation cells across all effect sizes for both cell types
mean_pert_cells_summary <- gasperini_power_simulation %>%
  group_by(grna_target, response_id, cell_type) %>%
  summarize(
    mean_sim_pert_cells = mean(num_pert_cells, na.rm = TRUE),
    .groups = "drop"
  )

# Merge both power calculations and mean_pert_cells back to final_pairs
final_pairs <- final_pairs %>%
  left_join(
    power_wo_pos_controls_gasperini,
    by = c("design_file_target_name" = "grna_target", "gene_id" = "response_id", "cell_type")
  ) %>%
  left_join(
    mean_pert_cells_summary,
    by = c("design_file_target_name" = "grna_target", "gene_id" = "response_id", "cell_type")
  )


### MAKE POWER_WO_POS_CONTROLS NA FOR ALL EFFECT SIZES =======================

# Power at all effect sizes without positive controls should be NA where significance wasn't calculated
# Because Positive controls weren't included in the FDR correction
for(es in effect_sizes) {
  power_col <- paste0("power_at_effect_size_", es, "_wo_pos_controls_20fdr")
  final_pairs[[power_col]] <- ifelse(is.na(final_pairs$significant_wo_pos_controls_20fdr), 
                                     NA, 
                                     final_pairs[[power_col]])
}

# Set mean_sim_pert_cells of pairs which power isn't calculated for to NA
final_pairs <- final_pairs %>% 
  mutate(mean_sim_pert_cells = ifelse(is.na(significant_wo_pos_controls_20fdr), NA, mean_sim_pert_cells)) %>%
  mutate(
  # Recall significant_with_positive controls based on FDR 10.
  significant = case_when(
    !is.na(sceptre_adj_p_value) ~ sceptre_adj_p_value <= 0.10,
    TRUE ~ NA
  ))


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(final_pairs, snakemake@output$results_wo_pos_controls)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)