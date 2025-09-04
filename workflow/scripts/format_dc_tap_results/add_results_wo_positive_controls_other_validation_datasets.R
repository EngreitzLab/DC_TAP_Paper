# Script: add_results_wo_positive_controls_other_validation_datasets.R

### SETUP =====================================================================

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
element_gene_pairs <- read_tsv(snakemake@input$validation_dataset_power_analysis_results)

# Read in power simulation data for ALL effect sizes
effect_sizes <- c(2, 3, 5, 10, 15, 20, 25, 50)

message("Loading combined power analysis files")
# Load each validation dataset power simulations for all effect sizes using map_dfr
validations_power_simulation <- map_dfr(seq_along(effect_sizes), function(i) {
    fread(snakemake@input$combined_power_analysis_output_validation_dataset[[i]]) %>%
      as_tibble() %>%
      mutate(effect_size = effect_sizes[i])
})

### CUSTOM FDR CORRECTION WITHOUT POSITIVE CONTROLS ==========================

# Note: Refactor this script or divide the responsibility of this script 
#       as there could be a case where the input validation dataset 
#       had not been pre-filtered of it's selfPromoter Pairs.

message("Applying custom FDR correction excluding positive controls")

# Apply FDR correction by cell type, excluding positive controls
final_pairs <- element_gene_pairs %>%
  mutate(
    # Create a flag for pairs to include in FDR correction
    # Note, for validation datasets pulled from ENCODE test Dataset Analysis pipeline
    # selfPromoter DE-G pairs have already been filtered
    include_in_fdr = !is.na(pValue),
    
    # Apply Benjamini-Hochberg correction only to the subset
    sceptre_adj_p_value_wo_pos_controls = ifelse(
      include_in_fdr & !is.na(pValue),
      p.adjust(ifelse(include_in_fdr & !is.na(pValue), pValue, NA), method = "BH"),
      NA_real_
    ),
    
    # Call significance at FDR specified padj_threshold (current default at 0.2)
    significant_wo_pos_controls = case_when(
      include_in_fdr & !is.na(sceptre_adj_p_value_wo_pos_controls) ~ sceptre_adj_p_value_wo_pos_controls <= snakemake@params$padj_threshold,
      TRUE ~ NA
    )
  ) %>%
  ungroup()

# Function to calculate power for a given effect size
calculate_power_for_effect_size <- function(df, effect_size_val, max_p_val) {
  df %>%
    filter(effect_size == effect_size_val, !is.na(log_2_fold_change)) %>%
    group_by(grna_target, response_id) %>%
    summarize(
      !!paste0("power_at_effect_size_", effect_size_val, "_wo_pos_controls") := 
        mean(p_value < max_p_val & log_2_fold_change < 0),
      .groups = "drop"
    )
}
recalculate_power_for_effect_size <- function(df, effect_size_val, max_p_val) {
  df %>%
    filter(effect_size == effect_size_val, !is.na(log_2_fold_change)) %>%
    group_by(grna_target, response_id) %>%
    summarize(
      !!paste0("power_at_effect_size_", effect_size_val) := 
        mean(p_value < max_p_val & log_2_fold_change < 0),
      .groups = "drop"
    )
}

### RECALCULATE POWER FOR WITH_POS_CONTROLS FOR ALL EFFECT SIZES ===============

message("Recalculating power for with_pos_controls significance threshold for all effect sizes")

# Update Significant calls based on specified FDR threshold
final_pairs <- final_pairs %>%
  mutate(Significant = pValueAdjusted <= snakemake@params$padj_threshold)

# Get max nominal p-value from with_pos_controls significant pairs
max_nom_p_val_validations_recal <- final_pairs %>% 
  filter(Significant == TRUE) %>% 
  pull(pValue) %>% 
  max(na.rm = TRUE)

message(paste("Max nominal p-value for validation dataset with_pos_controls:", max_nom_p_val_validations_recal))

# Calculate power for all effect sizes for Validation Dataset
power_wo_pos_controls_validations_recal <- map(effect_sizes, function(es) {
  recalculate_power_for_effect_size(validations_power_simulation, es, max_nom_p_val_validations_recal)
}) %>%
  reduce(left_join, by = c("grna_target", "response_id"))

# Merge both power calculations back to final_pairs
final_pairs <- final_pairs %>%
  # Format PerturbationTargetID to remove ":." artifact. 
  mutate(PerturbationTargetID = sub(":\\.$", "", PerturbationTargetID)) %>%
  # Drop existing power calculations
  select(-matches("^PowerAtEffectSize")) %>%
  left_join(
    power_wo_pos_controls_validations_recal,
    by = c("PerturbationTargetID" = "grna_target", "measuredEnsemblID" = "response_id")
  )

### RECALCULATE POWER FOR WO_POS_CONTROLS FOR ALL EFFECT SIZES ===============

message("Recalculating power for wo_pos_controls significance threshold for all effect sizes")

# Get max nominal p-value from wo_pos_controls significant pairs
max_nom_p_val_validations <- final_pairs %>% 
  filter(significant_wo_pos_controls == TRUE) %>% 
  pull(pValue) %>% 
  max(na.rm = TRUE)

message(paste("Max nominal p-value for validation dataset wo_pos_controls:", max_nom_p_val_validations))

# Calculate power for all effect sizes for Validation Dataset
power_wo_pos_controls_validations <- map(effect_sizes, function(es) {
  calculate_power_for_effect_size(validations_power_simulation, es, max_nom_p_val_validations)
}) %>%
  reduce(left_join, by = c("grna_target", "response_id"))

# Calculate mean perturbation cells across all effect sizes
mean_pert_cells_summary <- validations_power_simulation %>%
  group_by(grna_target, response_id) %>%
  summarize(
    mean_sim_pert_cells = mean(num_pert_cells, na.rm = TRUE),
    .groups = "drop"
  )

# Merge both power calculations and mean_pert_cells back to final_pairs
final_pairs <- final_pairs %>%
  left_join(
    power_wo_pos_controls_validations,
    by = c("PerturbationTargetID" = "grna_target", "measuredEnsemblID" = "response_id")
  ) %>%
  left_join(
    mean_pert_cells_summary,
    by = c("PerturbationTargetID" = "grna_target", "measuredEnsemblID" = "response_id")
  )


### MAKE POWER_WO_POS_CONTROLS NA FOR ALL EFFECT SIZES =======================

# Power at all effect sizes without positive controls should be NA where significance wasn't calculated
# Because Positive controls weren't included in the FDR correction
for(es in effect_sizes) {
  power_col <- paste0("power_at_effect_size_", es, "_wo_pos_controls")
  final_pairs[[power_col]] <- ifelse(is.na(final_pairs$significant_wo_pos_controls), 
                                     NA, 
                                     final_pairs[[power_col]])
}

# Set mean_sim_pert_cells of pairs which power isn't calculated for to NA
final_pairs <- final_pairs %>% 
  mutate(mean_sim_pert_cells = ifelse(is.na(significant_wo_pos_controls), NA, mean_sim_pert_cells))


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(final_pairs, snakemake@output$results_wo_pos_controls)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)