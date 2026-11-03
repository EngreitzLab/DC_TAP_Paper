## Annotate cis-results with element and pair annotations for DC-TAP-seq experiments

# save.image("pairs.rda")
# stop()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# load annotated results
results <- fread(snakemake@input[[1]])

# only retain data on the given sample
sample <- sub("_DC_TAPseq", "", snakemake@wildcards$sample)
results <- filter(results, cell_type == sample)

# add valid connection columns and distance to TSS to cis-results
cis_results <- results %>% 
  filter(include_in_fdr == TRUE) %>% 
  mutate(pair_uid = paste0(gene_id, "|", design_file_target_name), n_nonzero_trt = NA_integer_,
         n_nonzero_cntrl = NA_integer_, pass_qc = TRUE) %>% 
  select(response_id = gene_id, grna_target = design_file_target_name, n_nonzero_trt,
         n_nonzero_cntrl, pass_qc, p_value = sceptre_p_value,
         log_2_fold_change = log_2_FC_effect_size, significant = significant_wo_pos_controls_20fdr,
         pair_uid, gene_symbol, dist_to_tss = distance_to_gencode_gene_TSS, element_location,
         valid_connection = all_of(snakemake@params$valid_col))

# get all valid elements that do not overlap any annotated promoters
valid_elements <- results %>% 
  filter(include_in_fdr == TRUE, Random_DistalElement_Gene == TRUE) %>%
  pull(design_file_target_name) %>% 
  unique()

# create output with column labeling valid elements
output <- mutate(cis_results, valid_element = grna_target %in% valid_elements)

# save output to file
fwrite(output, file = snakemake@output[[1]], sep = ",", quote = FALSE, na = "NA")
