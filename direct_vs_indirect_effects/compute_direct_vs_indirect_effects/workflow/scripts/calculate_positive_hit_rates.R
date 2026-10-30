## Calculate cis and trans positive hit rates for a sample. Cis positive rate is calculated across
## specified distance bins

# save.image("calc_rate.rda")
# stop()

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# load cis and trans analysis results
trans_results <- readRDS(snakemake@input$trans_results)
cis_results <- read_csv(snakemake@input$cis_results, show_col_types = FALSE)

# load gene universe file
gene_univ <- read_tsv(snakemake@input$gene_univ, show_col_types = FALSE)

# significance threshold for cis effects
pval_threshold_cis <- cis_results %>% 
  filter(significant == TRUE) %>% 
  pull(p_value) %>% 
  max()

# add significance and regulated column to trans-effect results based on cis significance threshold
trans_results <- trans_results %>% 
  mutate(significant = p_value <= pval_threshold_cis,
         regulated_negative = significant & log_2_fold_change < 0,
         regulated_positive = significant & log_2_fold_change >= 0)

# add regulated columns to cis results
cis_results <- cis_results %>% 
  mutate(regulated_negative = significant & log_2_fold_change < 0,
         regulated_positive = significant & log_2_fold_change >= 0)

# get all valid elements from cis results
valid_elements <- cis_results %>% 
  filter(valid_element == TRUE) %>% 
  pull(grna_target) %>% 
  unique()

# add a valid connection column to trans-results based on whether perturbation is a valid element
trans_results <- mutate(trans_results, valid_connection = grna_target %in% valid_elements)

# filter cis and trans results based on whether they are valid connection for the analysis
cis_results_filt <- filter(cis_results, valid_connection == TRUE, pass_qc == TRUE)
trans_results_filt <- filter(trans_results, valid_connection == TRUE, pass_qc == TRUE)

# only retain data on genes in distal regulation genes universe
cis_results_filt <- filter(cis_results_filt, response_id %in% gene_univ$Ensembl_ID)
trans_results_filt <- filter(trans_results_filt, response_id %in% gene_univ$Ensembl_ID)

# distance bins for calculating positive hit rate
cis_results_filt <- mutate(cis_results_filt, abs_dist_to_tss = abs(dist_to_tss))
max_dist <- ceiling(max(cis_results_filt$abs_dist_to_tss, na.rm = TRUE) / 1e6) * 1e6
dist_bins <- seq(0, max_dist, by = snakemake@params$bin_size)

# bin cis pairs by distance
cis_results_filt <- cis_results_filt %>% 
  mutate(dist_bin = cut(abs_dist_to_tss, breaks = dist_bins, include.lowest = TRUE))

# calculate positive rates and mean distance per distance bin
cis_pos_rate <- cis_results_filt %>% 
  group_by(dist_bin) %>% 
  summarize(total_pairs = n(),
            significant_pairs = sum(significant, na.rm = TRUE),
            negative_pairs = sum(regulated_negative, na.rm = TRUE),
            positive_pairs = sum(regulated_positive, na.rm = TRUE),
            positive_rate_significant = significant_pairs / total_pairs,
            positive_rate_negative = negative_pairs / total_pairs,
            positive_rate_positive = positive_pairs / total_pairs,
            mean_dist = mean(abs_dist_to_tss))

# calculate trans-effects positive rate
trans_pos_rate <- trans_results_filt %>% 
  summarize(total_pairs = n(),
            significant_pairs = sum(significant, na.rm = TRUE),
            negative_pairs = sum(regulated_negative, na.rm = TRUE),
            positive_pairs = sum(regulated_positive, na.rm = TRUE),
            positive_rate_significant = significant_pairs / total_pairs,
            positive_rate_negative = negative_pairs / total_pairs,
            positive_rate_positive = positive_pairs / total_pairs)

# save output to files
write_tsv(cis_pos_rate, snakemake@output$cis_rate)
write_tsv(trans_pos_rate, snakemake@output$trans_rate)
