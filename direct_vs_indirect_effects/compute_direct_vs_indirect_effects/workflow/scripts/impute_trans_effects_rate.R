## Impute trans positive hit rates for CRISPR datasets where analysis couldn't be performed 

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# load trans positive hit rates and list of CRISPR datasets
trans_hit_rates <- read_tsv(snakemake@input$trans_rates, show_col_types = FALSE)
crispr_datasets <- read_tsv(snakemake@input$dataset_ids, show_col_types = FALSE)

# calculate average trans hit rates across all datasets
avg_trans_hit_rates <- trans_hit_rates %>% 
  summarize(positive_rate_significant = mean(positive_rate_significant),
            positive_rate_negative = mean(positive_rate_negative),
            positive_rate_positive = mean(positive_rate_positive))

# get all CRISPR datasets for which trans hit rate wasn't computed and needs to be imputed
impute_datasets <- setdiff(crispr_datasets$dataset, trans_hit_rates$dataset)

# create table with imputed trans hit rates for these datasets
imputed_trans_hit_rates <- tibble(dataset = impute_datasets, imputed_indirect_rate = TRUE) %>% 
  bind_cols(avg_trans_hit_rates)

# combine with calculated trans hit rates table
trans_hit_rates <- trans_hit_rates %>% 
  mutate(imputed_indirect_rate = FALSE) %>% 
  bind_rows(imputed_trans_hit_rates)

# save table with imputed hit rates to output file
write_tsv(trans_hit_rates, file = snakemake@output[[1]])
