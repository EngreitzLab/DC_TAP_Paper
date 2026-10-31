## Combine cis and trans positive hit rates from different samples into one dataframe

# save.image("combine_rates.rda")
# stop()

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Input files --------------------------------------------------------------------------------------

# all cis positive hit rates for held-out datasets
cis_rate_heldout_files <- snakemake@input$cis_rates_heldout
names(cis_rate_heldout_files) <- basename(dirname(cis_rate_heldout_files))

# all trans positive hit rates for held-out datasets
trans_rate_heldout_files <- snakemake@input$trans_rates_heldout
names(trans_rate_heldout_files) <- basename(dirname(trans_rate_heldout_files))

# all cis positive hit rates for training datasets
cis_rate_training_files <- snakemake@input$cis_rates_training
names(cis_rate_training_files) <- basename(dirname(dirname(cis_rate_training_files)))

# all trans positive hit rates for training datasets
trans_rate_training_files <- snakemake@input$trans_rates_training
names(trans_rate_training_files) <- basename(dirname(dirname(trans_rate_training_files)))

# Load and combine files ---------------------------------------------------------------------------

# load cis positive hit rates for held-out data and combine into one table
cis_rate_heldout <- cis_rate_heldout_files %>% 
  lapply(FUN = read_tsv, show_col_types = FALSE) %>% 
  bind_rows(.id = "dataset")

# load trans positive hit rates for held-out data and combine into one table
trans_rate_heldout <- trans_rate_heldout_files %>% 
  lapply(FUN = read_tsv, show_col_types = FALSE) %>% 
  bind_rows(.id = "dataset")

# load cis positive hit rates for training data and combine into one table
cis_rate_training <- cis_rate_training_files %>% 
  lapply(FUN = read_tsv, show_col_types = FALSE) %>% 
  bind_rows(.id = "dataset")

# load trans positive hit rates for training data and combine into one table
trans_rate_training <- trans_rate_training_files %>% 
  lapply(FUN = read_tsv, show_col_types = FALSE) %>% 
  bind_rows(.id = "dataset")

# combine cis and trans rates into two tables
cis_rate <- bind_rows(heldout = cis_rate_heldout, training = cis_rate_training, .id = "type")
trans_rate <- bind_rows(heldout = trans_rate_heldout, training = trans_rate_training,
                         .id = "type")

# save combined tables to output files
write_tsv(cis_rate, file = snakemake@output$cis_rates)
write_tsv(trans_rate, file = snakemake@output$trans_rates)
