## Fit powerlaw models to model direct effects hit rate as function of distance to TSS

# save.image("model_direct_effects_rate.rda")
# stop()

# attach required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

## Define functions --------------------------------------------------------------------------------

# fit a powerlaw to model the positive hit ratio as a function of distance to TSS
fit_powerlaw_model <- function(direct_rates, rate_type) {
  
  # get direct rates for given type and remove 0s
  direct_rates_type <- filter(direct_rates, type == rate_type)
  direct_rates_model <- filter(direct_rates_type, direct_rate  > 0)
  message("Discarding ", nrow(direct_rates_type) - nrow(direct_rates_model),
          " distance bins with a direct rate of 0")
  
  # fit a powerlaw using a linear model
  fit <- lm(log(direct_rate) ~ log(dist_to_tss), data = direct_rates_model)

  # extract fitted coefficients and report fit
  a <- exp(coef(fit)[1])
  b <- coef(fit)[2]
  message("Fitted model: direct_rate = ", round(a, 4), " * dist_to_tss^", round(b, 4))
  
  return(fit)
  
}

## For models for direct positive rate as function of distance to TSS ------------------------------

# set random seed for any RNG processes
set.seed(snakemake@params$seed)

# load cis and trans rates
cis_rates <- read_tsv(snakemake@input$cis_rates, show_col_types = FALSE)
trans_rates <- read_tsv(snakemake@input$trans_rates, show_col_types = FALSE)

# filter cis positive rate data to within specified distance
cis_rates <- filter(cis_rates, mean_dist <= snakemake@params$max_distance)

# add expected indirect effect rates to cis positive rates
cis_rates <- trans_rates %>% 
  select(dataset,
         indirect_rate_all = positive_rate_significant,
         indirect_rate_negative = positive_rate_negative,
         indirect_rate_positive = positive_rate_positive) %>% 
  left_join(cis_rates, ., by = "dataset")

# calculate expected direct effects rate by subtracting the indirect rate from the cis rates
cis_rates <- cis_rates %>% 
  rowwise() %>%
  mutate(direct_rate_all = max(positive_rate_significant - indirect_rate_all, 0),
         direct_rate_negative = max(positive_rate_negative - indirect_rate_negative, 0),
         direct_rate_positive = max(positive_rate_positive - indirect_rate_positive, 0))

# get direct rates per dataset in long format
direct_rates_datasets <- cis_rates %>% 
  select(dataset, dist_bin, dist_to_tss = mean_dist, direct_rate_all, direct_rate_negative,
         direct_rate_positive) %>% 
  pivot_longer(cols = -c(dataset, dist_bin, dist_to_tss), names_to = "type",
               values_to = "direct_rate")

# average direct rates across datasets to mitigate zeros and convert to long format
direct_rates_avg <- cis_rates %>% 
  group_by(dist_bin) %>% 
  summarize(direct_rate_all = mean(direct_rate_all),
            direct_rate_negative = mean(direct_rate_negative),
            direct_rate_positive = mean(direct_rate_positive),
            dist_to_tss = mean(mean_dist)) %>% 
  pivot_longer(cols = -c(dist_bin, dist_to_tss), names_to = "type", values_to = "direct_rate")

# fit powerlaw models for each type of direct rate
rate_types <- structure(unique(direct_rates_avg$type), names = unique(direct_rates_avg$type))
models <- lapply(rate_types, FUN = fit_powerlaw_model, direct_rates = direct_rates_avg)

# save direct rates to files for downstream analyses and plots
write_tsv(direct_rates_datasets, file = snakemake@output$direct_rates_datasets)
write_tsv(direct_rates_avg, file = snakemake@output$direct_rates_average)

# save all models to .rds file to later predict direct effect rates as a function of distance to TSS
saveRDS(models, file = snakemake@output$models)
