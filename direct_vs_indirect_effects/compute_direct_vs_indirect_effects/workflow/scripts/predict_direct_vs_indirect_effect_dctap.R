## Calculate probability of CRISPR hits to be direct vs. indirect effects based on predicted direct
## effects rate dependent on distance to TSS and estimated indirect effects rate

# save.image("direct_vs_indirect_dctap.rda")
# stop()

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(GenomicRanges)
})

## Define functions --------------------------------------------------------------------------------

# compute distance to TSS
compute_dist_tss <- function(crispr, tss_annot) {
  
  # only retain crispr E-G pairs with genes in TSS universe
  missing <- setdiff(crispr$measuredGeneSymbol, tss_annot$name)
  if (length(missing) > 0) {
    message("Filtering out E-G pairs for ", length(missing), " genes not in TSS annotations")
    crispr <- filter(crispr, measuredGeneSymbol %in% tss_annot$name)
  }
  
  # convert tss_annot to GRanges object
  tss_annot <- makeGRangesFromDataFrame(tss_annot, seqnames.field = "#chr", keep.extra.columns = TRUE)
  names(tss_annot) <- tss_annot$name
  
  # get tss annotations as 1bp coordinate
  tss_annot <- resize(tss_annot, width = 1, fix = "center")
  
  # GRanges object for E-G pairs in crispr using centers of candidate enhancers as coordinates
  eg_pairs <- makeGRangesFromDataFrame(crispr, seqnames.field = "chrom", start.field = "chromStart",
                                       end = "chromEnd", keep.extra.columns = TRUE)
  eg_pairs <- resize(eg_pairs, width = 1, fix = "center")
  
  # combine into a paired GRanges object containing enhancer and TSS/gene annotations for each pair
  eg_pairs <- Pairs(first = eg_pairs, second = tss_annot[crispr$measuredGeneSymbol])
  
  # compute distance to TSS and add to CRISPR E-G pairs
  crispr$dist_to_tss <- distance(eg_pairs)
  
  return(crispr)
  
}

# function to predict the direct hit rate for a given distance
predict_direct_hit_rate <- function(models, new_data, type) {
  
  # get model for specified type
  fit <- models[[type]]
  
  # predict direct hit rate for provided distances and exponentiate to get back to original scale
  predicted_log_direct_rate <- predict(fit, newdata = new_data)
  predicted_direct_rate <- exp(predicted_log_direct_rate)
  
  # cap predicted values at 1, because the rate can't be higher than that
  predicted_direct_rate <- pmin(predicted_direct_rate, 1)
  
  return(predicted_direct_rate)
  
}

# main function to compute probability of direct vs. indirect effect for E-G pairs in CRISPR data
predict_direct_vs_indirect <- function(crispr, type = c("negative", "positive", "all"), models,
                                       indirect_rates) {
  
  type <- match.arg(type)
  
  # predict the expected direct hit rate for each E-G pair based on distance to TSS
  crispr$direct_rate <- predict_direct_hit_rate(models, new_data = crispr,
                                                type = paste0("direct_rate_", type))
  
  # add corresponding indirect rate for the same type of hit
  indirect_col <- get_indirect_rate_column(type)
  crispr <- indirect_rates %>% 
    select(dataset_crispr_files, all_of(indirect_col)) %>% 
    left_join(crispr, ., by = c("Dataset" = "dataset_crispr_files"))
  
  # calculate direct vs. indirect effects probability
  crispr <- crispr %>% 
    mutate(direct_vs_indirect =  direct_rate / (direct_rate + indirect_rate))
  
  # add type to all columns added to crispr data for output
  rename_cols <- c("direct_rate", "indirect_rate", "direct_vs_indirect")
  names(rename_cols) <- paste(rename_cols, type, sep = "_")
  crispr <- dplyr::rename(crispr, all_of(rename_cols))
  
  return(crispr)
  
}

# helper function to get indirect rate rate column based on type
get_indirect_rate_column <- function(type) {
  if (type == "all") type <- "significant"
  indirect_col <- structure(paste0("positive_rate_", type), names = "indirect_rate")
  return(indirect_col)
}

## Predict direct vs. indirect hit probability -----------------------------------------------------

# load file containing fitted models
models <- readRDS(snakemake@input$models)

# load indirect rates for each dataset
indirect_rates <- read_tsv(snakemake@input$indirect_rates, show_col_types = FALSE)

# load crispr dataset id mapping file and add to indirect rates table
dataset_ids <- read_tsv(snakemake@input$dataset_ids, show_col_types = FALSE)
indirect_rates <- left_join(indirect_rates, dataset_ids, by = "dataset")

# load DC-TAP-seq processed CRISPR data and add 'Dataset' column
crispr <- read_tsv(snakemake@input$crispr, show_col_types = FALSE)
crispr <- mutate(crispr, Dataset = paste0(cell_type, "_DC_TAP"))

# compute distance to TSS for each CRISPR E-G pair based on the provided TSS universe if requested
if (is.null(snakemake@params$dist_col)) {
  tss <- read_tsv(snakemake@input$gene_tss, show_col_types = FALSE)
  crispr <- compute_dist_tss(crispr, tss_annot = tss)
} else {
  crispr <- mutate(crispr, dist_to_tss = !!sym(snakemake@params$dist_col))
}

# predict direct vs. indirect effects probability for each E-G pair based on distance to TSS
crispr <- predict_direct_vs_indirect(crispr, type = "negative", models = models,
                                     indirect_rates = indirect_rates)
crispr <- predict_direct_vs_indirect(crispr, type = "positive", models = models,
                                     indirect_rates = indirect_rates)
#crispr <- predict_direct_vs_indirect(crispr, type = "all", models = models,
#                                     indirect_rates = indirect_rates)

# remove temporary columns added to compute direct vs indirect effects probability
crispr <- select(crispr, -c(Dataset, dist_to_tss))

# save output to file
write_tsv(crispr, file = snakemake@output[[1]])
