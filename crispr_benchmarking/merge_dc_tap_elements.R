## Resize and merge DC-TAP-seq elements for CRISPR benchmarking

library(tidyverse)
library(here)

# parameters
sig_col <- "significant_wo_pos_controls_20fdr"

# load DC-TAP-seq data
crispr_file <- here("WIP/resize_final_dc_tap_elements_20250820/TableS3_directEffectRate_chromatinCategories.tsv")
crispr <- read_tsv(crispr_file, show_col_types = FALSE)

# select required columns
epbench <- crispr %>% 
  select(chrom = resized_merged_targeting_chr_hg38, chromStart = resized_merged_targeting_start_hg38,
         chromEnd = resized_merged_targeting_end_hg38, name = resized_merged_element_gene_pair_identifier_hg38,
         EffectSize = pct_change_effect_size, chrTSS = chrTSS_hg38, startTSS = startTSS_hg38,
         endTSS = endTSS_hg38, measuredGeneSymbol = gene_symbol, Significant = all_of(sig_col),
         pValueAdjusted = sceptre_adj_p_value_wo_pos_controls,
         PowerAtEffectSize15 = power_at_effect_size_15_wo_pos_controls,
         PowerAtEffectSize25 = power_at_effect_size_25_wo_pos_controls,
         ValidConnection = Random_DistalElement_Gene, CellType = cell_type,
         ubiq_category, direct_rate_negative, indirect_rate_negative, direct_vs_indirect_negative,
         direct_rate_positive, indirect_rate_positive, direct_vs_indirect_positive,
         element_category) %>% 
  mutate(Regulated = Significant & EffectSize <= 0, .after = CellType)

# apply power filter
epbench <- epbench %>% 
  mutate(ValidConnection = case_when(
    ValidConnection == FALSE ~ "Not random E-G pair",
    Regulated == FALSE & ValidConnection == TRUE & PowerAtEffectSize25 < 0.8 ~ "underpowered",
    TRUE ~ "TRUE"
  ))

# only retain valid element gene pairs
epbench <- filter(epbench, ValidConnection == "TRUE")

# function to resolve cases where elements of a pair where merged
resolve_merge <- function(pairs) {
  
  if (nrow(pairs) > 1) {
    pairs <- arrange(pairs, desc(Regulated), desc(pValueAdjusted))
    pairs <- pairs %>% 
      slice_head(n = 1) %>% 
      mutate(EffectSize = NA_real_, pValueAdjusted = NA_real_,
             direct_rate_negative = mean(pairs$direct_rate_negative),
             direct_vs_indirect_negative = mean(pairs$direct_vs_indirect_negative),
             direct_rate_positive = mean(pairs$direct_rate_positive),
             direct_vs_indirect_positive = mean(pairs$direct_vs_indirect_positive),
             element_category = paste(unique(pairs$element_category), collapse = ","))
  }
  
  return(pairs)
  
}

# split dataset by cell type
epbench_k562 <- filter(epbench, CellType == "K562")
epbench_wtc11 <- filter(epbench, CellType == "WTC11")

# resolve merges
epbench_k562 <- epbench_k562 %>% 
  group_by(name) %>% 
  group_split() %>% 
  map_dfr(resolve_merge)

epbench_wtc11 <- epbench_wtc11 %>% 
  group_by(name) %>% 
  group_split() %>% 
  map_dfr(resolve_merge)

# combine data from cell types
output <- bind_rows(epbench_k562, epbench_wtc11)

# save to output file
output_file <- here("WIP/resize_final_dc_tap_elements_20250820/EPBenchmark_DC_TAPseq_randomElementGenePairs_0.8pwrAt25effect.tsv.gz")
write_tsv(output, file = output_file)
