## Annotate cis-results with element and pair annotations from ENCODE output

# save.image("pairs.rda")
# stop()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# load annotated results
cis_results <- readRDS(snakemake@input$cis_results)
encode_results <- fread(snakemake@input$encode_file)

# add valid connection and distance to TSS to cis-results
cis_results <- encode_results %>% 
  mutate(grna_target = paste0(chrom, ":", chromStart, "-", chromEnd),
         pair_uid = paste0(measuredEnsemblID, "|", grna_target)) %>% 
  select(grna_target, response_id = measuredEnsemblID, pair_uid, gene_symbol = measuredGeneSymbol,
         dist_to_tss = distToTSS, ValidConnection) %>% 
  distinct() %>% 
  left_join(cis_results, ., by = c("grna_target", "response_id"))

# get all invalid elements that overlap promoters or TSS
element_filter <- c("overlaps potential promoter", "TSS targeting guide(s)")
invalid_elements <- cis_results %>% 
  filter(ValidConnection %in% element_filter) %>%
  pull(grna_target) %>% 
  unique()

# get set of elements of valid pairs that are not part of the invalid elements list
valid_elements <- cis_results %>% 
  filter(!is.na(ValidConnection), !grna_target %in% invalid_elements) %>% 
  pull(grna_target) %>% 
  unique()

# create output with column labeling valid elements and make new simplified valid connection column
output <- cis_results %>% 
  mutate(valid_element = grna_target %in% valid_elements,
         valid_connection = ValidConnection == "TRUE") %>% 
  select(-ValidConnection)

# save output to file
fwrite(output, file = snakemake@output[[1]], sep = ",", quote = FALSE, na = "NA")
