## Process massive trans analysis output from sceptre, and compute and analyze indirect effect rates

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
})

# load Gasperini et al., 2019 gene universe, including chromosome, TPM and DC_TAP-seq screen annotations
genes <- read_tsv("gasperini_genes_annotated.tsv.gz", show_col_types = FALSE)

# load elements overlapping between K562 DC-TAP-seq and Gasperini et al.
overlapping_elements <- read_tsv("overlapping_elements.tsv", show_col_types = FALSE)

# load cis-analysis results
cis_results <- read_tsv("/oak/stanford/groups/engreitz/Users/jgalante/DC_TAP_Paper/results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_element_gene_pair_categories.tsv",
                        show_col_types = FALSE)

# filter for distal element-gene pairs only (remove promoters)
cis_results <- filter(cis_results, DistalElement_Gene == TRUE)

# extract p-value cutoff for distal element-gene interactions
pval_threshold <- cis_results %>% 
  filter(DistalElement_Gene == TRUE, abs(distance_to_abc_canonical_TSS) <= 1e6) %>% 
  filter(significant == TRUE) %>% 
  slice_max(sceptre_p_value) %>% 
  pull(sceptre_p_value)

# Schema defining column formats for sceptre trans analysis results. Each parquet file carries its
# own dictionary, and Arrow cannot unify differing dictionaries across files during a
# compute/join. If column formats aren't set before any operation that materializes columns
# (e.g. the distinct() below), collect() or compute() fails with "Unifying differing dictionaries" 
sceptre_schema <- schema(
  response_id      = utf8(),
  grna_target      = utf8(),
  n_nonzero_trt    = int32(),
  n_nonzero_cntrl  = int32(),
  pass_qc          = bool(),
  p_value          = float64(),
  fold_change      = float64(),
  se_fold_change   = float64(),
  log_2_fold_change = float64()
)

# open sceptre results
trans_results <- open_dataset("sceptre_outputs/trans_results", schema = sceptre_schema)

# cast the dictionary-encoded id columns to plain strings. Each parquet file carries its
# own dictionary, and Arrow cannot unify differing dictionaries across files during a
# compute/join. This must happen before any operation that materializes these columns
# (e.g. the distinct() below), otherwise collect() fails with "Unifying differing dictionaries"
trans_results <- trans_results %>%
  mutate(grna_target = cast(grna_target, string()),
         response_id = cast(response_id, string()))

# only retain unique tested pairs that did pass sceptre QC
trans_results <- trans_results %>%
  filter(pass_qc == TRUE) %>%
  distinct()

# extract chromosome for each perturbation from grna_target id (format: "{chr}:{start}-{end}")
trans_results <- mutate(trans_results, pert_chr = sub("(chr.+):.+", "\\1", grna_target))

# add gene annotations to sceptre output
trans_results <- left_join(trans_results, genes, by = c("response_id" = "gene_id"))

# filter for trans-acting interactions and only retain pairs involving genes with TPM values
# (or DC-TAP-seq genes) or distal elements
distal_elements <- unique(pull(cis_results, intended_target_name_hg19))
trans_results <- trans_results %>%
  filter(pert_chr != gene_chr) %>%
  filter(!is.na(tpm) | k562_dc_tap == TRUE) %>% 
  filter(grna_target %in% distal_elements)

# label elements that overlap K562 DC-TAP-seq candidate elements
trans_results <- trans_results %>% 
  mutate(k562_dc_tap_gene = k562_dc_tap,
         k562_dc_tap_element = grna_target %in% overlapping_elements$gasperini_element)

# re-calculate significance based on nominal p-value cutoff from cis-analysis
trans_results <- mutate(trans_results, significant = p_value <= pval_threshold)

# perform all computations and materialize the filtered results into memory
trans_results <- collect(trans_results)

# compute indirect effects rate per gene
indirect_effects_genes <- trans_results %>% 
  group_by(response_id, gene_name, tpm, k562_dc_tap_gene) %>% 
  summarize(tested_trans_pairs = n(),
            significant_trans_pairs = sum(significant),
            trans_hit_rate = significant_trans_pairs / tested_trans_pairs,
            .groups = "drop")

# compute indirect effects rate per element
indirect_effects_elements <- trans_results %>% 
  group_by(grna_target, k562_dc_tap_element) %>% 
  summarize(tested_trans_pairs = n(),
            significant_trans_pairs = sum(significant),
            trans_hit_rate = significant_trans_pairs / tested_trans_pairs,
            .groups = "drop")

# save indirect effects tables to output file
write_tsv(indirect_effects_genes, file = "gasperini2019_indirect_effects_per_gene.tsv.gz")
write_tsv(indirect_effects_elements, file = "gasperini2019_indirect_effects_per_element.tsv.gz")
