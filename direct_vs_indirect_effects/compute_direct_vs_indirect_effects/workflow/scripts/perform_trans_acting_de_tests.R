## Run Sceptre for trans-acting effects

# save.image(paste0("trans_effects_", snakemake@wildcards$sample, ".rda"))
# stop()

# open log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# set seed for RNG
set.seed(snakemake@params$seed)

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(rtracklayer)
  library(sceptre)
})

# Prepare input data -------------------------------------------------------------------------------

message("\nLoading input data")

# load sceptre object containing all data, covariates and analysis parameters
sceptre_object <- readRDS(snakemake@input[[1]])

# load genome annotations and get chromosome of each gene
annot <- import(snakemake@input$genome)
genes <- annot %>% 
  as.data.frame() %>% 
  select(response_chr = seqnames, response_id = gene_id) %>% 
  mutate(response_id = sub("\\..+$", "", response_id)) %>%
  distinct()

# remove Y chromosome genes, since some of them are annotated to both X and Y
genes <- filter(genes, response_chr != "chrY")

# load list of cis E-G pairs and extract genes of valid cis interactions
cis_results <- fread(snakemake@input$cis_results)
cis_genes <- unique(pull(filter(cis_results, valid_connection == TRUE), response_id))

# Calculating average expression -------------------------------------------------------------------

# calculate average UMIs/gene and create table
avg_umi <- rowMeans(sceptre_object@response_matrix[[1]])
avg_umi <- data.table(gene = names(avg_umi), avg_expr = avg_umi)

# write average gene expression values to output file
fwrite(avg_umi, snakemake@output$gene_expr, quote = FALSE, na = "NA")

# Creating trans-acting E-G pairs to test  ---------------------------------------------------------

message("Creating trans pairs to test")

# create all perturbation-gene pairs
pairs <- construct_trans_pairs(sceptre_object)
         
# only retain pairs involving enhancers and genes on other chromosomes than the enhancer
trans_pairs <- pairs %>%
  filter(grepl("^chr.+", grna_target)) %>% 
  left_join(genes, by = "response_id") %>% 
  mutate(grna_target_chr = sub("^(chr.+):.+", "\\1", grna_target)) %>% 
  filter(grna_target_chr != response_chr)

# filter for genes in cis analysis only if specified
if (snakemake@params$include_genes == "cis") {
  message("Building trans-pairs using genes from cis analysis")
  trans_pairs <- filter(trans_pairs, response_id %in% cis_genes)
} else if (snakemake@params$include_genes == "all") {
  message("Building trans-pairs using all available genes")
} else {
  stop("Invalid 'include_genes' parameter", call. = FALSE)
}

# count how many genes are available per perturbation
available_genes_per_pert <- count(trans_pairs, grna_target, name = "available_genes")

# pick up to 100 trans genes per enhancer perturbation
message("Sampling ", snakemake@params$genes_per_pert, " trans genes per perturbation")
trans_pairs_sampled <- trans_pairs %>% 
  group_by(grna_target) %>% 
  sample_n(size = snakemake@params$genes_per_pert) %>% 
  ungroup()

# count how many genes are tested per perturbation and add to available gene counts
tested_genes_per_pert <- count(trans_pairs_sampled, grna_target, name = "tested_genes")
genes_per_pert <- left_join(available_genes_per_pert, tested_genes_per_pert, by = "grna_target")

# raise warning if for some perturbations less than 100 genes are available
if (any(genes_per_pert$tested_genes < 100)) {
  warning("Fewer than 100 genes available for testing for some perturbations", call. = FALSE)
}

# save available and tested genes per perturbation to output file
fwrite(genes_per_pert, file = snakemake@output$genes_per_pert, quote = FALSE, na = "NA")

# Performing differential expression tests ---------------------------------------------------------

message("\nRunning sceptre:\n")

# set analysis parameters
sceptre_object <- set_analysis_parameters(
  sceptre_object,
  discovery_pairs = distinct(select(trans_pairs_sampled, grna_target, response_id)),
  side = snakemake@params$side,
  grna_integration_strategy = snakemake@params$grna_integration_strategy
)

# assign guides to cells using specified approach
sceptre_object <- assign_grnas(
  sceptre_object,
  method = snakemake@params$guide_assignment_method,
  threshold = snakemake@params$guide_threshold
)

# run QC and filtering
sceptre_object <- run_qc(sceptre_object)

# run calibration check using negative control perturbations
guide_targets <- unique(sceptre_object@grna_target_data_frame$grna_target)
if ("non-targeting" %in% guide_targets) {
  sceptre_object <- run_calibration_check(sceptre_object, parallel = FALSE)
}

# run discovery analysis to test for trans-acting effects
sceptre_object <- run_discovery_analysis(sceptre_object, parallel = FALSE)

# write all results to output directory
message("\nSaving outputs to files")
invisible(write_outputs_to_directory(sceptre_object, dirname(snakemake@output$results)))
message("All done!")

# close log file connection
sink()
sink(type = "message")
close(log)
