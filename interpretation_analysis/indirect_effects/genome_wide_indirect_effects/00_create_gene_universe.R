## create gene universe for all genes in Gasperini et al., 2019

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(sceptre)
})

# load table with main screen results 
screen_results_file <- "/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/analyses/dc_tapseq_paper/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv"  
screen_results <- read_tsv(screen_results_file, show_col_types = FALSE)

# filter for K562 results only
k562_screen_results <- filter(screen_results, cell_type == "K562")

# load genome annotations
annot <- import("../../../results/genome_annotation_files/gencode.v26lift37.annotation.gtf.gz",
                format = "gtf")

# extract gene symbols, ids (without version) and chromosomes for all genes
genes <- annot[annot$type == "gene"] %>% 
  as.data.frame() %>% 
  select(gene_id, gene_name, gene_type, gene_chr = seqnames) %>% 
  mutate(gene_id = sub("\\..+", "", gene_id)) %>% 
  distinct()
  
# only retain annotations on regular chromosomes and from protein-coding or lincRNA genes
regular_chrs <- paste0("chr", c(1:22, "X", "M"))
genes_filt <- genes %>% 
  filter(gene_chr %in% regular_chrs) %>% 
  filter(gene_type %in% c("protein_coding", "lincRNA"))

# load Gasperini et al., 2019 sceptre object containing UMI counts used to compute TPM values
gasp_sceptre_file <- "/oak/stanford/groups/engreitz/Users/jgalante/DC_TAP_Paper/results/main_figure_1_and_2/duplicate_pairs_analysis/differential_expression/sceptre_diffex_input.rds"
gasp_sceptre <- readRDS(gasp_sceptre_file)  

# calculate TPM for each gene
umi_per_gene <- rowSums(gasp_sceptre@response_matrix[[1]])
tpm_per_gene <- umi_per_gene * 1e6 / sum(umi_per_gene)
tpm_per_gene <- enframe(tpm_per_gene, name = "gene_id", value = "tpm")

# add TPM to genes table
genes_filt <- left_join(genes_filt, tpm_per_gene, by = "gene_id")

# label genes that were part of the K562 DC-TAP-seq enhancer screen
k562_dc_tap_genes <- unique(k562_screen_results$gene_id)
dc_tap_genes_in_gasp_genes <- sum(k562_dc_tap_genes %in% genes_filt$gene_id)
message(dc_tap_genes_in_gasp_genes, " out of ", length(k562_dc_tap_genes),
        " K562 DC-TAP-seq genes found in Gasperini et al., 2019 data")
genes_filt <- mutate(genes_filt, k562_dc_tap = gene_id %in% k562_dc_tap_genes)

# save filtered and annotated gene list to output file
write_tsv(genes_filt, file = "gasperini_genes_annotated.tsv.gz")
  