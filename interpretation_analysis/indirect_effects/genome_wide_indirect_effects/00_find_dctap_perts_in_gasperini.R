## Intersect coordinates of K562 DC-TAP-seq candidate elements with Gasperini et al., 2019
## candidate elements to find elements perturbed in both screens

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(sceptre)
  library(GenomicRanges)
})

# load DC-TAP-seq cis-analysis results
screen_results_file <- "/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/analyses/dc_tapseq_paper/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv"  
screen_results <- fread(screen_results_file)

# filter for K562 results only
k562_screen_results <- filter(screen_results, cell_type == "K562")

# load Gasperini et al., 2019 sceptre object and extract grna targets table with candidate elements
gasp <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp = "sceptre_object_ondisc/sceptre_object.rds",
  response_odm_file_fp = "sceptre_object_ondisc/response.odm",
  grna_odm_file_fp = "sceptre_object_ondisc/grna.odm"
)
gasp_grna_targets <- gasp@grna_target_data_frame

# create GenomicRanges for DC-TAP-seq candidate elements
dctap_elements <- k562_screen_results %>% 
  select(intended_target_name_hg38, chr = targeting_chr_hg19, start = targeting_start_hg19,
         end = targeting_end_hg19) %>% 
  distinct() %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)

# create GenomicRanges for Gasperini candidate elements
gasp_elements <- gasp_grna_targets %>% 
  select(grna_target) %>% 
  distinct() %>% 
  separate(grna_target, into = c("chr", "start", "end"), sep = ":|-", remove = FALSE) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)

# intersect K562 DC-TAP-seq and Gasperini elements
ovl <- as.data.frame(findOverlapPairs(dctap_elements, gasp_elements))

# reformat table with DC TAP-seq elements overlapping Gasperini elements
output <- ovl %>% 
  select(dc_tapseq_element = `first.intended_target_name_hg38`,
         gasperini_element = `second.grna_target`, gasp_chr = second.X.seqnames,
         gasp_start = second.X.start, gasp_end = second.X.end)

# write to output file
fwrite(output, file = "overlapping_elements.tsv", sep = "\t")
