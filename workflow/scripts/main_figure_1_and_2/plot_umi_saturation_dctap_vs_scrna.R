## Figure 1e: mean UMIs per cell per gene vs sequencing depth, DC-TAP-seq vs scRNA-seq with guide capture.
## Adapted from J. Ray's scRNAvsdc-tap_reanalysis_2.R. The input tables are read-downsampling results from
##   DC-TAP-seq: ejagoda/multi_moi/10x_analysis/moi3_fixed
##   scRNA-seq:  jray/220210_CRISPRi_direct_capture/cellranger_outputs/P16-scRNA-seq-guide-capture
## Usage: Rscript plot_umi_saturation_dctap_vs_scrna.R <input_dir> <output_pdf>

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_pdf <- args[2]

dc_tap <- read.csv(file.path(input_dir, "complexity_df_dctap_moi3.csv"))
scrna <- read.csv(file.path(input_dir, "complexity_df_smaller_scrna_w_guide.csv"))

# Mean UMIs per cell per gene at each downsampled depth
average_over_genes <- function(df, method) {
  df %>%
    group_by(mean.total.reads.per.cell) %>%
    summarise(mean.umi.per.gene = mean(mean.umis.per.cell.per.gene), .groups = "drop") %>%
    mutate(method = method)
}
dc_tap_avg <- average_over_genes(dc_tap, "DC-TAP")
scrna_avg <- average_over_genes(scrna, "scRNA-seq")
combined_df <- bind_rows(dc_tap_avg, scrna_avg)

p <- ggplot(combined_df, aes(x = mean.total.reads.per.cell, y = mean.umi.per.gene, color = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(x = "Mean Total Reads per Cell", y = "Mean UMIs per Cell per Gene", color = "Method") +
  theme_classic(base_size = 14) +
  scale_color_manual(values = c("DC-TAP" = "red", "scRNA-seq" = "blue")) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "top"
  )
ggsave(output_pdf, p, width = 6, height = 5)

# UMI capture fold difference at the largest shared depth (depths rounded to the nearest 100 reads)
dc_tap_avg <- mutate(dc_tap_avg, rounded_reads = round(mean.total.reads.per.cell, -2))
scrna_avg <- mutate(scrna_avg, rounded_reads = round(mean.total.reads.per.cell, -2))
common_depths <- intersect(dc_tap_avg$rounded_reads, scrna_avg$rounded_reads)
if (length(common_depths) == 0) stop("No common sequencing depths found after rounding.")
max_shared_depth <- max(common_depths)

umi_at_depth <- function(df) {
  df %>% filter(rounded_reads == max_shared_depth) %>% pull(mean.umi.per.gene) %>% mean() %>% round(2)
}
dc_val <- umi_at_depth(dc_tap_avg)
scrna_val <- umi_at_depth(scrna_avg)
cat("Max shared depth (reads/cell):", max_shared_depth, "\n")
cat(" - DC-TAP UMI/cell/gene:", dc_val, "\n")
cat(" - scRNA-seq UMI/cell/gene:", scrna_val, "\n")
cat(" => Fold difference (DC-TAP / scRNA-seq):", round(dc_val / scrna_val, 2), "\n")
