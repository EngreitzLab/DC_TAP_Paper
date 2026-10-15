# Fig. 2e: per-gRNA effect sizes measured by DC-TAP-seq vs qPCR (K562 and WTC11), with qPCR 95% CIs.
# Adapted from J. Huang's 250505_JR97b_wtc11_qpcr_validation_graphs_jredit.Rmd ("final check" version).
# usage: Rscript plot_dctap_vs_qpcr_effect_sizes.R <effect_size_table.tsv> <output.pdf>
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(ggrepel)
})
args <- commandArgs(trailingOnly = TRUE)
effect_sizes <- read_tsv(args[1], show_col_types = FALSE, na = c("", "NA", "N/A"))

dctap_vs_qpcr_plot <- ggplot(effect_sizes, aes(x = qpcr_pctChange, y = dctap_pctChange, color = target_type, shape = cell_type)) +
  geom_point(size = 6) +
  geom_errorbarh(aes(xmin = qpcr_pctChange_95ci_lo, xmax = qpcr_pctChange_95ci_up), linewidth = 0.8, height = 2, alpha = 0.6) +
  geom_text_repel(aes(label = gene_symbol), size = 5, max.overlaps = 20, force = 2, force_pull = 0.1,
                  min.segment.length = 0, point.padding = 20, show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  scale_color_manual(values = c("Enhancer" = "#F8766D", "TSS" = "#00BFC4")) +
  scale_shape_manual(values = c("K562" = 16, "WTC11" = 17)) +
  coord_cartesian(xlim = c(-100, 0), ylim = c(-100, 0)) +
  stat_cor(method = "pearson", aes(label = after_stat(r.label), group = 1), color = "black", size = 8, show.legend = FALSE) +
  labs(title = "DC-TAP-seq vs qPCR Effect Sizes", x = "qPCR Effect Size (%)", y = "DC-TAP Effect Size (%)",
       color = "Target Type", shape = "Cell Type") +
  theme_classic(base_size = 28) +
  theme(legend.position = "right", plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20), axis.text = element_text(size = 20),
        legend.text = element_text(size = 20), legend.title = element_text(size = 20))

ggsave(args[2], plot = dctap_vs_qpcr_plot, width = 16, height = 12)
