# Fig. 2e: per-gRNA effect sizes measured by DC-TAP-seq vs qPCR (K562 and WTC11), with 95% CIs for both.
# Adapted from J. Ray's jr97a_qPCR_plot.R.
# usage: Rscript plot_dctap_vs_qpcr_effect_sizes.R <effect_size_table.tsv> <output.pdf>
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
effect_sizes <- read_tsv(args[1], show_col_types = FALSE, na = c("", "NA", "N/A")) %>%
  filter(gene_symbol != "SMIM1") %>%
  mutate(qpcr_xmin = pmin(qpcr_pctChange_95ci_lo, qpcr_pctChange_95ci_up),
         qpcr_xmax = pmax(qpcr_pctChange_95ci_lo, qpcr_pctChange_95ci_up))

pearson_r <- cor(effect_sizes$qpcr_pctChange, effect_sizes$dctap_pctChange, use = "pairwise.complete.obs")

dctap_vs_qpcr_plot <- ggplot(effect_sizes, aes(x = qpcr_pctChange, y = dctap_pctChange, color = target_type, shape = cell_type)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = qpcr_xmin, xmax = qpcr_xmax), width = 0.5, alpha = 0.9, orientation = "y") +
  geom_errorbar(aes(ymin = dctap_pctChange_95ci_lo, ymax = dctap_pctChange_95ci_up), width = 0.5, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = -95, y = -5, label = paste0("R = ", round(pearson_r, 2)), hjust = 0, size = 5) +
  labs(x = "qPCR Effect Size (%)", y = "DC-TAP-seq Effect Size (%)", color = "Target Type", shape = "Cell Type") +
  scale_color_grey(start = 0.2, end = 0.6, labels = c(Enhancer = "Distal Element", TSS = "Self Promoter")) +
  coord_cartesian(xlim = c(-100, 0), ylim = c(-100, 0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.line = element_line(color = "black"), axis.title = element_text(size = 18),
        axis.text = element_text(size = 16), legend.text = element_text(size = 14),
        legend.title = element_text(size = 14), legend.justification = c("left", "center"))

ggsave(args[2], plot = dctap_vs_qpcr_plot, width = 8, height = 6)
