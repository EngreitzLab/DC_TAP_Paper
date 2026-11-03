# Fig. 2i: qPCR effect of SLC2A3 distal-element guide #5742 vs a non-targeting guide (n = 4 replicates each).
# From J. Ray's SLC2A3_only.R; values come from jr97a_redos_k562_qpcr_validation_percentkd_qced_JRedit.xlsx
# (sheet "SLC2A3"). The element is labelled in hg38 coordinates, as in the manuscript.
# usage: Rscript plot_slc2a3_qpcr.R <output.pdf>
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
element_label <- "chr12:7960676-7960977 (guide 5742)"

replicate_values <- data.frame(
  Sample = rep(c(element_label, "Non-targeting"), each = 4),
  Value = c(-24.25684088, -38.52006691, -24.94929652, -37.97432223,   # guide 5742
            -13.55002211, 6.3770917, -1.32972055, 10.2048334))        # non-targeting
summary_values <- data.frame(
  Sample = c(element_label, "Non-targeting"),
  Mean = c(-31.42513164, 0.42554561),
  CI_low = c(-39.06, -9.995617519),
  CI_up = c(-23.6, 11.1057009))

# guide first, as in the manuscript
sample_order <- c(element_label, "Non-targeting")
replicate_values$Sample <- factor(replicate_values$Sample, levels = sample_order)
summary_values$Sample <- factor(summary_values$Sample, levels = sample_order)

slc2a3_qpcr_plot <- ggplot(summary_values, aes(x = Sample, y = Mean, fill = Sample)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width = 0.2, linewidth = 0.8) +
  geom_jitter(data = replicate_values, aes(x = Sample, y = Value), width = 0.1, size = 2.5,
              shape = 21, fill = "black", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(-50, 25) +
  labs(y = "Effect on SLC2A3 expression by qPCR (%)", x = "") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"))
ggsave(args[1], slc2a3_qpcr_plot, width = 4, height = 5)
