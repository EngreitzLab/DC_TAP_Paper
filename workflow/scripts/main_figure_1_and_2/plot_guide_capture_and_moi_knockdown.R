# Fig. 1c: guide UMIs per guide per cell per 1,000 reads for CROP-seq, mU6 direct capture (pBA904) and
#          hU6 direct capture (sgOpti-CS) pilot experiments.
# Fig. 1d: knockdown at increasing MOI in K562 (GATA1 gRNAs), relative to the mean knockdown at MOI 1.
# Adapted from E. Jagoda's make_figure_1c_d_sfigs_v3.R (github.com/EngreitzLab/DC_TAP_paper_Evvie).
# usage: Rscript plot_guide_capture_and_moi_knockdown.R <input_dir> <fig1c.pdf> <fig1d.pdf>
suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]

# --- Fig. 1c ---
# mean reads per cell and total cells per pilot experiment (from the Cell Ranger summaries)
crop_reads_per_cell <- 1677.587
hu6_reads_per_cell <- 1091.170592310356772; hu6_total_cells <- 17322
mu6_reads_per_cell <- 1363; mu6_total_cells <- 25669

# CROP-seq: not every guide is detected per cell, so use UMIs of the top guide (max UMIs x fraction top guide)
crop <- read_tsv(file.path(input_dir, "P10_tap_all_guide_info.txt"), show_col_types = FALSE) %>%
  transmute(experiment = "Crop", umis_per_guide = replace_na(max_umis * fraction_top_guide, 0),
            reads_per_cell = crop_reads_per_cell)

# direct capture: per-cell protospacer calls; cells with no guide call count as 0 UMIs.
# Cells with several guides get the mean UMIs per guide.
read_protospacer_calls <- function(file, experiment, total_cells, reads_per_cell) {
  calls <- read_csv(file, show_col_types = FALSE)
  umis <- calls[[4]]  # 4th column: UMI count(s) per called guide, "|"-separated
  umis_per_guide <- map_dbl(str_split(as.character(umis), "\\|"), ~ sum(as.numeric(.x))) / calls$num_features
  tibble(experiment = experiment, umis_per_guide = c(rep(0, total_cells - nrow(calls)), umis_per_guide),
         reads_per_cell = reads_per_cell)
}
hu6 <- read_protospacer_calls(file.path(input_dir, "protospacer_calls_per_cell_052322_analysis_L-SG-TAP-Whole.csv"), "hU6", hu6_total_cells, hu6_reads_per_cell)
mu6 <- read_protospacer_calls(file.path(input_dir, "protospacer_calls_per_cell_052322_analysis_H-DC-TAP-Whole.csv"), "mU6", mu6_total_cells, mu6_reads_per_cell)

guide_capture <- bind_rows(crop, mu6, hu6) %>%
  mutate(umis_per_guide_per_1k_reads = umis_per_guide / reads_per_cell * 1000,
         experiment = factor(experiment, levels = c("Crop", "mU6", "hU6")))
print(guide_capture %>% group_by(experiment) %>% summarize(mean_umis_per_guide_per_1k_reads = mean(umis_per_guide_per_1k_reads), n = n()))

fig1c <- ggplot(guide_capture, aes(x = experiment, y = log10(umis_per_guide_per_1k_reads))) +
  geom_boxplot(fill = "grey80") +
  scale_x_discrete(labels = c(Crop = "CROP", mU6 = "pBA904\n(mU6)", hU6 = "sgOpti-CS\n(hU6)")) +
  labs(x = NULL, y = "Guide UMIs/guide/cell/1k reads (log10)") +
  theme_classic(base_size = 14)
ggsave(args[2], fig1c, width = 4, height = 4.5)

# --- Fig. 1d ---
moi_order <- c("moi1_combined_reps", "moi6_multi_moi", "moi9_multi_moi", "moi14_multi_moi", "moi20_hsg")
# this table was written with R row names, so read it with read.table
moi_knockdown <- read.table(file.path(input_dir, "input_file_for_k562_multi_moi1_panel.txt"), header = TRUE, sep = "\t") %>%
  filter(experiment_use %in% moi_order) %>%
  mutate(experiment_use = factor(experiment_use, levels = moi_order))
fig1d <- ggplot(moi_knockdown, aes(x = experiment_use, y = as.numeric(scaled_to_mean_exp_moi1))) +
  geom_boxplot(fill = "grey80") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_discrete(labels = c("1", "6", "9", "14", "20")) +
  coord_cartesian(ylim = c(0, 1.75)) +
  labs(x = "Multiplicity of infection (MOI)", y = "Relative KD (versus MOI 1)") +
  theme_classic(base_size = 14)
ggsave(args[3], fig1d, width = 4.5, height = 4.5)
