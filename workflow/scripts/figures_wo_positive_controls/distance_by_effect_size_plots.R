# Script: distance_by_effect_size_plots.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/distance_by_effect_size_plots_FDR_wo_pos_controls.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggExtra)
  library(cowplot)
})

message("Loading input files")
screen_results <- read_tsv(snakemake@input$screen_results)

# Filter for Random_DistalElement_Gene
# significance and power from the final calls (FDR 20%, positive controls excluded)
results_with_element_gene_pair_categories_modified <- screen_results %>% 
  filter(Random_DistalElement_Gene) %>%
  mutate(significant = significant_wo_pos_controls_20fdr) %>%
  filter(significant | power_at_effect_size_15_wo_pos_controls_20fdr >= 0.8)


###  MAKE PLOTS ===============================================================

# Fig. 2b layout: Helvetica/Arial 7 pt tick labels, 8 pt titles, 3.5 pt points, marginal densities on the
# same y scale as the scatter (published range -40 to 50%)
pt <- 72.27 / 25.4
stroke_pt <- function(points) points / (pt * 0.75)  # ggplot linewidth for a stroke width in points
effect_colors <- c("Negative effect" = "red", "Positive effect" = "blue", "Not significant" = "grey")

plot_data <- results_with_element_gene_pair_categories_modified %>%
  mutate(percent_change = pct_change_effect_size,
         distance = abs(distance_to_gencode_gene_TSS),
         effect_color = case_when(significant & percent_change > 0 ~ "Positive effect",
                                  significant & percent_change < 0 ~ "Negative effect",
                                  TRUE ~ "Not significant"),
         effect_color = factor(effect_color, levels = names(effect_colors))) %>%
  arrange(significant)  # significant points drawn on top

plot_distance_by_es <- function(cell_type_to_plot, show_legend = FALSE) {
  p <- ggplot(filter(plot_data, cell_type == cell_type_to_plot),
              aes(x = distance / 1e6, y = percent_change, color = effect_color)) +
    geom_point(size = 1.64, stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = stroke_pt(0.95)) +
    # "\u00ad" draws a hyphen on the pdf device, where "-" is drawn as a minus sign
    scale_color_manual(values = effect_colors, name = "DE\u00adG pairs:") +
    scale_x_continuous(breaks = seq(0, 2, 0.5)) +
    scale_y_continuous(limits = c(-40, 50), breaks = c(-25, 0, 25, 50)) +
    labs(x = "Distance to TSS (Mb)", y = "Effect Size (%)") +
    theme_classic(base_size = 8) +
    theme(axis.title = element_text(size = 8, color = "black"),
          axis.title.y = element_text(margin = margin(r = 1)),
          axis.title.x = element_text(margin = margin(t = 4.5)),
          axis.text = element_text(size = 7, color = "black"),
          axis.text.y = element_text(margin = margin(r = 1.5)),
          axis.ticks.length = unit(2.5, "pt"),
          axis.line = element_line(linewidth = stroke_pt(0.95), color = "black"),
          axis.ticks = element_line(linewidth = stroke_pt(0.95), color = "black"),
          legend.position = if (show_legend) "top" else "none",
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 8),
          legend.key.size = unit(6, "pt"),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(8, 2, 2, 0))
  if (show_legend) return(p)
  pdf(NULL)  # ggMarginal opens a graphics device
  p <- ggMarginal(p, type = "density", margins = "both", groupColour = TRUE, groupFill = FALSE,
                  size = 7.6, linewidth = stroke_pt(0.75))
  dev.off()
  p
}

add_title <- function(p, title) {
  ggdraw(p) + draw_label(title, x = 0.53, y = 0.9, fontface = "bold", size = 8)
}

p1 <- plot_distance_by_es("K562")
p2 <- plot_distance_by_es("WTC11")
legend <- get_plot_component(plot_distance_by_es("K562", show_legend = TRUE) +
                               guides(color = guide_legend(override.aes = list(size = 1.5))), "guide-box-top")

# K562 and WTC11 side by side with a shared legend (Fig. 2b)
side_by_side <- plot_grid(legend, plot_grid(p1, p2, nrow = 1), ncol = 1, rel_heights = c(24.5, 179)) +
  draw_label("K562", x = 0.265, y = 0.87, fontface = "bold", size = 8) +
  draw_label("WTC11", x = 0.765, y = 0.87, fontface = "bold", size = 8)


### SAVE OUTPUT ===============================================================

message("Saving output files")
ggsave(plot = add_title(p1, "K562"), filename = snakemake@output$k562_distance_by_es,
       device = "pdf", width = 200 / 72, height = 180 / 72)
ggsave(plot = add_title(p2, "WTC11"), filename = snakemake@output$wtc11_distance_by_es,
       device = "pdf", width = 200 / 72, height = 180 / 72)
ggsave(plot = side_by_side, filename = snakemake@output$k562_wtc11_distance_by_es,
       device = "pdf", width = 400 / 72, height = 203.5 / 72)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)