# Script: effect_size_plots.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/effect_size_plots_FDR_wo_pos_controls.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages(
  library(tidyverse)
)

message("Loading input files")
screen_results <- read_tsv(snakemake@input$screen_results)


### MAKE PLOT =================================================================

# Fig. 2c: significant negative-effect pairs, random distal elements (FDR 20%, positive controls excluded)
# and self-promoters (SCEPTRE significance). Adapted from J. Ray's effect_size_plots_RDE-SP_w_20_fdr.R.
effect_size_pairs <- screen_results %>%
  filter((Random_DistalElement_Gene & !is.na(sceptre_adj_p_value_wo_pos_controls) &
            sceptre_adj_p_value_wo_pos_controls < 0.2) |
           (selfPromoter & !is.na(significant) & significant)) %>%
  filter(pct_change_effect_size < 0) %>%
  mutate(Element_Type = factor(if_else(Random_DistalElement_Gene, "Random Distal Element", "Self Promoter"),
                               levels = c("Random Distal Element", "Self Promoter")))

effect_size_boxplot <- ggplot(effect_size_pairs, aes(x = cell_type, y = pct_change_effect_size)) +
  geom_boxplot(outlier.shape = NA, fill = "gray80", color = "black") +
  geom_jitter(shape = 21, size = 1.5, alpha = 0.8, stroke = 0.2, fill = "black", color = "black") +
  labs(x = "Cell Type", y = "Effect Size (%)") +
  coord_cartesian(ylim = c(-70, 0)) +
  facet_wrap(~ Element_Type) +
  theme_classic() +
  theme(legend.position = "none", strip.text = element_text(size = 12), axis.text.x = element_text(size = 12))

message("Pairs per element type and cell type:")
print(count(effect_size_pairs, Element_Type, cell_type))
message("Median effect size per element type:")
print(summarise(group_by(effect_size_pairs, Element_Type), median_effect_size = median(pct_change_effect_size), n = n()))

# All self-promoter pairs with a significant negative effect (positive controls have no FDR 20% call,
# so these use the original SCEPTRE significance)
positive_control_self_promoters_effect_size_boxplot <- screen_results %>%
  filter(selfPromoter, significant, pct_change_effect_size < 0) %>%
  ggplot(aes(x = cell_type, y = pct_change_effect_size)) +
  geom_boxplot(aes(fill = cell_type), outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(aes(color = significant), width = 0.2, size = 0.9, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "darkgrey", "TRUE" = "black")) +
  labs(
    x = "Cell Type",
    y = "Effect Size (%)"
  ) +
  theme_classic()


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(plot = effect_size_boxplot, filename = snakemake@output$effect_size_boxplot,
       device = "pdf", useDingbats = FALSE, height = 7, width = 7)
ggsave(plot = positive_control_self_promoters_effect_size_boxplot, filename = snakemake@output$positive_control_self_promoters_effect_size_boxplot,
       device = "pdf", height = 2.5, width = 3)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)