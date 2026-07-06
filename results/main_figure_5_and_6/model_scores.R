suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggdist)
  library(cowplot)
})

format_number <- function(values) {
  vapply(values, function(value) {
    if (is.na(value)) return(NA_character_)
    
    if (value >= 100000) {
      # Scientific notation
      formatC(value, digits = 1, format = "E") %>% trimws()
    } else if (value %% 1 == 0) {
      # Integer < 100k
      format(value, big.mark = ",") %>% trimws()
    } else {
      # Float < 100k
      formatC(round(value, 1), format = "f", big.mark = ",", digits = 1)
    }
  }, character(1))
}

plot_score_by_category <- function(crispr, group_var, significance_var, e2g_threshold, e2g_score_name, out_prefix, scale_log10 = FALSE) {
	if (group_var == "element_category") {
		category_names <- c("H3K27me3 element", "CTCF element", "High H3K27ac", "H3K27ac", "No H3K27ac")
		cp <- c("#429130", "#49bcbc", "#c5373d", "#d9694a", "#c5cad7")
		names(cp) <- category_names
		label <- "Element category"
		order_levels <- c("High H3K27ac", "H3K27ac", "H3K27me3 element", "CTCF element", "No H3K27ac")
	} else if (group_var == "enhancerness") {
		cp <- c(`H3K27ac+ element` = "#D9694A", `Other element` = "#435369")
		label <- "Element type"
		order_levels <- rev(names(cp))
	} else if (group_var == "CTCF_category") {
		cp <- c(`CTCF element` = "#49bcbc", `Other element` = "#435369")
		label <- "Element type"
		order_levels <- rev(names(cp))
    } else if (group_var == "ubiq_category") {
		cp <- c("#792374", "#b778b3")
		names(cp) <- c("Ubiq. expr. gene", "Other gene")
		label <- "Promoter class"
		order_levels <- names(cp)
	} else if (group_var == "EffectSize_5pct") {
		cp <- c("#9b241c", "#006eae")
		names(cp) <- c("Over 5%", "Under 5%")
		label <- "Effect size"
		order_levels <- names(cp)
	} else if (group_var == "pDirect_category") {
		cp <- c("#0e3716", "#429130", "#96a0b3")
		names(cp) <- c("Over 90%", "50-90%", "Under 50%")
		label <- "P(Direct effect)"
		order_levels <- names(cp)
	} else if (group_var == "CellType") {
		cp <- c("#8f1426", "#c90072")
		names(cp) <- c("K562", "WTC11")
		label <- "Cell type"
		order_levels <- names(cp)
	} else {
		stop("Unsupported group_var")
	}

  ## all pairs
  res_all <- crispr %>%
    mutate(category = !!sym(group_var),
      positive = !!sym(significance_var),
      category = factor(category, levels = order_levels, ordered = TRUE))

  # adjust max scores for AG
  res_all <- res_all %>% 
    mutate(pred_value = ifelse(pred_value > 30, 30, pred_value))

    # just positives
	res <- res_all %>% 
		filter(positive == 1) %>% 
        group_by(category) %>% mutate(n_category = n()) %>% ungroup() %>% filter(n_category > 0)

	# summary labels
	smry <- res_all %>% 
		group_by(category) %>% 
		summarize(n_total = n(),
      n_positive = sum(positive == 1),
      n_true_positive = sum(pred_value >= e2g_threshold & positive == 1),
      n_false_positive = sum(pred_value >= e2g_threshold & positive == 0),
      n_true_negative = sum(pred_value < e2g_threshold & positive == 0),
      n_false_negative = sum(pred_value < e2g_threshold & positive == 1),
      mean_score_positives = mean(pred_value[positive == 1], na.rm = TRUE),
      mean_effect_size_positives = mean(EffectSize[positive == 1], na.rm = TRUE),

      effect_sizes_tp_str = paste(format(round(na.omit(EffectSize[pred_value >= e2g_threshold & positive == 1])), nsmall = 0),
                                   collapse = ", "),
      effect_sizes_fn_str = paste(format(round(na.omit(EffectSize[pred_value < e2g_threshold & positive == 1])), nsmall = 0),
                                   collapse = ", "),
      promoter_class_tp_str = paste(na.omit(ubiq_category[pred_value >= e2g_threshold & positive == 1]),
                                   collapse = ", "),
      promoter_class_fn_str = paste(na.omit(ubiq_category[pred_value < e2g_threshold & positive == 1]),
                                   collapse = ", "),
      .groups = "drop") %>% 
    filter(n_positive > 0) %>% 
    mutate(
      precision = n_true_positive / (n_true_positive + n_false_positive),
      recall = n_true_positive / n_positive,
      recall_pct = round(100 * recall, 1),
      recall_label = paste0(recall_pct, "%"),
		  n_label = paste0(n_positive))

	# Compute p-values if only 2 groups
	pval_df <- res %>%
		filter(n_distinct(category) == 2) %>%
		summarize(
			pval = tryCatch({
				wilcox.test(pred_value ~ category)$p.value
			}, error = function(e) NA_real_),
			.groups = "drop"
		)

	# Annotate p-values
	pval_df <- pval_df %>%
		mutate(pval_label = ifelse(!is.na(pval), paste0("p = ", signif(pval, 2)), NA),
			y = 1)

	# Plot
	pos_dodge <- 0.9
	pos_jitter <- position_jitterdodge(jitter.width = 0.2, dodge.width = pos_dodge)

    if (max(res$pred_value > 2)) {
        # alphagenome
        y_n_label <- 27
        y_recall <- 30
    } else if(max(res$pred_value < 1)) {
        # abc
        y_n_label <- 0.3
        y_recall <- 0.35
    } else {
        # encode-re2g
        y_n_label <- 0.95
        y_recall <- 1
    }
    
	p <- ggplot(res, aes(x = category, y = pred_value)) + 
		geom_boxplot(aes(fill = category), color = 'black', width = 0.7, outlier.shape = NA, position = position_dodge(pos_dodge)) +
		geom_jitter(aes(fill = category), color = 'black', size = 1.5, shape = 16, alpha = 0.7, position = pos_jitter) +
		geom_hline(yintercept = e2g_threshold, linewidth = 0.5, linetype = "dashed", color = "#435369") +
		geom_text(data = smry, aes(y = y_n_label, group = category, label = n_label), size = 4, color = "black", position = position_dodge(pos_dodge)) +
        geom_text(data = smry, aes(y = y_recall, group = category, label = recall_label), size = 4, color = "black", position = position_dodge(pos_dodge)) +
		#geom_text(data = pval_df, aes(x = Dataset, y = y, label = pval_label), inherit.aes = FALSE, size = 3.5, vjust = 0) +
		scale_fill_manual(values = cp) + scale_size_identity() +
		labs(x = "", y = e2g_score_name, color = label) + 
		theme_classic() +
		theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
			axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
			legend.position = "none")

    if (scale_log10) {
        p <- p + scale_y_log10()
    }

	w <- ifelse(length(order_levels) < 4, 3, 4.5)
  log_addn <- ifelse(scale_log10, ".log10", "")

  out_file_bp <- paste0(out_prefix, group_var, "_by_score_", significance_var, log_addn, ".pdf")
  out_file_pval <- paste0(out_prefix, group_var, "_by_score_", significance_var, ".pvalues.tsv")
  out_file_smry <- paste0(out_prefix, group_var, "_by_score_", significance_var, ".tsv")

  fwrite(pval_df, out_file_pval, sep = "\t", quote = FALSE, col.names = TRUE)
  fwrite(smry, out_file_smry, sep = "\t", quote = FALSE, col.names = TRUE)
  ggsave(out_file_bp, p, height = 5, width = w)
}


### MAIN
project_dir <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/2025_0227_CTCF_and_H3K27ac"
results_dir <- file.path(project_dir, "results", "2026_0626_dc_tap_revisions"); dir.create(results_dir, showWarnings = FALSE)

#  benchmarking version of elements: has ENCODE_rE2G, ABCdnase, AlphaGenome_RNAScorer
crispr_benchmark_file <- "/home/groups/engreitz/Users/emattei/git/CRISPR_comparison/results/DC_TAPseq_filtered/expt_pred_merged_annot.txt.gz"
models <- c("ENCODE_rE2G", "ABCdnase", "AlphaGenome_RNAScorer")
score_cols <- c("Score", "ABC.Score", "normalized_score")
model_score_names <- c("ENCODE-rE2G score", "ABC score", "AlphaGenome RNA gradient score")
thresholds <- c(0.201, 0.018, 11.837)

# --- plot ENCODE-rE2G score comparisons --- #
if (TRUE) {
  index <- 2
  model_id <- paste0(models[index], ".", score_cols[index])
  e2g_threshold <- thresholds[index]
  score_name <- model_score_names[index]
  model_name <- models[index]
  
  crispr_annot <- fread(crispr_benchmark_file) %>%
    mutate(EffectSize = case_when((name == "SAT2|chr17:7632382-7633443" & ExperimentCellType == "WTC11") ~ -4.733108597,
                                  (name == "EIF3K|chr19:38556711-38557710" & ExperimentCellType == "K562") ~ -3.420984314,
                                  (name == "FOXH1|chr8:144472292-144473578" & ExperimentCellType == "WTC11") ~ -12.03837989,
                                  TRUE ~ EffectSize)) %>% 
    filter(pred_uid == model_id) %>% 
      mutate(Dataset = "DC_TAP") %>%
      mutate(CellType = ExperimentCellType) %>% 
      mutate(EffectSize_5pct = ifelse(abs(EffectSize) >= 5, "Over 5%", "Under 5%")) %>% 
      mutate(pDirect_category = ifelse(direct_vs_indirect_negative > 0.9, "Over 90%",
        ifelse(direct_vs_indirect_negative > 0.5, "50-90%", "Under 50%")))
  
  crispr_annot_direct <- filter(crispr_annot, direct_vs_indirect_negative >= 0.5)
  crispr_annot_direct_enh <- filter(crispr_annot_direct, element_category %in% c("High H3K27ac", "H3K27ac"))

  out_dir <- file.path(results_dir, "e2g_model_scores"); dir.create(out_dir, showWarnings = FALSE)
  this_out_dir <- file.path(out_dir, model_name); dir.create(this_out_dir, showWarnings = FALSE)
  out_prefix <- paste0(this_out_dir, "/")
  out_prefix_direct <- paste0(this_out_dir, "/filter50pct_")
  out_prefix_direct_enh <- paste0(this_out_dir, "/filter50pct_enh_")

  plot_score_by_category(crispr_annot, "pDirect_category", "Regulated", e2g_threshold, score_name, out_prefix)
  plot_score_by_category(crispr_annot_direct, "EffectSize_5pct", "Regulated", e2g_threshold, score_name, out_prefix_direct)
  plot_score_by_category(crispr_annot_direct, "element_category", "Regulated", e2g_threshold, score_name, out_prefix_direct)
  plot_score_by_category(crispr_annot_direct, "ubiq_category", "Regulated", e2g_threshold, score_name, out_prefix_direct)
  plot_score_by_category(crispr_annot_direct_enh, "CellType", "Regulated", e2g_threshold, score_name, out_prefix_direct_enh)

}
