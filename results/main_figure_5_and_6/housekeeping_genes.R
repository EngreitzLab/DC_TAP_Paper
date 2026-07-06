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

plot_percent_positive_enrichment_effect_size_combined <- function(crispr, group_var, direct_effect_weighted, direct_effect_threshold, all_power, out_prefix) {
  sig_vars <- c("Significant", "Regulated")

  if (direct_effect_weighted & !is.null(direct_effect_threshold)) {stop("Choose at most one of direct effect weighting or filtering!")}

  plots <- list()
  combined_pct_pos <- list()
  combined_enr <- list()
  power_enr <- list()
  combined_es <- list()
  max_prop <- list()
  max_prop_pairs <- list()
  max_es <- list()

  ds_cp <- c(K562_DC_TAP = "#006eae", WTC11_DC_TAP = "#00488d", DC_TAP = "#005a9d", Combined = "#435369",
    Gasperini2019 = "#d3a9ce", Nasser2021 = "#b778b3", Schraivogel2020 = "#a64791",
    Klann = "#d3a9ce", Morris = "#d3a9ce", Xie = "#d3a9ce") 

  for (significance_var in sig_vars) {
    z <- qnorm(0.05/2, lower.tail=FALSE)
    pos_dodge <- 0.9

    summarize_pairs_vars <- c("elementName", "measuredGeneSymbol")
    summarize_pairs_label <- "E-G pairs"

    if (group_var == "element_category") {
      category_names <- c("H3K27me3 element", "CTCF element", "High H3K27ac", "H3K27ac", "No H3K27ac")
      cp <- c("#429130", "#49bcbc", "#c5373d", "#d9694a", "#c5cad7")
      names(cp) <- category_names
      label <- "Element category"
      order_levels <- rev(names(cp))
      summarize_vars <- c("elementName")
      summarize_label <- "elements"

    } else if (group_var == "enhancerness") {
      cp <- c(`H3K27ac+ element` = "#d9694a", `Other element` = "#435369")
      label <- "Element type"
      order_levels <- rev(names(cp))
      summarize_vars <- c("elementName")
      summarize_label <- "elements" 

    } else if (group_var == "ubiq_category") {
      cp <- c("#792374", "#b778b3")
      names(cp) <- c("Ubiq. expr. gene", "Other gene")
      label <- "Promoter class"
      order_levels <- names(cp)
      summarize_vars <- c("measuredGeneSymbol")
      summarize_label <- "genes"

    } else if (group_var == "distance_category") {
      cp <- c("#002359", "#00488d", "#006eae", "#5496ce", "#9bcae9")
      names(cp) <- c("0-10 kb", "10-100 kb", "100-250 kb", "250 kb-1 Mb", "1 Mb-2 Mb")
      label <- "Distance to TSS"
      order_levels <- names(cp)
      summarize_vars <- summarize_pairs_vars
      summarize_label <- summarize_pairs_label

    } else {
      stop("Unsupported group_var")
    }

    ## prepare data
    if (direct_effect_weighted) {
        crispr_this <- crispr %>%
            mutate(positive_indicator = ifelse(!!sym(significance_var) == TRUE, direct_vs_indirect, 0))
        if (all_power) {
            crispr_powered <- crispr_this
        } else {
            crispr_powered <- crispr_this %>% filter(WellPowered == TRUE)
        }
        

    } else { # filter by threshold (not significant OR significant with high direct rate
        crispr_this <- crispr %>%
            filter(!(!!sym(significance_var) == TRUE & direct_vs_indirect < direct_effect_threshold)) %>%
            mutate(positive_indicator = ifelse(!!sym(significance_var) == TRUE, 1, 0))
        crispr_powered <- crispr_this #%>% filter(WellPowered == TRUE)
    }

    # hit rate
    pct_pos <- crispr_powered %>% 
      group_by(Dataset, pick(all_of(summarize_vars))) %>%
      #mutate(anyRegulated = any(!!sym(significance_var))) %>%
      mutate(anyRegulated = max(positive_indicator)) %>% 
      ungroup() %>%
      select(all_of(summarize_vars), Dataset, anyRegulated, category = !!sym(group_var)) %>% 
      distinct() %>% 
      group_by(Dataset, category) %>% 
      summarize(n_tested_category = n(), n_positive_category = sum(anyRegulated), .groups = "drop") %>%
      mutate(prop_positive_elements = n_positive_category / n_tested_category,
             prop_adjust = (n_positive_category + 2) / (n_tested_category + 4),
             SE_prop = sqrt(prop_adjust * (1 - prop_adjust) / (n_tested_category + 4)),
             CI_prop_low = pmax(prop_positive_elements - z * SE_prop, 0),
             CI_prop_high = prop_positive_elements + z * SE_prop,
             n_label = paste0("(", format_number(n_positive_category), ")"),
             category = factor(category, levels = order_levels, ordered = TRUE),
             Dataset = factor(Dataset, levels = names(ds_cp), ordered = TRUE),
             significance = significance_var,
             grouped_by = summarize_label) %>%
      filter(!is.na(category))
    pct_pos$y_label_pos <- max(pct_pos$CI_prop_high, na.rm = TRUE) * 1.1

    # hit rate in terms of pairs
    pct_pos_pairs <- crispr_powered %>% 
      group_by(Dataset, pick(all_of(summarize_pairs_vars))) %>%
      #mutate(anyRegulated = any(!!sym(significance_var))) %>%
      mutate(anyRegulated = max(positive_indicator)) %>% 
      ungroup() %>%
      select(all_of(summarize_pairs_vars), Dataset, anyRegulated, category = !!sym(group_var)) %>% 
      distinct() %>% 
      group_by(Dataset, category) %>% 
      summarize(n_tested_category = n(), n_positive_category = sum(anyRegulated), .groups = "drop") %>%
      mutate(prop_positive_elements = n_positive_category / n_tested_category,
             prop_adjust = (n_positive_category + 2) / (n_tested_category + 4),
             SE_prop = sqrt(prop_adjust * (1 - prop_adjust) / (n_tested_category + 4)),
             CI_prop_low = pmax(prop_positive_elements - z * SE_prop, 0),
             CI_prop_high = prop_positive_elements + z * SE_prop,
             n_label = paste0("(", format_number(n_positive_category), ")"),
             category = factor(category, levels = order_levels, ordered = TRUE),
             Dataset = factor(Dataset, levels = names(ds_cp), ordered = TRUE),
             significance = significance_var,
             grouped_by = summarize_pairs_label) %>%
      filter(!is.na(category))
    pct_pos_pairs$y_label_pos <- max(pct_pos_pairs$CI_prop_high, na.rm = TRUE) * 1.1

    # enrichment
    enr <- crispr_powered %>% 
      group_by(Dataset) %>%
      #mutate(n_tested = n(), n_positive = sum(!!sym(significance_var))) %>%
      mutate(n_tested = n(), n_positive = sum(positive_indicator)) %>% 
      group_by(Dataset, category = !!sym(group_var), n_tested, n_positive) %>%
      summarize(n_tested_category = n(),
        n_positive_category = sum(positive_indicator), # sum(!!sym(significance_var)),
        .groups = "drop") %>% 
      mutate(prop_tested = n_tested_category/n_tested,
             prop_pos = ifelse(n_positive > 0, n_positive_category/n_positive, 0),
             enrichment = ifelse(prop_tested > 0, prop_pos/prop_tested, NA),
             enrichment_label = round(enrichment, 2),
             SE_log_enr = ifelse(n_positive_category > 0 & n_tested_category > 0,
                                 sqrt(((n_positive - n_positive_category) / n_positive_category) / n_positive) +
                                 ((n_tested - n_tested_category) / n_tested_category) / n_tested,
                                 NA),
             CI_enr_low = exp(log(enrichment) - z * SE_log_enr),
             CI_enr_high = exp(log(enrichment) + z * SE_log_enr),
             p_enr = phyper(n_positive_category, n_positive, n_tested, n_positive_category + n_tested_category, log.p = FALSE, lower.tail = FALSE),
             p_adjust_enr = p.adjust(p_enr, method = "bonferroni"),
             sign_label = case_when(p_adjust_enr < 0.001 ~ "***", p_adjust_enr < 0.01 ~ "**", p_adjust_enr < 0.05 ~ "*", TRUE ~ ""),
             category = factor(category, levels = order_levels, ordered = TRUE),
             Dataset = factor(Dataset, levels = names(ds_cp), ordered = TRUE),
             significance = significance_var) %>%
      filter(!is.na(category))

    # effect size
    es_col <- ifelse(significance_var == "Regulated", "EffectSize", "absEffectSize")
    es <- crispr_powered %>% 
        filter(!!sym(significance_var) == 1) %>% 
        mutate(category = !!sym(group_var),
            absEffectSize = abs(EffectSize),
            category = factor(category, levels = order_levels, ordered = TRUE),
            Dataset = factor(Dataset, levels = names(ds_cp), ordered = TRUE),
            plotEffectSize = !!sym(es_col))

    es_smry <- es %>% 
        group_by(Dataset, category) %>% 
        summarize(n_pairs = n(),
            n_label = paste0("(", n_pairs, ")"),
            mean_EffectSize = mean(EffectSize, na.rm = TRUE), 
            mean_absEffectSize = mean(absEffectSize, na.rm = TRUE),
            sd_absEffectSize = sd(absEffectSize, na.rm = TRUE)) %>% ungroup() %>%
        mutate(se_absEffectSize = sd_absEffectSize / sqrt(n_pairs),
            CI_ES_low = mean_absEffectSize - z * se_absEffectSize,
            CI_ES_high = mean_absEffectSize + z * se_absEffectSize,
            significance = significance_var)
    es_smry$y_label_pos <- min(100, max(es$plotEffectSize, na.rm = TRUE) + 3)

    ## plot
    # params
    max_prop[[significance_var]] <- max(pct_pos$CI_prop_high, na.rm = TRUE)
    max_prop_pairs[[significance_var]] <- max(pct_pos_pairs$CI_prop_high, na.rm = TRUE)

    y_label1 <- ifelse(significance_var == "Regulated",
                       paste0("Fraction tested ", summarize_label, "\nwith 1+ downregulated hit (# hits)"),
                       paste0("Fraction tested ", summarize_label, "\nwith 1+ up/downregulated hit (# hits)"))
    
    y_label4 <- ifelse(significance_var == "Regulated",
                       paste0("Fraction tested ", summarize_pairs_label, "\nwith 1+ downregulated hit (# hits)"),
                       paste0("Fraction tested ", summarize_pairs_label, "\nwith 1+ up/downregulated hit (# hits)"))

    y_label2 <- ifelse(significance_var == "Regulated",
                       "Enrichment of downregulated hits\n(% significant / % tested)",
                       "Enrichment of up/downregulated hits\n(% significant / % tested)")
                    
    max_es[[significance_var]] <- min(100, max(es$plotEffectSize, na.rm = TRUE))
    y_label3 <- ifelse(significance_var == "Regulated",
                       "% effect size of downregulated hits (# hits)",
                       "Abs (% effect size) of up/downregulated hits (# hits)")
    # hit rate
    pct_pos_plot <- filter(pct_pos, prop_positive_elements > 0)
    p1 <- ggplot(pct_pos_plot, aes(x = Dataset, y = prop_positive_elements, fill = category)) + 
      geom_col(position = position_dodge(width = pos_dodge)) +
      geom_linerange(aes(ymin = CI_prop_low, ymax = CI_prop_high), position = position_dodge(width = pos_dodge), linewidth = 0.5, color = "black") +
      geom_text(aes(y = y_label_pos, label = n_label), size = 2, position = position_dodge(width = pos_dodge)) +
      ylim(c(0, max_prop[[significance_var]] * 1.15)) +
      scale_fill_manual(values = cp) +
      labs(x = "Dataset", y = y_label1, fill = label) +
      theme_classic() +
      theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
        axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
        legend.position = "none")
    
    # hit rate, pairs
    pct_pos_pairs_plot <- filter(pct_pos_pairs, prop_positive_elements > 0)
    p4 <- ggplot(pct_pos_pairs_plot, aes(x = Dataset, y = prop_positive_elements, fill = category)) + 
      geom_col(position = position_dodge(width = pos_dodge)) +
      geom_linerange(aes(ymin = CI_prop_low, ymax = CI_prop_high), position = position_dodge(width = pos_dodge), linewidth = 0.5, color = "black") +
      geom_text(aes(y = y_label_pos, label = n_label), size = 2, position = position_dodge(width = pos_dodge)) +
      ylim(c(0, max_prop_pairs[[significance_var]] * 1.15)) +
      scale_fill_manual(values = cp) +
      labs(x = "Dataset", y = y_label4, fill = label) +
      theme_classic() +
      theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
        axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
        legend.position = "none")

    # enrichment
    enr_plot <- enr %>% filter(n_positive_category > 0, CI_enr_high < 100)
    p2 <- ggplot(enr_plot, aes(x = Dataset, y = enrichment, fill = category)) + 
      geom_hline(yintercept = 1, linetype = "dashed", color = "#c5cad7") +
      geom_col(position = position_dodge(width = pos_dodge)) +
      geom_linerange(aes(ymin = CI_enr_low, ymax = CI_enr_high), position = position_dodge(width = pos_dodge), linewidth = 0.5, color = "black") +
      geom_text(aes(y = CI_enr_high + 0.2, label = sign_label), size = 2, position = position_dodge(width = pos_dodge)) +
      scale_fill_manual(values = cp) +
      labs(x = "Dataset", y = y_label2, fill = label) +
      theme_classic() +
      theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
        axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
        legend.position = "none")
    
    # effect size
    p3 <- ggplot(es, aes(x = Dataset, y = plotEffectSize)) + 
        geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", color = "#c5cad7") +
        geom_boxplot(aes(color = category), fill = NA, width = 0.7, outlier.shape = 16, outlier.size = 0.75, position = position_dodge(pos_dodge)) +
        geom_text(data = es_smry, aes(y = y_label_pos, group = category, label = n_label), size = 2, color = "black", position = position_dodge(pos_dodge)) +
        ylim(c(NA, max_es[[significance_var]]) + 6) +
		scale_color_manual(values = cp) + 
		labs(x = "Dataset", y = y_label3) + 
      theme_classic() +
      theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
        axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
        legend.position = "none")

    plots[[significance_var]] <- list(p1, p2, p3, p4)
    combined_pct_pos[[significance_var]] <- rbind(pct_pos, pct_pos_pairs) %>% select(-n_label, -prop_adjust)
    combined_enr[[significance_var]] <- enr %>% select(-sign_label)
    combined_es[[significance_var]] <- es_smry %>% select(-n_label)
  }

    ## power enrichment
    enr_power <- crispr %>% 
        group_by(Dataset) %>%
        mutate(n_tested = n(), n_powered = sum(WellPowered)) %>% 
        group_by(Dataset, category = !!sym(group_var), n_tested, n_powered) %>%
        summarize(n_tested_category = n(), n_powered_category = sum(WellPowered), .groups = "drop") %>% 
        mutate(prop_tested = n_tested_category/n_tested,
            prop_pos = ifelse(n_powered > 0, n_powered_category / n_powered, 0),
            enrichment = ifelse(prop_tested > 0, prop_pos/prop_tested, NA),
            enrichment_label = round(enrichment, 2),
            SE_log_enr = ifelse(n_powered_category > 0 & n_tested_category > 0,
                                sqrt(((n_powered - n_powered_category) / n_powered_category) / n_powered) +
                                ((n_tested - n_tested_category) / n_tested_category) / n_tested,
                                NA),
            CI_enr_low = exp(log(enrichment) - z * SE_log_enr),
            CI_enr_high = exp(log(enrichment) + z * SE_log_enr),
            p_enr = phyper(n_powered_category, n_powered, n_tested, n_powered_category + n_tested_category, log.p = FALSE, lower.tail = FALSE),
            p_adjust_enr = p.adjust(p_enr, method = "bonferroni"),
            sign_label = case_when(p_adjust_enr < 0.001 ~ "***", p_adjust_enr < 0.01 ~ "**", p_adjust_enr < 0.05 ~ "*", TRUE ~ ""),
            category = factor(category, levels = order_levels, ordered = TRUE),
            Dataset = factor(Dataset, levels = names(ds_cp), ordered = TRUE),
            significance = significance_var) %>%
    filter(!is.na(category))

    # make plt for power enrichment
    p5 <- ggplot(enr_power, aes(x = Dataset, y = enrichment, fill = category)) + 
      geom_hline(yintercept = 1, linetype = "dashed", color = "#c5cad7") +
      geom_col(position = position_dodge(width = pos_dodge)) +
      geom_linerange(aes(ymin = CI_enr_low, ymax = CI_enr_high), position = position_dodge(width = pos_dodge), linewidth = 0.5, color = "black") +
      geom_text(aes(y = CI_enr_high + 0.2, label = sign_label), size = 2, position = position_dodge(width = pos_dodge)) +
      scale_fill_manual(values = cp) +
      labs(x = "Dataset", y = "Enrichment\n(% well-powered / % tested)", fill = label) +
      theme_classic() +
      theme(axis.text = element_text(size = 8, color = "#000000"), axis.title = element_text(size = 9),
        axis.ticks = element_line(color = "#000000"), axis.ticks.x = element_blank(),
        legend.position = "top")

    binary_groups <- c("ubiq_category", "enhancerness")
    n_dataset <- length(unique(crispr$Dataset))
    w <- ifelse(group_var %in% binary_groups, 12 + (n_dataset - 2) * 3, 16 + (n_dataset - 2) * 4)

    ## save plots
    ggsave(paste0(out_prefix, group_var, "_power_enrichment.pdf"), p5, height = 4, w = w/4)

    # plot combined others
  legend <- cowplot::get_plot_component(plots[["Significant"]][[1]] + theme(legend.position = "top"), 'guide-box-top')

  combined <- plot_grid(
    legend,
    plot_grid(plots[["Significant"]][[3]], plots[["Significant"]][[2]], plots[["Significant"]][[1]], plots[["Significant"]][[4]], nrow = 1, rel_widths = c(1, 1, 1, 1)),
    plot_grid(plots[["Regulated"]][[3]], plots[["Regulated"]][[2]],  plots[["Regulated"]][[1]], plots[["Regulated"]][[4]], nrow = 1, rel_widths = c(1, 1, 1, 1, 1)),
    ncol = 1, rel_heights = c(0.15, 1, 1))

  ggsave(paste0(out_prefix, group_var, "_combined_summary.pdf"), combined, height = 8, width = w)

  fwrite(enr_power, paste0(out_prefix, group_var, "_power_enrichment.tsv"), sep = "\t")
  fwrite(bind_rows(combined_pct_pos), paste0(out_prefix, group_var, "_crispr_results_combined_percent_positives.tsv"), sep = "\t")
  fwrite(bind_rows(combined_enr), paste0(out_prefix, group_var, "_crispr_results_combined_enrichment.tsv"), sep = "\t")
  fwrite(bind_rows(combined_es), paste0(out_prefix, group_var, "_crispr_results_combined_effect_sizes.tsv"), sep = "\t")
}

compare_metrics_across_groups <- function(crispr, group_var, out_prefix, pairs = c("Category", "Dataset"),
    direct_effect_weighted, direct_effect_threshold, all_power) {
  if (direct_effect_weighted & !is.null(direct_effect_threshold)) {stop("Choose either direct effect weighting or filtering!")}

  pairs <- match.arg(pairs)
  sig_vars <- c("Significant", "Regulated")
  summarize_pairs_vars <- c("elementName", "measuredGeneSymbol")
  summarize_pairs_label <- "pairs"

  if (group_var == "element_category") {
    category_names <- c("H3K27me3 element", "CTCF element", "High H3K27ac", "H3K27ac", "No H3K27ac")
    cp <- c("#429130", "#49bcbc", "#c5373d", "#d9694a", "#c5cad7")
    names(cp) <- category_names
    label <- "Element category"
    order_levels <- rev(names(cp))
    summarize_vars <- c("elementName")
    summarize_label <- "elements"
  } else if (group_var == "enhancerness") {
    cp <- c(`H3K27ac+ element` = "#D9694A", `Other element` = "#435369")
    label <- "Element type"
    order_levels <- rev(names(cp))
    summarize_vars <- c("elementName")
    summarize_label <- "elements"
  } else if (group_var == "ubiq_category") {
    cp <- c("#792374", "#b778b3")
    names(cp) <- c("Ubiq. expr. gene", "Other gene")
    label <- "Promoter class"
    order_levels <- names(cp)
    summarize_vars <- c("measuredGeneSymbol")
    summarize_label <- "genes"
  } else if (group_var == "distance_category") {
    cp <- c("#002359", "#00488d", "#006eae", "#5496ce", "#9bcae9")
    names(cp) <- c("0-10 kb", "10-100 kb", "100-250 kb", "250 kb-1 Mb", "1 Mb-2 Mb")
    label <- "Distance to TSS"
    order_levels <- names(cp)
    summarize_vars <- summarize_pairs_vars
    summarize_label <- summarize_pairs_label
  } else {
    stop("Unsupported group_var")
  }

  all_results <- list()

  for (significance_var in sig_vars) {

    ## prepare data
    if (direct_effect_weighted) {
        crispr_this <- crispr %>%
            mutate(positive_indicator = ifelse(!!sym(significance_var) == TRUE, direct_vs_indirect, 0))

    } else { # filter by threshold (not significant OR significant with high direct rate
        crispr_this <- crispr %>%
            filter(!(!!sym(significance_var) == TRUE & direct_vs_indirect < direct_effect_threshold)) %>%
            mutate(positive_indicator = ifelse(!!sym(significance_var) == TRUE, 1, 0))
        
    }

    if (all_power) {
        crispr_powered <- crispr_this
    } else {
        crispr_powered <- crispr_this %>% filter(WellPowered == TRUE)
    }

    crispr_sub <- crispr_powered %>%
      mutate(category = !!sym(group_var),
             absEffectSize = abs(EffectSize)) %>%
      filter(!is.na(category))

    hit_data <- crispr_sub %>%
      group_by(Dataset, pick(all_of(summarize_vars))) %>%
      #mutate(anyRegulated = any(!!sym(significance_var))) %>%
      mutate(anyRegulated = max(positive_indicator)) %>%
      ungroup() %>%
      select(all_of(summarize_vars), Dataset, anyRegulated, category) %>%
      distinct() %>%
      group_by(Dataset, category) %>%
      summarize(n_total = n(), n_hit = sum(anyRegulated), hit_rate = n_hit / n_total, .groups = "drop")
    
    hit_data_pairs <- crispr_sub %>%
      group_by(Dataset, pick(all_of(summarize_pairs_vars))) %>%
      mutate(anyRegulated = max(positive_indicator)) %>%
      ungroup() %>%
      select(all_of(summarize_vars), Dataset, anyRegulated, category) %>%
      distinct() %>%
      group_by(Dataset, category) %>%
      summarize(n_total = n(), n_hit = sum(anyRegulated), hit_rate = n_hit / n_total, .groups = "drop")

    hit_data_list <- list(smry = hit_data, pairs = hit_data_pairs)
    hit_pairs_list <- list()
    for (hd in names(hit_data_list)) {
        if (pairs == "Category") {
            hit_pairs_list[[hd]] <- hit_data_list[[hd]] %>%
                group_by(Dataset) %>%
                filter(n_distinct(category) >= 2) %>%
                group_modify(~{
                combs <- combn(unique(.x$category), 2, simplify = FALSE)
                do.call(rbind, lapply(combs, function(pair) {
                    d1 <- .x[.x$category == pair[1], ]
                    d2 <- .x[.x$category == pair[2], ]
                    if (nrow(d1) == 0 | nrow(d2) == 0) return(NULL)
                    x <- c(d1$n_hit, d2$n_hit)
                    n <- c(d1$n_total, d2$n_total)
                    p <- tryCatch(prop.test(x, n)$p.value, error = function(e) NA_real_)
                    data.frame(group1 = pair[1], group2 = pair[2],
                            group1_value = d1$hit_rate, group2_value = d2$hit_rate,
                            p_value = p, test = "prop.test")
                }))
                }, .groups = "drop") %>%
                mutate(metric = ifelse(hd == "smry", paste0("HitRate_", summarize_label),
                    paste0("HitRate_", summarize_pairs_label)))
        } else {
            hit_pairs_list[[hd]] <- hit_data_list[[hd]] %>%
                group_by(category) %>%
                filter(n_distinct(Dataset) >= 2) %>%
                group_modify(~{
                combs <- combn(unique(.x$Dataset), 2, simplify = FALSE)
                do.call(rbind, lapply(combs, function(pair) {
                    d1 <- .x[.x$Dataset == pair[1], ]
                    d2 <- .x[.x$Dataset == pair[2], ]
                    if (nrow(d1) == 0 | nrow(d2) == 0) return(NULL)
                    x <- c(d1$n_hit, d2$n_hit)
                    n <- c(d1$n_total, d2$n_total)
                    p <- tryCatch(prop.test(x, n)$p.value, error = function(e) NA_real_)
                    data.frame(group1 = pair[1], group2 = pair[2],
                            group1_value = d1$hit_rate, group2_value = d2$hit_rate,
                            p_value = p, test = "prop.test")
                }))
                }, .groups = "drop") %>%
                mutate(metric = ifelse(hd == "smry", paste0("HitRate_", summarize_label),
                    paste0("HitRate_", summarize_pairs_label)))
        }
    }

    enr_data <- crispr_sub %>%
        group_by(Dataset) %>%
        summarize(n_total = n(), n_sig = sum(positive_indicator), .groups = "drop") %>%
        right_join(crispr_sub %>%
            group_by(Dataset, category) %>%
            summarize(n_cat = n(), n_sig_cat = sum(positive_indicator), .groups = "drop"),
            by = "Dataset") %>%
        mutate(enrichment = (n_sig_cat / n_cat) / (n_sig / n_total),
            SE_log_enr = sqrt(1 / n_sig_cat - 1 / n_cat + 1 / n_sig - 1 / n_total))

    if (pairs == "Category") {
        message("enrichment pairs category")
      enr_pairs <- enr_data %>%
        group_by(Dataset) %>%
        filter(n_distinct(category) >= 2) %>%
        group_modify(~{
          combs <- combn(unique(.x$category), 2, simplify = FALSE)
          do.call(rbind, lapply(combs, function(pair) {
            d1 <- .x[.x$category == pair[1], ]
            d2 <- .x[.x$category == pair[2], ]
            d <- log(d1$enrichment / d2$enrichment)
            SE_d <- sqrt(d1$SE_log_enr^2 + d2$SE_log_enr^2)
            z <- d / SE_d
            p <- pnorm(-abs(z)) * 2
            data.frame(group1 = pair[1], group2 = pair[2],
                       group1_value = d1$enrichment, group2_value = d2$enrichment,
                       p_value = p, test = "log_enrichment_z")
          }))
        }, .groups = "drop") %>%
        mutate(metric = "Enrichment")
    } else {
        message("enrichment pairs dataset")
      enr_pairs <- enr_data %>%
        group_by(category) %>%
        filter(n_distinct(Dataset) >= 2) %>%
        group_modify(~{
          combs <- combn(unique(.x$Dataset), 2, simplify = FALSE)
          do.call(rbind, lapply(combs, function(pair) {
            d1 <- .x[.x$Dataset == pair[1], ]
            d2 <- .x[.x$Dataset == pair[2], ]
            d <- log(d1$enrichment / d2$enrichment)
            SE_d <- sqrt(d1$SE_log_enr^2 + d2$SE_log_enr^2)
            z <- d / SE_d
            p <- pnorm(-abs(z)) * 2
            data.frame(group1 = pair[1], group2 = pair[2],
                       group1_value = d1$enrichment, group2_value = d2$enrichment,
                       p_value = p, test = "log_enrichment_z")
          }))
        }, .groups = "drop") %>%
        mutate(metric = "Enrichment")
    }

    es_data <- crispr_sub %>%
      filter(!!sym(significance_var) == 1) %>%
      mutate(plotEffectSize = if (significance_var == "Regulated") EffectSize else abs(EffectSize)) %>%
      filter(!is.na(plotEffectSize))

    if (pairs == "Category") {
      es_pairs <- es_data %>%
        group_by(Dataset, category) %>%
        summarize(effect_values = list(plotEffectSize), median_effect = median(plotEffectSize), .groups = "drop") %>%
        group_by(Dataset) %>%
        filter(n_distinct(category) >= 2) %>%
        group_modify(~{
          combs <- combn(unique(.x$category), 2, simplify = FALSE)
          do.call(rbind, lapply(combs, function(pair) {
            d1 <- .x$effect_values[.x$category == pair[1]][[1]]
            d2 <- .x$effect_values[.x$category == pair[2]][[1]]
            p <- wilcox.test(d1, d2)$p.value
            data.frame(group1 = pair[1], group2 = pair[2],
                       group1_value = median(d1), group2_value = median(d2),
                       p_value = p, test = "wilcox")
          }))
        }, .groups = "drop") %>%
        mutate(metric = "EffectSize")
    } else {
      es_pairs <- es_data %>%
        group_by(Dataset, category) %>%
        summarize(effect_values = list(plotEffectSize), median_effect = median(plotEffectSize), .groups = "drop") %>%
        group_by(category) %>%
        filter(n_distinct(Dataset) >= 2) %>%
        group_modify(~{
          combs <- combn(unique(.x$Dataset), 2, simplify = FALSE)
          do.call(rbind, lapply(combs, function(pair) {
            d1 <- .x$effect_values[.x$Dataset == pair[1]][[1]]
            d2 <- .x$effect_values[.x$Dataset == pair[2]][[1]]
            p <- wilcox.test(d1, d2)$p.value
            data.frame(group1 = pair[1], group2 = pair[2],
                       group1_value = median(d1), group2_value = median(d2),
                       p_value = p, test = "wilcox")
          }))
        }, .groups = "drop") %>%
        mutate(metric = "EffectSize")
    }


    hit_pairs <- rbindlist(hit_pairs_list) %>% as.data.frame()
    combined <- bind_rows(hit_pairs, enr_pairs, es_pairs) %>%
      mutate(significance = significance_var) %>%
      group_by(metric, significance) %>%
      mutate(p_value_adj = p.adjust(p_value, method = "BH")) %>%
      ungroup()

    all_results[[significance_var]] <- combined
  }

  final <- bind_rows(all_results)
  out_file <- paste0(out_prefix, group_var, ".signficance_table.by_", pairs, ".tsv")
  fwrite(final, out_file, sep = "\t")

}


read_dc_tap_table_s3 <- function(crispr_path, filter_to_random_DEG = TRUE) {
    res <- fread(crispr_path, sep = "\t")

    if (filter_to_random_DEG) {
        res <- res %>% 
            filter(Random_DistalElement_Gene == TRUE) %>% 
            mutate(Dataset = ifelse(cell_type == "K562", "K562_DC_TAP", "WTC11_DC_TAP"))
    }

    # format columns
    res <- res %>% 
        select(chr = targeting_chr_hg38, start = targeting_start_hg38, end = targeting_end_hg38, cell_type,
            measuredGeneSymbol = gene_symbol, distanceToTSS = distance_to_gencode_gene_TSS,
            direct_vs_indirect_negative, direct_vs_indirect_positive,
            ubiq_category, element_category, Dataset,
            pct_change_effect_size, Significant = significant_wo_pos_controls_20fdr,
            (starts_with("power_at") & ends_with("wo_pos_controls_20fdr"))) %>% 
        mutate(
            elementName = paste0(cell_type, "|", chr, ":", start, "-", end),
            WellPowered = (Significant | power_at_effect_size_15_wo_pos_controls_20fdr >= 0.8), 
            EffectSize = pct_change_effect_size / 100,
            Regulated = (Significant & EffectSize < 0)) %>% 
        distinct() %>% 
        select(-pct_change_effect_size)

    return(res)
}


### MAIN
project_dir <- "/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/2025_0227_CTCF_and_H3K27ac"
results_dir <- file.path(project_dir, "results", "2026_0626_dc_tap_revisions"); dir.create(results_dir, showWarnings = FALSE)

table_s3_file <- file.path(project_dir, "reference", "Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv")

cell_types_cr <- c("GM12878", "HCT116", "Jurkat", "K562_training", "K562_validation", "WTC11")#, "H1_for_WTC11", "H9_for_WTC11")
ct_cp_cr <- c("#e9c54e", "#f29742", "#5eb342", "#792374", "#b778b3", "#5496ce")#, "#00488d", "#00488d")
ds_cp <- c(Gasperini2019 = "#d3a9ce", Nasser2021 = "#b778b3", Schraivogel2020 = "#a64791", K562_DC_TAP = "#006eae", WTC11_DC_TAP = "#00488d")

names(ct_cp_cr) <- cell_types_cr

# --- plot housekeeping gene comparisons, split by cell type --- #
if (TRUE) {
    crispr_merged <- read_dc_tap_table_s3(table_s3_file, filter_to_random_DEG = TRUE) %>% 
        #mutate(Dataset = ifelse(Dataset %in% c("K562_DC_TAP", "WTC11_DC_TAP"), "DC_TAP", Dataset)) %>% 
        mutate(enhancerness = ifelse(element_category %in% c("High H3K27ac", "H3K27ac"), "H3K27ac+ element", "Other element")) %>% 
        mutate(EffectSize = EffectSize * 100) %>% 
        mutate(direct_vs_indirect = ifelse(EffectSize < 0, direct_vs_indirect_negative, direct_vs_indirect_positive))

    group_vars <- c("ubiq_category")
    this_results_dir <- file.path(results_dir, "hkg_comparison"); dir.create(this_results_dir, showWarnings = FALSE)
    
    for (g in group_vars) {
        for (use_allpower in c(TRUE)){
            if (use_allpower) {
                out_dir <- file.path(this_results_dir, "splitDC_TAP_all_metrics_comparison_allPower"); dir.create(out_dir, showWarnings = FALSE)
                out_dir2 <- file.path(this_results_dir, "splitDC_TAP_all_metrics_significance_allPower"); dir.create(out_dir2, showWarnings = FALSE)
            } else {
                out_dir <- file.path(this_results_dir, "splitDC_TAP_all_metrics_comparison_wellPowered"); dir.create(out_dir, showWarnings = FALSE)
                out_dir2 <- file.path(this_results_dir, "splitDC_TAP_all_metrics_significance_wellPowered"); dir.create(out_dir2, showWarnings = FALSE)
            }

            # direct effect weighted
            out_prefix <- paste0(out_dir, "/weighted_")
            out_prefix2 <- paste0(out_dir2, "/weighted_")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = TRUE, direct_effect_threshold = NULL)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = TRUE, direct_effect_threshold = NULL)

            # direct effect filtered
            out_prefix <- paste0(out_dir, "/filter50_")
            out_prefix2 <- paste0(out_dir2, "/filter50_")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)

            # no adjustments for p(direct)
            out_prefix <- paste0(out_dir, "/")
            out_prefix2 <- paste0(out_dir2, "/")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0)
        }
    }
}

# --- plot housekeeping gene comparisons, combined cell type --- #
if (TRUE) {
    crispr_merged <- read_dc_tap_table_s3(table_s3_file, filter_to_random_DEG = TRUE) %>% 
        mutate(Dataset = ifelse(Dataset %in% c("K562_DC_TAP", "WTC11_DC_TAP"), "DC_TAP", Dataset)) %>% 
        mutate(enhancerness = ifelse(element_category %in% c("High H3K27ac", "H3K27ac"), "H3K27ac+ element", "Other element")) %>% 
        mutate(EffectSize = EffectSize * 100) %>% 
        mutate(direct_vs_indirect = ifelse(EffectSize < 0, direct_vs_indirect_negative, direct_vs_indirect_positive))

    group_vars <- c("ubiq_category")
    this_results_dir <- file.path(results_dir, "hkg_comparison"); dir.create(this_results_dir, showWarnings = FALSE)
    
    for (g in group_vars) {
        for (use_allpower in c(TRUE)){
            if (use_allpower) {
                out_dir <- file.path(this_results_dir, "combinedDC_TAP_all_metrics_comparison_allPower"); dir.create(out_dir, showWarnings = FALSE)
                out_dir2 <- file.path(this_results_dir, "combinedDC_TAP_all_metrics_significance_allPower"); dir.create(out_dir2, showWarnings = FALSE)
            } else {
                out_dir <- file.path(this_results_dir, "combinedDC_TAP_all_metrics_comparison_wellPowered"); dir.create(out_dir, showWarnings = FALSE)
                out_dir2 <- file.path(this_results_dir, "combinedDC_TAP_all_metrics_significance_wellPowered"); dir.create(out_dir2, showWarnings = FALSE)
            }

            # direct effect weighted
            out_prefix <- paste0(out_dir, "/weighted_")
            out_prefix2 <- paste0(out_dir2, "/weighted_")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = TRUE, direct_effect_threshold = NULL)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = TRUE, direct_effect_threshold = NULL)

            # direct effect filtered
            out_prefix <- paste0(out_dir, "/filter50_")
            out_prefix2 <- paste0(out_dir2, "/filter50_")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)

            # no adjustments for p(direct)
            out_prefix <- paste0(out_dir, "/")
            out_prefix2 <- paste0(out_dir2, "/")
            plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0)
            compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
                direct_effect_weighted = FALSE, direct_effect_threshold = 0)
        }
    }
}