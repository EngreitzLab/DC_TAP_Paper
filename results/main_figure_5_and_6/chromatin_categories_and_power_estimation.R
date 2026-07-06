### --- DEFINE CHROMATIN CATEGORIES ---
# get thresholds based on genome-wide elements
get_category_thresholds <- function(enh, quantiles) {
    enh_ctcf <- enh %>% filter(CTCF_peak_overlap == 1) %>% 
        select(cell_type, CTCF.H3K27ac.ratio.CTCF_peak = CTCF.H3K27ac.ratio) %>% 
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))
    
    enh_h3k27ac <- enh %>% filter(H3K27ac_peak_overlap == 1) %>% 
        select(cell_type, H3K27ac.RPM.H3K27ac_peak = H3K27ac.RPM) %>% 
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    enh_h3k27me3 <- enh %>% filter(H3K27me3_peak_overlap == 1) %>% 
        select(cell_type, H3K27me3.RPM.expandedRegion.H3K27me3_peak = H3K27me3.RPM.expandedRegion) %>% 
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    enh_other <- enh %>% 
        select(cell_type, CTCF.RPM, H3K27ac.RPM, H3K27ac.RPM.expandedRegion, DHS.RPM) %>%
        pivot_longer(cols = -cell_type, names_to = "feature", values_to = "value") %>% 
        group_by(cell_type, feature) %>% 
        reframe(quantile = quantiles, value = quantile(value, probs = quantiles, na.rm = TRUE))

    res <- rbind(enh_ctcf, enh_h3k27ac, enh_h3k27me3, enh_other)
    
    return(res)
}

# get table of quantile values from genome-wide elements
get_threshold_key <- function(thresholds, feature_col, quantile_this) {
    filt <- thresholds %>% filter(feature == feature_col, quantile == quantile_this)
    key <- setNames(filt$value, filt$cell_type)
    return(key)
}

# categorize elements
categorize_elements <- function(enh, thresholds, H3K27ac_q_high = 0.9, H3K27ac_q_low = 0.5) {
    ### CATEGORIZATION LOGIC ###
    # if element overlaps H3K27ac peak:
        # if H3K27ac.RPM.expandedRegion > 90% --> High H3K27ac
        # if H3K27ac.RPM.expandedRegion < 90% --> H3K27ac
    # else:
        # if H3K27ac.RPM.expandedRegion > 90% --> High H3K27ac
        # if H3K27ac.RPM.expandedRegion > 50% --> H3K27ac
        # if element overlaps CTCF peak --> CTCF element
        # if element overlaps H3K27me3 peak --> H3K27me3 element
        # else: No H3K27ac
        

    key_high <- get_threshold_key(thresholds, "H3K27ac.RPM.expandedRegion", H3K27ac_q_high)
    key_low <- get_threshold_key(thresholds, "H3K27ac.RPM.expandedRegion", H3K27ac_q_low)

    enh <- enh %>%
         mutate(element_category = case_when(
                H3K27ac_peak_overlap == 1 & H3K27ac.RPM.expandedRegion >= key_high[cell_type] ~ "High H3K27ac",
                H3K27ac_peak_overlap == 1 ~ "H3K27ac",
                H3K27ac.RPM.expandedRegion >= key_high[cell_type] ~ "High H3K27ac",
                H3K27ac.RPM.expandedRegion >= key_low[cell_type] ~ "H3K27ac",
                CTCF_peak_overlap == 1 ~ "CTCF element",
                H3K27me3_peak_overlap == 1 ~ "H3K27me3 element",
                TRUE ~ "No H3K27ac"))
    return(enh)
}

# summarize properties of genome-wide element-gene pairs
annotate_genomewide_pairs <- function(enh, e2g_files, cell_types, remove_promoters, distance_threshold) {
    res_list <- vector("list", length(cell_types))
    res_list_genes <- vector("list", length(cell_types))

    for (i in seq_along(cell_types)) {
        ct <- cell_types[i]
        print(ct)

        enh_ct <- filter(enh, cell_type == ct)
        pairs_file <- e2g_files[ct]; print(pairs_file)

        pairs <- fread(pairs_file, sep = "\t") %>% 
            select(chr, start, end, class, TargetGene, ubiquitousExpressedGene, distanceToTSS, CellType) %>% 
            filter(distanceToTSS < distance_threshold) %>%
            mutate(distance_category = case_when(distanceToTSS < 10e3 ~  "0-10 kb",
                                             distanceToTSS < 100e3 ~ "10-100 kb",
                                             distanceToTSS < 250e3 ~ "100-250 kb",
                                             distanceToTSS < 1000e3 ~ "250 kb-1 Mb",
                                             distanceToTSS < 2000e3 ~ "1 Mb-2 Mb",
                                             TRUE ~ ">2 Mb"),
                    ubiq_category = ifelse(ubiquitousExpressedGene %in% c("True", TRUE), "Ubiq. expr. gene", "Other gene"))

        if (remove_promoters) {
            pairs <- filter(pairs, class != "promoter")
        }

        res_list[[i]] <- left_join(pairs, enh_ct, by = c("chr", "start", "end")) %>% 
            group_by(cell_type, element_category, distance_category, ubiq_category) %>% 
            summarize(n_pairs = n())

        res_list_genes[[i]] <- pairs %>% select(cell_type = CellType, TargetGene, ubiq_category) %>% distinct() %>% 
            group_by(cell_type, ubiq_category) %>%
            summarize(n_e2g_genes = n())
    }

    res <- rbindlist(res_list) %>% as.data.frame()
    res_genes <- rbindlist(res_list_genes) %>% as.data.frame()

    return(list(res, res_genes))
}

### --- PLOT PROPORTIONS OF ELEMENT-GENE PAIRS ---
# helper function to format numbers nicely
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

# plot proportions 
plot_category_proportion_by_dataset_sceptre <- function(enh, crispr, category_col = "element_category",
    direct_effect_weighted = FALSE, direct_effect_threshold = 0.5) {

    # plotting params
    if (category_col == "element_category") {
        category_names <- c("H3K27me3 element", "CTCF element", "High H3K27ac", "H3K27ac", "No H3K27ac")
        cp <- c("#429130", "#49bcbc", "#c5373d", "#d9694a", "#c5cad7")
        names(cp) <- category_names
     } else if (category_col == "distance_category") {
        cp <- c("#002359", "#00488d", "#006eae", "#5496ce", "#9bcae9", "#9bcae9")
        names(cp) <- c("0-10 kb", "10-100 kb", "100-250 kb", "250 kb-1 Mb", "> 1 Mb", "1 Mb-2 Mb")
    } else if (category_col == "ubiq_category") {
        cp <- c("#792374", "#b778b3")
        names(cp) <- c("Ubiq. expr. gene", "Other gene")
    } 

    ct_include <- c("K562_WTC11", unique(crispr$cell_type))
    ds_include <- c("Gasperini2019", "Nasser2021", "Schraivogel2020", "K562_DC_TAP", "WTC11_DC_TAP", "DC_TAP")

    cat_order <- c(gw = "Genome-wide\nE-G pairs", crispr_tested = "All CRISPR\ntested pairs", crispr_powered = "Well-powered\nCRISPR pairs",
        crispr_pos_downreg = "CRISPR+\ndownregulated", crispr_pos_downreg_weighted = "Weighted CRISPR+\ndownregulated",
        crispr_pos_upreg ="CRISPR+\nupregulated", crispr_pos_upreg_weighted ="Weighted CRISPR+\nupregulated")
    ct_order <- c(K562_WTC11 = "K562 + WTC11", K562_WTC11_validation = "K562 + WTC11 (DC-TAP)",
        K562 = "K562", K562_training = "K562 (previous)", K562_validation = "K562 (DC-TAP)",
        WTC11 = "WTC11",  WTC11_validation = "WTC11 (DC-TAP)")

    text_angle <- 45

    # prepare genome-wide prediciton data
    enh <- enh %>%
        filter(cell_type %in% ct_include) %>% 
        mutate(cell_type_label = cell_type, category_use = !!sym(category_col)) %>% 
        group_by(cell_type_label, category_use) %>% 
        summarize(n_pairs = sum(n_pairs), .groups = "drop") %>%
        mutate(plot_category = cat_order["gw"])

    # all crispr
    crispr_all <- crispr %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        select(cell_type_label, category_use) %>% 
        mutate(plot_category = cat_order["crispr_tested"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = n(), .groups = "drop")
    
    # crispr well-powered
    crispr_wp <- crispr %>% 
        filter(WellPowered == TRUE) %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        select(cell_type_label, category_use) %>% 
        mutate(plot_category = cat_order["crispr_powered"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = n(), .groups = "drop")

    # crispr pos
    crispr_pos_downreg <- crispr %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        filter(WellPowered == TRUE, Regulated == 1, direct_vs_indirect > direct_effect_threshold) %>%
        select(cell_type_label, category_use) %>% 
        mutate(plot_category = cat_order["crispr_pos_downreg"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = n(), .groups = "drop")

    crispr_pos_downreg_weighted <- crispr %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        filter(WellPowered == TRUE, Regulated == 1) %>%
        select(cell_type_label, category_use, direct_vs_indirect, EffectSize) %>% 
        mutate(plot_category = cat_order["crispr_pos_downreg_weighted"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = sum(direct_vs_indirect), .groups = "drop")
    
    crispr_pos_upreg <- crispr %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        filter(WellPowered == TRUE, Significant == TRUE, EffectSize > 0, direct_vs_indirect > direct_effect_threshold) %>%
        select(cell_type_label, category_use) %>% 
        mutate(plot_category = cat_order["crispr_pos_upreg"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = n(), .groups = "drop")
    
    crispr_pos_upreg_weighted <- crispr %>% 
        mutate(cell_type_label = paste0(cell_type, "_", data_category), category_use = !!sym(category_col)) %>%
        filter(WellPowered == TRUE, Significant == TRUE, EffectSize > 0) %>%
        select(cell_type_label, category_use, direct_vs_indirect, EffectSize) %>% 
        mutate(plot_category = cat_order["crispr_pos_upreg_weighted"]) %>% 
        group_by(cell_type_label, category_use, plot_category) %>% 
        summarize(n_pairs = sum(direct_vs_indirect), .groups = "drop")
    
    if (direct_effect_weighted) {
        crispr_pos_downreg <- crispr_pos_downreg_weighted
        crispr_pos_upreg <- crispr_pos_upreg_weighted
    }

    res <- rbind(enh, crispr_all, crispr_wp, crispr_pos_downreg) %>% #, crispr_pos_upreg) %>% 
        mutate(category_use = factor(category_use, levels = names(cp), ordered = TRUE),
            plot_category = factor(plot_category, levels = cat_order, ordered = TRUE),
            cell_type_label = ct_order[cell_type_label],
            cell_type_label = factor(cell_type_label, levels = ct_order, ordered = TRUE))

    # calculate totals per category
    dataset_totals <- res %>%
        group_by(plot_category, cell_type_label) %>%
        summarise(total_count = sum(n_pairs), .groups = "drop")

    # calculate proportions and prepare data for plotting

    prop_data <- res %>%
        #group_by(plot_category, cell_type_label, element_category) %>%
        mutate(category_count = n_pairs) %>%
        left_join(dataset_totals, by = c("plot_category", "cell_type_label")) %>%
        mutate(proportion = category_count / total_count,
            y_position = cumsum(proportion) - proportion / 2,
            category_use = factor(category_use, levels = names(cp), ordered = TRUE),
            plot_category = factor(plot_category, levels = cat_order, ordered = TRUE),
            cell_type_label = factor(cell_type_label, levels = ct_order, ordered = TRUE),
            category_count = format_number(category_count),
            total_count = format_number(total_count))
            #category_count = ifelse(category_count < 100e3, format_int(category_count), format_sci(category_count)),
            #total_count = ifelse(total_count < 100e3, format_int(total_count), format_sci(total_count)))

    # create stacked bar plot
    x_lab <- ""

    p <- ggplot(prop_data, aes(x = plot_category, y = proportion, fill = category_use)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = category_count), position = position_stack(vjust = 0.5), color = "#000000", size = 2.5) +
        geom_text(aes(y = -0.05, label = paste0("(", total_count, ")")), size = 3, color = "#000000", vjust = 0) +
        labs(x = x_lab, fill = "Element category", y = "Proportion of pairs with element in category") +
        facet_grid(. ~ cell_type_label, scales = "free", space = "free") +
        scale_fill_manual(values = cp) +
        theme_classic() + theme(strip.background = element_blank(), panel.grid = element_blank(),
            axis.text = element_text(size = 10, color = "#000000"), axis.text.x = element_text(angle = text_angle, hjust = 1, vjust = 1),
            axis.title = element_text(size = 12), axis.ticks = element_line(color = "#000000"), legend.position = "right",
            plot.title = element_text(size = 14, color = "#000000"), plot.subtitle = element_text(size = 12, color = "#000000"))

    return(p)
}

### --- ESTIMATE UNDETECTED VERSUS DETECTED POSITIVE SIGNIFICANT PAIRS --- 
plot_positives_by_effect_size <- function(crispr, es_bins, bin_min, bin_max, significance_var, out_prefix) {
    # stacked barplot with x-axis = effect size bins, y axis = # positive DE-G pairs
    # categories for bars: "statistically significant" vs "not statistically significant"
    # second category estimated based on rate of positives and sum of power for all tested pairs for low end of effect size
    
    ## params
    power_cols <- paste0("power_at_effect_size_", bin_min); names(power_cols) <- es_bins
    category_key <- c("Statistically significant" = "#435369", "Not detected as significant" = "#c5cad7")
    dataset_key <- c(all_DC_TAP = "DC-TAP", Gasperini2019 = "Gasperini et al.")

    ## calculate total power and pairs per dataset
    bin_sumTestedPower <- crispr %>%
        mutate(!!sym(power_cols[1]) := 0) %>%
        group_by(Dataset) %>%
        mutate(nTested_total = n()) %>% 
        group_by(Dataset, nTested_total) %>%
        summarize(across(all_of(power_cols), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
        pivot_longer(cols = all_of(names(power_cols)),
            names_to = "EffectSize_bin",
            values_to = "sum_testedPower")

    ## summarize results
    res <- crispr %>% 
        mutate(abs_EffectSize = abs(EffectSize),
            detectedPositive = !!sym(significance_var),
            EffectSize_bin = cut(abs_EffectSize, breaks = c(bin_min[1], bin_max), labels = es_bins,
                right = TRUE, include.lowest = TRUE),
            EffectSize_bin = as.character(EffectSize_bin)) %>%
        group_by(Dataset, EffectSize_bin) %>% 
        summarize(nDetectedPositive = sum(detectedPositive),
            nTested_bin = n(),
            .groups = "drop") %>%
        left_join(bin_sumTestedPower, by = c("Dataset", "EffectSize_bin")) %>% 
        mutate(positiveRate = nDetectedPositive / sum_testedPower,
            nTotalPositive = positiveRate * nTested_total,
            nUndetectedPositive = ifelse(EffectSize_bin != es_bins[1], nTotalPositive - nDetectedPositive, 0),
            EffectSize_bin = factor(EffectSize_bin, levels = es_bins, ordered = TRUE)) %>%
        ungroup() %>%
        arrange(Dataset, EffectSize_bin)

    fwrite(res, paste0(out_prefix, "positives_by_effect_size.", significance_var, ".tsv"), sep = "\t", quote = FALSE)

    ## reformat for plotting
    res_long <- res %>%
        select(Dataset, EffectSize_bin, nDetectedPositive, nUndetectedPositive) %>%
    pivot_longer(cols = c(nDetectedPositive, nUndetectedPositive),
        names_to = "DetectionCategory", values_to = "n") %>%
    mutate(DetectionCategory = ifelse(DetectionCategory == "nDetectedPositive", names(category_key)[1], names(category_key)[2]),
        DetectionCategory = factor(DetectionCategory, levels = rev(names(category_key)), ordered = TRUE),
        DatasetLabel = dataset_key[Dataset],
        DatasetLabel = factor(DatasetLabel, levels = unique(dataset_key), ordered = TRUE))

    # plot
    p1 <- ggplot(res_long, aes(x = EffectSize_bin, y = n, fill = DetectionCategory)) +
        geom_bar(stat = "identity") +
        facet_wrap(~DatasetLabel, ncol = 2, scales = "free_y") +
        scale_fill_manual(values = category_key, name = NULL) +
        labs(x = "", y = "# of positive enhancer-gene pairs", fill = "Detected?") +
        theme_classic() + theme(
            strip.background = element_blank(), panel.grid = element_blank(), strip.text = element_text(size = 12, color = "#000000"),
            axis.text = element_text(size = 10, color = "#000000"), axis.text.x = element_blank(), #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 10), axis.ticks = element_line(color = "#000000"), 
            legend.position = "top")

    res <- res %>% 
        mutate(DatasetLabel = dataset_key[Dataset],
            DatasetLabel = factor(DatasetLabel, levels = unique(dataset_key), ordered = TRUE),
            positiveRate = ifelse(sum_testedPower == 0, 0, positiveRate))

    p2 <- ggplot(res, aes(x = EffectSize_bin, y = positiveRate)) +
        geom_bar(stat = "identity", fill = "#1c2a43") +
        facet_wrap(~DatasetLabel, ncol = 2, scales = "fixed") +
        labs(x = "Effect size bin", y = "Positive rate\n(Detected positives / sum (power of tested pairs))") +
        theme_classic() + theme(
            strip.background = element_blank(), panel.grid = element_blank(), strip.text = element_blank(),
            axis.text = element_text(size = 10, color = "#000000"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 10), axis.ticks = element_line(color = "#000000"),
            legend.position = "none")

    gr <- plot_grid(p1, p2, nrow = 2, align = "hv", rel_heights = c(1, 1))
    ggsave2(paste0(out_prefix, "positives_by_effect_size_and_rate.", significance_var, ".pdf"), gr, height = 6, width = 7)
}

### --- DATA I/O ---
# read in crispr data processed by dc-tap pipeline and annotated with chrom-annotate pipeline
read_sceptre_crispr_data <- function(crispr_path, direct_effect_path, filter_to_random_DEG = TRUE, annot = TRUE) {
    de <- fread(direct_effect_path) %>% 
        select(element_gene_pair_identifier_hg38, cell_type,
            direct_vs_indirect_negative, direct_vs_indirect_positive, distance_category, old_element_category = element_category) %>%
        distinct()

    res <- fread(crispr_path, sep = "\t") %>%
        mutate(cell_type = ifelse(cell_type == "WTC11_for_WTC11", "WTC11", cell_type)) %>%
        left_join(de, by = c("cell_type", "element_gene_pair_identifier_hg38")) %>% 
        filter(!is.na(direct_vs_indirect_negative))
    print(nrow(res))


    if (length(unique(res$cell_type)) == 2) {
        # DC-TAP
            temp_init <- res %>% filter(Random_DistalElement_Gene == TRUE) %>% filter(significant | power_at_effect_size_15 >= 0.8)
            print(nrow(temp_init))

        if (filter_to_random_DEG) {
            res <- res %>% 
                filter(Random_DistalElement_Gene == TRUE) %>% 
                mutate(Dataset = ifelse(cell_type == "K562", "K562_DC_TAP", "WTC11_DC_TAP"))
        }
        res <- res %>% mutate(data_category = "DC_TAP_SCEPTRE") 
    } else {
        # Gasperini
        if (filter_to_random_DEG) {
            res <- res %>% 
                filter(DistalElement_Gene == TRUE) 
        }
        res <- res %>% 
            mutate(Dataset = "Gasperini2019",
                data_category = "Gasperini_SCEPTRE")
    }

    # format columns
    res <- res %>% 
        select(chr = targeting_chr_hg38, start = targeting_start_hg38, end = targeting_end_hg38, cell_type,
            measuredGeneSymbol = gene_symbol, distanceToTSS = distance_to_gencode_gene_TSS,
            direct_vs_indirect_negative, direct_vs_indirect_positive,
            distance_category, old_element_category, #any_of(c("ubiq_category")),
            pct_change_effect_size, Significant = significant, Dataset, starts_with("power_at"),
            ends_with("peak_overlap"), ends_with(".RPM"), ends_with(".RPM.expandedRegion")) %>% 
        mutate(
            elementName = paste0(cell_type, "|", chr, ":", start, "-", end),
            WellPowered = (Significant | power_at_effect_size_15 >= 0.8), 
            EffectSize = pct_change_effect_size / 100,
            Regulated = (Significant & EffectSize < 0)) %>% 
        distinct() %>% 
        select(-pct_change_effect_size)

    if (annot) {
        res <- res %>% mutate(CTCF.H3K27ac.ratio = (CTCF.RPM) / (H3K27ac.RPM + 0.001))
    }


    return(res)
}

# read in genome-wide element lists annotated with chrom-annotate pipeline
read_enh_lists <- function(file_path, remove_promoters) {
    enh <- fread(file_path, sep = "\t") 
    if (remove_promoters) {
        enh <- filter(enh, class != "promoter")
    }

    return(enh)
}



### --- DEFINE FILE PATHS AND PARAMS ---

# download from: https://github.com/EngreitzLab/ENCODE_rE2G/blob/dev/resources/external_features/gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv
promoter_class_path <- "ENCODE_rE2G/resources/external_features/gene_promoter_class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.tsv"

# download from: https://github.com/EngreitzLab/ENCODE_rE2G/blob/dev/reference/CollapsedGeneBounds.hg38.TSS500bp.bed"
abc_genes_path <- "ENCODE_rE2G/reference/CollapsedGeneBounds.hg38.TSS500bp.bed"

## inputs for combining all annotations
# parameters for categorization
remove_promoters <- TRUE
quantiles <- c(0.1, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.75, 0.85, 0.9, 0.95)
cell_types <- c("K562", "GM12878", "HCT116", "Jurkat", "WTC11")

# results from chrom-annotate pipeline
chrom_annotate_results <- "chrom_annotate_results"
enh_lists <- c(file.path(chrom_annotate_results, cell_types, "EnhancerList.extended.tsv"))
names(enh_lists) <- cell_types

crispr_annot_files <- c(other_crispr = file.path(chrom_annotate_results, "CRISPR_data", "other_crispr_screens.tsv"),
    dc_tap_gasperini = file.path(chrom_annotate_results, "CRISPR_data", "dc_tap_gasperini.random_distal_element_pairs.tsv"))

# results from direct effect rate annotations
direct_effect_results <- "direct_effect_results"
crispr_annot_files <- c(other_crispr = file.path(direct_effect_results, "other_crispr_screens.direct_effect.tsv.gz")
    dc_tap_gasperini = file.path(direct_effect_results, "dc_tap_gasperini.random_distal_element_pairs.direct_effect.tsv.gz"))

# results from encode_re2g pipeline to get all candidate element-gene pairs
e2g_base <- "ENCODE_rE2G_results"
e2g_res <- c(GM12878 = "2025_0227_validation_new_inputs", HCT116 = "2025_0227_validation_new_inputs", K562 = "2025_0227_validation_new_inputs",
    Jurkat = "2025_0227_validation_new_inputs", WTC11 = "2025_0227_validation_new_inputs")
e2g_files <- lapply(cell_types, function(ct) file.path(e2g_base, e2g_res[ct], paste0(ct, "_H3K27ac_megamap"),
                                                "dhs_h3k27ac_megamap", "encode_e2g_predictions.tsv.gz")) %>% unlist() %>% setNames(cell_types)

## categorized pairs file paths 
recat_dir <- file.path("categorized_pairs")
crispr_encode_categorized_file <- file.path(recat_dir, "other_crispr_screens.categorized.tsv.gz")
crispr_sceptre_categorized_file <- file.path(recat_dir, "dc_tap_gasperini.random_distal_element_pairs.categorized.tsv.gz")
thresholds_file <- file.path(recat_dir, "thresholds.tsv")


### --- ANNOTATE DATA WITH CHROMATIN CATEGORIES ---
enh <- lapply(enh_lists, read_enh_lists, remove_promoters) %>% 
    rbindlist(idcol = "cell_type") %>% as.data.frame() %>%
    mutate(CTCF.H3K27ac.ratio = (CTCF.RPM) / (H3K27ac.RPM + 0.001)) %>% 
    replace(is.na(.), 0)
message("Read in candidate elements!")

thresholds <- get_category_thresholds(enh, quantiles) %>% arrange(cell_type, feature)
fwrite(thresholds, thresholds_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
message("Saved thresholds!")

dc_tap_gasperini <- read_sceptre_crispr_data(crispr_annot_files[["dc_tap_gasperini"]], direct_effect_files[["dc_tap_gasperini"]])
other_crispr <- read_sceptre_crispr_data(crispr_annot_files[["other_crispr"]], direct_effect_files[["other_crispr"]])
message("Read in and formatted CRISPR data!")

ubiq_expr_genes <- fread(promoter_class_path, sep = "\t") %>% filter(is_ubiquitous_uniform %in% c("True", TRUE)) %>% pull(TargetGene)
crispr_sceptre <- dc_tap_gasperini %>%
        categorize_elements(thresholds, H3K27ac_q_high = 0.9, H3K27ac_q_low = 0.5) %>% 
        mutate(ubiq_category = ifelse(measuredGeneSymbol %in% ubiq_expr_genes, "Ubiq. expr. gene", "Other gene"))
fwrite(crispr_sceptre, crispr_sceptre_categorized_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
other_crispr <- other_crispr %>% 
        categorize_elements(thresholds, H3K27ac_q_high = 0.9, H3K27ac_q_low = 0.5) %>% 
        mutate(ubiq_category = ifelse(measuredGeneSymbol %in% ubiq_expr_genes, "Ubiq. expr. gene", "Other gene"))
fwrite(crispr_encode, crispr_encode_categorized_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
message("Categorized CRISPR data!")

enh_cat <- enh %>% 
    select(chr, start, end, cell_type, H3K27me3_peak_overlap, CTCF_peak_overlap, H3K27ac_peak_overlap, H3K27ac.RPM.expandedRegion) %>% 
    categorize_elements(thresholds, H3K27ac_q_high = 0.9, H3K27ac_q_low = 0.5)

## summarize by category
enh_pairs_cat <- annotate_genomewide_pairs(enh_cat, e2g_files, cell_types, remove_promoters, 2e6)
gw_pairs <- enh_pairs_cat[[1]]
fwrite(gw_pairs, file.path(recat_dir, "all_gw_pairs.categorized.summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
message("Categorized genome-wide pairs!")

### --- PLOT CHROMATIN CATEGORY PROPORTIONS ---
crispr_sceptre_use <- crispr_sceptre %>%
        mutate(data_category = ifelse(Dataset == "Gasperini2019", "training", "validation")) %>% 
        mutate(direct_vs_indirect = ifelse(EffectSize < 0, direct_vs_indirect_negative, direct_vs_indirect_positive))

# DC-TAP by cell type (Fig. S5)
z <- plot_category_proportion_by_dataset_sceptre(gw_pairs, crispr_sceptre_use, "element_category",
        direct_effect_weighted = TRUE, direct_effect_threshold = 0)
    ggsave(file.path(recat_dir, "proportion_by_cell_type.element_category.de_weighted.for_dc_tap.pdf"), z, height = 6, width = 8)

# DC-TAP combined (Fig. 5b)
crispr_dc_combined <- crispr_sceptre %>%
    mutate(Dataset = ifelse(Dataset %in% c("K562_DC_TAP", "WTC11_DC_TAP"), "DC_TAP", Dataset)) %>% 
    mutate(cell_type = ifelse(Dataset == "DC_TAP", "K562_WTC11", cell_type))
gw_k562 <- gw_pairs %>% filter(cell_type == "K562")    
gw_pairs_combined <- gw_pairs %>% mutate(cell_type = ifelse(cell_type %in% c("K562", "WTC11"), "K562_WTC11", cell_type)) %>% 
    rbind(gw_k562)
z <- plot_category_proportion_by_dataset_sceptre(gw_pairs_combined, crispr_dc_combined, "element_category",
    direct_effect_weighted = TRUE, direct_effect_threshold = 0)
ggsave(file.path(recat_dir, "proportion_by_cell_type.element_category.de_weighted.combine_dc_tap.pdf"), z, height = 6, width = 7)


### -- COMPUTE RESULTS BY PROMOTER CLASS (Fig. 5c,d) ---
this_results_dir <- "housekeeping_genes"
crispr_sceptre_use <- crispr_sceptre %>%
    mutate(Dataset = ifelse(Dataset %in% c("K562_DC_TAP", "WTC11_DC_TAP"), "DC_TAP", Dataset)) %>% 
    mutate(EffectSize = EffectSize * 100) %>% 
    mutate(direct_vs_indirect = ifelse(EffectSize < 0, direct_vs_indirect_negative, direct_vs_indirect_positive))
group_var <- "ubiq_category"
for (use_allpower in c(TRUE, FALSE)){
    if (use_allpower) {
        out_dir <- file.path(this_results_dir, "v3_combinedDC_TAP_all_metrics_comparison_allPower"); dir.create(out_dir, showWarnings = FALSE)
        out_dir2 <- file.path(this_results_dir, "v3_combinedDC_TAP_all_metrics_significance_allPower"); dir.create(out_dir2, showWarnings = FALSE)
    } else {
        out_dir <- file.path(this_results_dir, "v3_combinedDC_TAP_all_metrics_comparison_wellPowered"); dir.create(out_dir, showWarnings = FALSE)
        out_dir2 <- file.path(this_results_dir, "v3_combinedDC_TAP_all_metrics_significance_wellPowered"); dir.create(out_dir2, showWarnings = FALSE)
    }

    # direct effect weighted
    out_prefix <- paste0(out_dir, "/weighted_"); out_prefix2 <- paste0(out_dir2, "/weighted_")
    plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
        direct_effect_weighted = TRUE, direct_effect_threshold = NULL)
    compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
        direct_effect_weighted = TRUE, direct_effect_threshold = NULL)

    # direct effect filtered
    out_prefix <- paste0(out_dir, "/filter50_"); out_prefix2 <- paste0(out_dir2, "/filter50_")
    plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
        direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)
    compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
        direct_effect_weighted = FALSE, direct_effect_threshold = 0.5)

    # no adjustments for p(direct)
    out_prefix <- paste0(out_dir, "/"); out_prefix2 <- paste0(out_dir2, "/")
    plot_percent_positive_enrichment_effect_size_combined(crispr = crispr_merged, group_var = g, out_prefix = out_prefix, all_power = use_allpower,
        direct_effect_weighted = FALSE, direct_effect_threshold = 0)
    compare_metrics_across_groups(crispr = crispr_merged, group_var = g, out_prefix = out_prefix2, pairs = "Category", all_power = use_allpower,
        direct_effect_weighted = FALSE, direct_effect_threshold = 0)
}

### --- ESTIMATE UNDETECTED VS DETECTED SIGNIFICANT POSITIVES ---
this_out <- file.path(this_results_dir, "detected_positives"); dir.create(this_out, showWarnings = FALSE)
out_prefix <- paste0(this_out, "/5_10_")
es_bins <- c("<=5%", "(5%, 10%]", "(10%, 15%]", "(15%, 25%]", "(25%, 50%]", "(50%, Inf)")
bin_min <- c(0, 5, 10, 15, 25, 50); names(bin_min) <- es_bins
bin_max <- c(5, 10, 15, 25, 50, Inf); names(bin_max) <- es_bins        
plot_positives_by_effect_size(crispr_sceptre, es_bins, bin_min, bin_max, "Regulated", out_prefix)

