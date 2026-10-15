list_filtered_h5_files <- function(h5_dir) {
	h5_files <- list.files(
		h5_dir,
		pattern = "_filtered_feature_bc_matrix\\.h5$",
		full.names = TRUE
	)

	stopifnot(length(h5_files) > 0)

	sample_ids <- stringr::str_remove(
		basename(h5_files),
		"_filtered_feature_bc_matrix\\.h5$"
	)
	names(h5_files) <- sample_ids

	h5_files
}

load_terra_metadata <- function(terra_metadata_path, sample_ids) {
	terra_metadata <- data.table::fread(terra_metadata_path)
	terra_metadata[, sample_id := `entity:validation_id`]

	missing_samples <- setdiff(sample_ids, terra_metadata$sample_id)
	if (length(missing_samples) > 0) {
		stop(
			"Missing Terra metadata for sample(s): ",
			paste(missing_samples, collapse = ", ")
		)
	}

	if (anyDuplicated(terra_metadata$sample_id) > 0) {
		duplicated_samples <- unique(terra_metadata$sample_id[duplicated(terra_metadata$sample_id)])
		stop(
			"Duplicate Terra metadata rows for sample(s): ",
			paste(duplicated_samples, collapse = ", ")
		)
	}

	terra_metadata
}

extract_modality_counts <- function(h5_file, modality_pattern, use_names) {
	counts <- Seurat::Read10X_h5(
		filename = h5_file,
		use.names = use_names,
		unique.features = TRUE
	)

	if (!is.list(counts)) {
		if (tolower(modality_pattern) == "gene expression") {
			return(counts)
		}

		return(NULL)
	}

	assay_name <- names(counts)[tolower(names(counts)) == tolower(modality_pattern)]
	if (length(assay_name) == 0) {
		return(NULL)
	}

	counts[[assay_name[1]]]
}

normalize_feature_names <- function(feature_names) {
	gsub("_", "-", feature_names, fixed = TRUE)
}

prepare_feature_metadata <- function(feature_metadata, feature_names) {
	normalized_metadata <- feature_metadata
	rownames(normalized_metadata) <- normalize_feature_names(rownames(normalized_metadata))

	if (anyDuplicated(rownames(normalized_metadata)) > 0) {
		stop("Feature metadata names are not unique after underscore normalization.")
	}

	missing_features <- setdiff(feature_names, rownames(normalized_metadata))
	if (length(missing_features) > 0) {
		stop(
			"Cannot align feature metadata for feature(s): ",
			paste(head(missing_features, 10), collapse = ", ")
		)
	}

	normalized_metadata[feature_names, , drop = FALSE]
}

build_feature_metadata <- function(template_counts, template_counts_named, id_name, symbol_name) {
	if (is.null(template_counts) || is.null(template_counts_named)) {
		stop("Missing template counts needed to build feature metadata.")
	}

	if (!identical(dim(template_counts), dim(template_counts_named))) {
		stop("Feature IDs and names do not align in the template H5 file.")
	}

	feature_metadata <- data.frame(
		row.names = rownames(template_counts),
		stringsAsFactors = FALSE
	)
	feature_metadata[[id_name]] <- rownames(template_counts)
	feature_metadata[[symbol_name]] <- rownames(template_counts_named)

	feature_metadata
}

load_readout_gene_ids <- function(table_s4_path) {
	table_s4 <- data.table::fread(table_s4_path)
	unique(table_s4[["ensemble_id"]])
}

load_gtf_gene_types <- function(gtf_path, cache_path = NULL) {
	if (!requireNamespace("rtracklayer", quietly = TRUE)) {
		stop("Package rtracklayer is required to import the GTF annotation.")
	}

	if (!is.null(cache_path) && file.exists(cache_path)) {
		cached_gene_types <- data.table::fread(cache_path)
	} else {
		gtf <- rtracklayer::import(gtf_path)
		gene_rows <- gtf[S4Vectors::mcols(gtf)$type == "gene"]
		cached_gene_types <- data.table::data.table(
			gene_id = as.character(S4Vectors::mcols(gene_rows)$gene_id),
			gene_type = as.character(S4Vectors::mcols(gene_rows)$gene_type)
		)
		cached_gene_types <- unique(cached_gene_types[!is.na(gene_id) & nzchar(gene_id)])

		if (!is.null(cache_path)) {
			dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
			data.table::fwrite(cached_gene_types, cache_path, sep = "\t")
		}
	}

	if (!all(c("gene_id", "gene_type") %in% names(cached_gene_types))) {
		stop("Gene type cache must contain gene_id and gene_type columns.")
	}

	unique(cached_gene_types)
}

fetch_gene_types_from_gtf <- function(gene_ids, gtf_path, cache_path = NULL) {
	gene_ids <- unique(stats::na.omit(as.character(gene_ids)))
	gene_ids <- gene_ids[nzchar(gene_ids)]
	gene_type_table <- load_gtf_gene_types(gtf_path = gtf_path, cache_path = cache_path)
	gene_type_lookup <- as.data.frame(gene_type_table)
	rownames(gene_type_lookup) <- gene_type_lookup$gene_id
	gene_type_lookup[gene_ids, c("gene_id", "gene_type"), drop = FALSE]
}

load_guide_reference <- function(guide_reference_path) {
	guide_reference <- data.table::fread(guide_reference_path)
	unique(guide_reference[, .(
		guide_id = id,
		guide_symbol = name,
		guide_sequence = sequence,
		feature_type
	)])
}

load_guide_element_map <- function(table_s3_path) {
	table_s3 <- data.table::fread(
		table_s3_path,
		select = c("intended_target_name_hg38", "guide_ids")
	)

	guide_element_map <- table_s3[
		!is.na(guide_ids) & nzchar(guide_ids),
		.(guide_id = gsub("_", "-", trimws(unlist(strsplit(guide_ids, ",", fixed = TRUE))), fixed = TRUE)),
		by = .(element = intended_target_name_hg38)
	]
	guide_element_map <- unique(guide_element_map[nzchar(guide_id)])

	guide_conflicts <- guide_element_map[, .(n_elements = data.table::uniqueN(element)), by = guide_id]
	guide_conflicts <- guide_conflicts[n_elements > 1]
	if (nrow(guide_conflicts) > 0) {
		stop(
			"Some guides map to multiple elements in Table S3, e.g.: ",
			paste(utils::head(guide_conflicts$guide_id, 10), collapse = ", ")
		)
	}

	guide_element_map[, .(element = element[1]), by = guide_id]
}

load_validation_seurat_object <- function(
	h5_dir,
	terra_metadata_path,
	guide_reference_path,
	table_s3_path,
	table_s4_path,
	gtf_path,
	gene_type_cache_path = NULL
) {
	h5_files <- list_filtered_h5_files(h5_dir)
	sample_ids <- names(h5_files)
	terra_metadata <- load_terra_metadata(terra_metadata_path, sample_ids)
	readout_gene_ids <- load_readout_gene_ids(table_s4_path)

	template_counts <- extract_modality_counts(
		h5_files[[1]],
		modality_pattern = "Gene Expression",
		use_names = FALSE
	)
	template_counts_named <- extract_modality_counts(
		h5_files[[1]],
		modality_pattern = "Gene Expression",
		use_names = TRUE
	)

	if (is.null(template_counts) || is.null(template_counts_named)) {
		stop("Could not find a Gene Expression assay in the template H5 file.")
	}

	gene_metadata <- build_feature_metadata(
		template_counts,
		template_counts_named,
		id_name = "gene_id",
		symbol_name = "gene_symbol"
	)
	gene_types <- fetch_gene_types_from_gtf(
		gene_ids = gene_metadata$gene_id,
		gtf_path = gtf_path,
		cache_path = gene_type_cache_path
	)
	gene_metadata$gene_type <- gene_types$gene_type
	gene_metadata$is_protein_coding <- gene_metadata$gene_type %in% "protein_coding"
	gene_metadata$is_readout_gene <- gene_metadata$gene_id %in% readout_gene_ids

	template_grna_counts <- extract_modality_counts(
		h5_files[[1]],
		modality_pattern = "CRISPR Guide Capture",
		use_names = FALSE
	)
	template_grna_counts_named <- extract_modality_counts(
		h5_files[[1]],
		modality_pattern = "CRISPR Guide Capture",
		use_names = TRUE
	)

	guide_metadata <- NULL
	if (!is.null(template_grna_counts) && !is.null(template_grna_counts_named)) {
		guide_metadata <- build_feature_metadata(
			template_grna_counts,
			template_grna_counts_named,
			id_name = "guide_id",
			symbol_name = "guide_symbol"
		)

		guide_element_map <- load_guide_element_map(table_s3_path)
		guide_reference <- load_guide_reference(guide_reference_path)
		guide_reference_lookup <- as.data.frame(guide_reference)
		rownames(guide_reference_lookup) <- guide_reference_lookup$guide_id
		guide_element_lookup <- as.data.frame(guide_element_map)
		rownames(guide_element_lookup) <- guide_element_lookup$guide_id
		guide_metadata$guide_sequence <- guide_reference_lookup[
			guide_metadata$guide_id,
			"guide_sequence"
		]
		guide_metadata$feature_type <- guide_reference_lookup[
			guide_metadata$guide_id,
			"feature_type"
		]
		guide_metadata$element <- guide_element_lookup[
			guide_metadata$guide_id,
			"element"
		]
		guide_metadata$element[is.na(guide_metadata$element)] <- "non-targeting"
	}

	counts_list <- vector("list", length(h5_files))
	names(counts_list) <- sample_ids
	grna_counts_list <- vector("list", length(h5_files))
	names(grna_counts_list) <- sample_ids
	cell_metadata_list <- vector("list", length(h5_files))

	for (index in seq_along(h5_files)) {
		sample_id <- sample_ids[[index]]
		counts <- extract_modality_counts(
			h5_files[[index]],
			modality_pattern = "Gene Expression",
			use_names = FALSE
		)
		grna_counts <- extract_modality_counts(
			h5_files[[index]],
			modality_pattern = "CRISPR Guide Capture",
			use_names = FALSE
		)

		if (!identical(rownames(counts), gene_metadata$gene_id)) {
			stop("Gene ordering differs for sample: ", sample_id)
		}

		rownames(counts) <- normalize_feature_names(rownames(counts))

		original_barcodes <- colnames(counts)
		prefixed_barcodes <- paste(sample_id, original_barcodes, sep = "_")
		colnames(counts) <- prefixed_barcodes

		if (!is.null(guide_metadata)) {
			if (is.null(grna_counts)) {
				stop("Missing CRISPR Guide Capture assay for sample: ", sample_id)
			}

			if (!identical(rownames(grna_counts), guide_metadata$guide_id)) {
				stop("Guide ordering differs for sample: ", sample_id)
			}

			if (!identical(colnames(grna_counts), original_barcodes)) {
				stop("RNA and gRNA barcodes differ for sample: ", sample_id)
			}

			rownames(grna_counts) <- normalize_feature_names(rownames(grna_counts))
			colnames(grna_counts) <- prefixed_barcodes
			grna_counts_list[[index]] <- grna_counts
		}

		counts_list[[index]] <- counts
		cell_metadata_list[[index]] <- data.table::data.table(
			cell_barcode = prefixed_barcodes,
			barcode = original_barcodes,
			sample_id = sample_id
		)
	}

	combined_counts <- do.call(cbind, counts_list)

	cell_metadata <- data.table::rbindlist(cell_metadata_list)
	cell_metadata <- merge(
		cell_metadata,
		terra_metadata,
		by = "sample_id",
		all.x = TRUE,
		sort = FALSE
	)

	cell_metadata <- as.data.frame(cell_metadata)
	rownames(cell_metadata) <- cell_metadata$cell_barcode
	cell_metadata <- cell_metadata[colnames(combined_counts), , drop = FALSE]
	cell_metadata$cell_barcode <- NULL

	combined_cdna <- Seurat::CreateSeuratObject(
		counts = combined_counts,
		assay = "RNA",
		meta.data = cell_metadata,
		project = "dc_tap_validation_cdna"
	)

	combined_cdna[["RNA"]] <- SeuratObject::AddMetaData(
		object = combined_cdna[["RNA"]],
		metadata = prepare_feature_metadata(
			feature_metadata = gene_metadata,
			feature_names = rownames(combined_cdna[["RNA"]])
		)
	)

	if (!is.null(guide_metadata)) {
		combined_grna_counts <- do.call(cbind, grna_counts_list)
		combined_cdna[["gRNA"]] <- Seurat::CreateAssayObject(counts = combined_grna_counts)
		combined_cdna[["gRNA"]] <- SeuratObject::AddMetaData(
			object = combined_cdna[["gRNA"]],
			metadata = prepare_feature_metadata(
				feature_metadata = guide_metadata,
				feature_names = rownames(combined_cdna[["gRNA"]])
			)
		)
	}

	combined_cdna
}

get_sample_colors <- function(groups, palette = "Dark 3") {
	group_levels <- unique(as.character(groups))
	group_levels <- group_levels[!is.na(group_levels) & nzchar(group_levels)]

	if (length(group_levels) == 0) {
		return(setNames(character(), character()))
	}

	setNames(grDevices::hcl.colors(length(group_levels), palette = palette), group_levels)
}

plot_knee_by_sample <- function(
	seu,
	assay = "RNA",
	sample_col = "sample",
	label_col = "sample_label_long"
) {
	stopifnot(sample_col %in% colnames(seu[[]]), label_col %in% colnames(seu[[]]))

	counts <- Seurat::GetAssayData(seu, assay = assay, layer = "counts")
	plot_df <- data.table::data.table(
		barcode = colnames(counts),
		total_umis = Matrix::colSums(counts),
		sample = seu[[sample_col]][, 1],
		sample_label = seu[[label_col]][, 1]
	)
	plot_df <- plot_df[!is.na(sample) & !is.na(sample_label) & total_umis > 0]

	if (nrow(plot_df) == 0) {
		stop("No cells with positive counts are available for the knee plot.")
	}

	sample_order <- unique(plot_df$sample_label)
	plot_df[, barcode_rank := seq_len(.N), by = sample]
	data.table::setorder(plot_df, sample, -total_umis)
	plot_df[, barcode_rank := seq_len(.N), by = sample]
	plot_df[, sample_label := factor(sample_label, levels = sample_order)]

	sample_colors <- get_sample_colors(sample_order)

	ggplot2::ggplot(
		plot_df,
		ggplot2::aes(x = barcode_rank, y = total_umis, color = sample_label)
	) +
		ggplot2::geom_line(linewidth = 0.6, alpha = 0.9) +
		ggplot2::facet_grid(rows = ggplot2::vars(sample_label), scales = "free_y") +
		ggplot2::scale_x_log10() +
		ggplot2::scale_y_log10() +
		ggplot2::scale_color_manual(values = sample_colors, guide = "none") +
		ggplot2::labs(
			title = "Knee Plot by Sample",
			x = "Barcode Rank",
			y = "Total RNA UMIs"
		) +
		ggplot2::theme_minimal(base_size = 14) +
		ggplot2::theme(
			strip.text.y = ggplot2::element_text(angle = 0),
			panel.grid.minor = ggplot2::element_blank()
		)
}

# Plot histogram for a value with one median marker per facet.
plot_hist_with_group_medians <- function(df, value_col, group_col = "day_factor",
	title, xlab, bins = 50) {
	stopifnot(value_col %in% names(df), group_col %in% names(df))

	plot_df <- data.frame(
		value = df[[value_col]],
		group = df[[group_col]],
		stringsAsFactors = FALSE
	)
	plot_df <- plot_df[!is.na(plot_df$value) & !is.na(plot_df$group), , drop = FALSE]

	if (nrow(plot_df) == 0) {
		stop("No non-missing values available for plotting.")
	}

	plot_df$group <- factor(plot_df$group, levels = unique(plot_df$group))

	median_stats <- stats::aggregate(
		value ~ group,
		data = plot_df,
		FUN = median
	)
	median_stats$label <- paste0("Median: ", round(median_stats$value, 2))

	ggplot2::ggplot(plot_df, ggplot2::aes(x = value)) +
		ggplot2::geom_histogram(
			bins = bins,
			color = "black",
			fill = "grey70",
			alpha = 0.8
		) +
		ggplot2::geom_vline(
			data = median_stats,
			ggplot2::aes(xintercept = value),
			linetype = "dashed",
			color = "red",
			linewidth = 0.8,
			inherit.aes = FALSE
		) +
		ggplot2::geom_text(
			data = median_stats,
			ggplot2::aes(x = value, y = Inf, label = label),
			vjust = 1.2,
			hjust = -0.1,
			color = "black",
			size = 3.5,
			inherit.aes = FALSE
		) +
		ggplot2::labs(title = title, x = xlab, y = "Barcode Count") +
		ggplot2::facet_grid(rows = ggplot2::vars(group), scales = "free_y") +
		ggplot2::theme_minimal(base_size = 14) +
		ggplot2::theme(
			strip.text.y = ggplot2::element_text(angle = 0),
			panel.grid.minor = ggplot2::element_blank()
		)
}

compute_feature_percentage <- function(seu, features, assay = "RNA", layer = "counts") {
	assay_counts <- Seurat::GetAssayData(seu, assay = assay, layer = layer)
	features <- intersect(features, rownames(assay_counts))

	if (length(features) == 0) {
		stop("No requested features were found in the assay.")
	}

	feature_totals <- Matrix::colSums(assay_counts[features, , drop = FALSE])
	cell_totals <- Matrix::colSums(assay_counts)

	as.numeric((feature_totals / cell_totals) * 100)
}
