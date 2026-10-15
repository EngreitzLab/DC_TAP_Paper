#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)])

if (length(script_path) != 1 || !nzchar(script_path)) {
	stop("Could not determine script path.")
}

repo_root <- normalizePath(
	file.path(dirname(script_path), "..", ".."),
	winslash = "/",
	mustWork = TRUE
)

terra_metadata_path <- file.path(repo_root, "data", "terra_table_fixed_pairing.tsv")
output_dir <- file.path(repo_root, "results", "cellranger")
output_path <- file.path(output_dir, "cellranger_metrics_summary.tsv")
plot_dir <- file.path(output_dir, "pngs")

system_bins <- c("/usr/bin", "/bin", "/usr/sbin", "/sbin")
current_path <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
Sys.setenv(PATH = paste(unique(c(system_bins, current_path)), collapse = .Platform$path.sep))

gsutil_candidates <- unique(c(
	Sys.which("gsutil"),
	path.expand("~/Programs/google-cloud-sdk/bin/gsutil"),
	"/opt/homebrew/bin/gsutil",
	"/usr/local/bin/gsutil"
))
gsutil_candidates <- gsutil_candidates[nzchar(gsutil_candidates)]
gsutil_path <- gsutil_candidates[file.exists(gsutil_candidates)][1]

if (is.na(gsutil_path) || !nzchar(gsutil_path)) {
	stop("Could not find gsutil. Install Google Cloud SDK or add gsutil to PATH.")
}

field_map <- c(
	cells_recovered = "Estimated Number of Cells",
	mRNA_reads = "Number of Reads",
	mean_mRNA_reads_per_cell = "Mean Reads per Cell",
	Reads_mapped_to_targeted_transcriptome = "Reads Mapped Confidently to Transcriptome",
	Median_targeted_mRNA_UMI_per_cell = "Median UMI Counts per Cell",
	mRNA_seq_saturation = "Sequencing Saturation",
	guide_reads = "CRISPR: Number of Reads",
	percent_guide_reads_with_protospacer = "CRISPR: Fraction Reads with Putative Protospacer Sequence",
	median_guide_UMIs_per_cell = "CRISPR: Median UMIs per Cell (summed over all recognized protospacers)",
	guide_sequencing_saturation = "CRISPR: Sequencing Saturation"
)

parse_metric_value <- function(value) {
	value <- trimws(as.character(value))

	if (!nzchar(value) || value %in% c("NA", "NaN", "NAN")) {
		return(NA_real_)
	}

	value <- gsub(",", "", value, fixed = TRUE)
	value <- sub("%$", "", value)
	suppressWarnings(as.numeric(value))
}

format_metric_label <- function(metric_name) {
	gsub("_", " ", metric_name, fixed = TRUE)
}

write_metric_barplot <- function(metrics_summary, metric_name, sample_colors) {
	plot_data <- metrics_summary[, c("sample", "sample_label_long", metric_name)]
	names(plot_data)[3] <- "metric_value"

	png(
		filename = file.path(plot_dir, paste0(metric_name, ".png")),
		width = 2200,
		height = 1400,
		res = 200
	)

	op <- par(mar = c(5, 18, 5, 12), xpd = FALSE)
	on.exit({
		par(op)
		dev.off()
	}, add = TRUE)

	barplot(
		height = plot_data$metric_value,
		names.arg = plot_data$sample_label_long,
		horiz = TRUE,
		las = 1,
		col = unname(sample_colors[plot_data$sample]),
		border = NA,
		main = format_metric_label(metric_name),
		xlab = format_metric_label(metric_name),
		cex.names = 0.9
	)

	legend(
		"topright",
		legend = names(sample_colors),
		fill = unname(sample_colors),
		title = "sample",
		bty = "n",
		inset = c(-0.18, 0),
		xpd = NA
	)
}

read_metrics_summary <- function(metrics_uri) {
	metrics_text <- system2(gsutil_path, c("cat", metrics_uri), stdout = TRUE, stderr = TRUE)
	status <- attr(metrics_text, "status")

	if (!is.null(status) && status != 0) {
		stop(
			"Failed to read metrics from ",
			metrics_uri,
			"\n",
			paste(metrics_text, collapse = "\n")
		)
	}

	read.csv(
		text = paste(metrics_text, collapse = "\n"),
		check.names = FALSE,
		stringsAsFactors = FALSE
	)
}

extract_sample_metrics <- function(sample_row) {
	metrics_uri <- paste0(sub("/+$", "", sample_row$output_folder), "/outs/metrics_summary.csv")
	metrics <- read_metrics_summary(metrics_uri)

	if (nrow(metrics) != 1) {
		stop("Expected exactly one row in metrics_summary.csv for ", sample_row$sample_id)
	}

	metric_values <- lapply(field_map, function(source_name) {
		if (!source_name %in% names(metrics)) {
			stop("Missing metric '", source_name, "' for sample ", sample_row$sample_id)
		}

		parse_metric_value(metrics[[source_name]][1])
	})

	data.frame(
		sample_id = sample_row$sample_id,
		sample = sample_row$sample,
		sample_label_long = sample_row$sample_label_long,
		metrics_uri = metrics_uri,
		as.data.frame(metric_values, check.names = FALSE),
		check.names = FALSE,
		stringsAsFactors = FALSE
	)
}

terra_metadata <- read.delim(
	terra_metadata_path,
	sep = "\t",
	check.names = FALSE,
	stringsAsFactors = FALSE
)

sample_table <- unique(terra_metadata[, c("entity:validation_id", "sample", "sample_label_long", "output_folder")])
names(sample_table)[names(sample_table) == "entity:validation_id"] <- "sample_id"

metrics_rows <- lapply(seq_len(nrow(sample_table)), function(index) {
	extract_sample_metrics(sample_table[index, , drop = FALSE])
})

metrics_summary <- do.call(rbind, metrics_rows)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
	metrics_summary,
	file = output_path,
	sep = "\t",
	quote = FALSE,
	row.names = FALSE
)

sample_levels <- unique(metrics_summary$sample)
sample_colors <- setNames(
	grDevices::hcl.colors(length(sample_levels), palette = "Dark 3"),
	sample_levels
)

for (metric_name in names(field_map)) {
	write_metric_barplot(metrics_summary, metric_name, sample_colors)
}

message("Wrote Cell Ranger metrics summary to ", output_path)
message("Wrote metric PNGs to ", plot_dir)