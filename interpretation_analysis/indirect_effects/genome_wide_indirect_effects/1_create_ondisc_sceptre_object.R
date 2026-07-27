## Create sceptre on-disc object for Gasperini et al., 2019 data for use with Nextflow pipeline

suppressPackageStartupMessages({
  library(sceptre)
  library(ondisc)
})

# load in-memory sceptre object from main analyses of Gasperini et al dataset
sceptre_object_file <- "/oak/stanford/groups/engreitz/Users/jgalante/DC_TAP_Paper/results/main_figure_1_and_2/duplicate_pairs_analysis/differential_expression/sceptre_diffex_input.rds"
sceptre_object <- readRDS(sceptre_object_file)

# extract used formula and save to formula object file for sceptre Nextflow pipeline
formula <- sceptre_object@formula_object
saveRDS(formula, file = "sceptre_formula.rds")

# get used MOI setting
moi <- ifelse(isTRUE(sceptre_object@low_moi), yes = "low", no = "high")

# extract gene (response) and grna count matrices
response_matrix <- sceptre_object@response_matrix[[1]]
grna_matrix <- sceptre_object@grna_matrix[[1]]

# extract required columns from gRNA targets data frame
grna_target_df <- sceptre_object@grna_target_data_frame[, c("grna_id", "grna_target")]

# extract covariates data frame
covariate_df <- sceptre_object@covariate_data_frame

# remove any covariates computed internally when creating the new sceptre object
computed_covars <- c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis",
                     "response_p_mito")
extra_cols <- setdiff(colnames(covariate_df), computed_covars)
if (length(extra_cols) > 0L) {
  extra_covars <- covariate_df[, extra_cols, drop = FALSE]
} else {
  extra_covars <- data.frame()
}

# create new on-disc sceptre object
output_dir <- "sceptre_object_ondisc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ondisc_sceptre_object <- import_data(
  response_matrix = response_matrix, 
  grna_matrix = grna_matrix,
  grna_target_data_frame = grna_target_df,
  moi = moi,
  extra_covariates = extra_covars,
  use_ondisc = TRUE,
  directory_to_write = output_dir
)

# get side parameter for differential expression testsfrom input sceptre object
side <- switch(as.character(sceptre_object@side_code), "-1" = "left", "0" = "both", "1" = "right",
               stop("Unexpected side_code: ", sceptre_object@side_code))

# set analysis parameters
ondisc_sceptre_object <- set_analysis_parameters(
  sceptre_object = ondisc_sceptre_object,
  side = side,
  formula_object = formula,
  resampling_mechanism = "permutations"
)

# save ondisc sceptre object to file
write_ondisc_backed_sceptre_object(ondisc_sceptre_object, directory_to_write = output_dir)
