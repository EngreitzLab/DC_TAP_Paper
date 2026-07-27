## =================================================================================================
## Confidence intervals for the per-pair "probability of direct effect"
##
## Addresses reviewer point (2): the published "probability of direct effect"
## framework reports a per-pair point estimate but no uncertainty. This script
## adds confidence intervals via a NON-PARAMETRIC CLUSTER (block) BOOTSTRAP that
## resamples perturbations (not by re-running sceptre/MAST).
##
## It faithfully re-implements the published estimation chain
## (see the upstream workflow CRISPR_indirect_effects_FDR_check):
##   1. indirect rate (per dataset) = fraction of significant trans pairs
##   2. cis rate (per distance-to-TSS bin) = fraction of significant cis pairs
##   3. direct rate = max(cis_rate - indirect_rate, 0), averaged across datasets,
##      then a power law direct_rate = a * dist^b is fit by lm(log ~ log)
##   4. per-pair directness = direct_rate(dist) / (direct_rate(dist) + indirect_rate)
##
## The bootstrap resamples each dataset's perturbations with replacement (a
## perturbation's cis and trans pairs move together), recomputes 1-4, and records
## the power-law coefficients and indirect rates for every replicate. Percentiles
## of the resulting directness distribution give the confidence intervals.
##
## Run on a compute node (NOT the login node), e.g.:
##   Rscript bootstrap_directness_CIs.R
## with the `sceptre` conda env (R 4.2.0 + tidyverse). See README.md.
## =================================================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

## Define inputs and parameters --------------------------------------------------------------------

# upstream source workflows that hold the differential-expression (DE) results
SRC <- "/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/CRISPR_indirect_effects_FDR_check"
ENC <- "/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/ENCODE_CRISPR_data"

# shared inputs
GENE_UNIV_FILE <- "/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/CRISPR_benchmarks/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed"
DCTAP_FILE     <- file.path(SRC, "results/annotated_crispr_data/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_fdr20_direct_effects.tsv")

# full direct effects models and indirect effect rates from DC-TAP-seq paper, used to validate that
# bootstrapping results match previous results
PUB_MODEL_FILE <- file.path(SRC, "results/direct_effect_models/direct_effects_model.rds")
PUB_TRANS_FILE <- file.path(SRC, "results/trans_positive_hit_rates.tsv")

# the 6 "held-out" datasets processed by the sceptre pipeline
SCEPTRE_DATASETS <- c("K562_DC_TAPseq", "WTC11_DC_TAPseq", "Klann2021",
                      "Xie2019", "Morris2023_v1", "Morris2023_v2")

# the 1 "training" dataset (Gasperini et al., 2019) processed by the (separate) MAST pipeline
GAS <- list(
  dataset = "Gasperini2019",
  cis     = file.path(ENC, "results/Gasperini2019/diff_expr/output_MAST_perCRE.tsv.gz"),
  trans   = file.path(ENC, "results/Gasperini2019/trans_effects/output_trans_effects_MAST_perCRE.tsv.gz"),
  encode  = file.path(ENC, "results/ENCODE/ENCODE_Gasperini2019_0.13gStd_MAST_perCRE_GRCh38.tsv.gz")
)

# DC-TAP cell type -> dataset whose indirect rate applies (mirrors dctap script)
CELLTYPE_TO_DATASET <- c(K562 = "K562_DC_TAPseq", WTC11 = "WTC11_DC_TAPseq")

# parameters (bin_size / max_distance match the upstream workflow)
B           <- as.integer(Sys.getenv("BOOT_B", "1000"))   # bootstrap replicates
SEED        <- as.integer(Sys.getenv("BOOT_SEED", "20250716"))
BIN_SIZE    <- 5e4
MAX_DIST    <- 1e6
CI_LEVEL    <- 0.95
CORES       <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
DIST_COL    <- "distance_to_gencode_gene_TSS"   # DC-TAP distance column

# effect types: "all" (=significant), "positive", "negative". DC-TAP directness
# is reported for positive & negative (matches predict_direct_vs_indirect_effect_dctap.R)
RATE_TYPES  <- c("all", "positive", "negative")
DCTAP_TYPES <- c("positive", "negative")

# create output directory if needed
OUTDIR <- file.path(getwd(), "results")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# lower, median and upper boundaries of confidence intervals
ci_probs <- c((1 - CI_LEVEL) / 2, 0.5, 1 - (1 - CI_LEVEL) / 2)   # lower, median, upper

set.seed(SEED)
message("Bootstrap CIs for per-pair probability of direct effects | B = ", B, ", seed = ", SEED, 
        ", cores = ", CORES)

## Section 1: loaders -> one normalized per-pair table ---------------------------------------------
## Normalized columns: dataset, perturbation, arm ("cis"/"trans"),
##   significant, regulated_positive, regulated_negative, dist_to_tss (cis only)
## Significance labels are computed once (frozen); the bootstrap resamples
## perturbations, not the significance threshold.

# gene universe (Ensembl IDs) used to filter the sceptre datasets
gene_univ <- read_tsv(GENE_UNIV_FILE, show_col_types = FALSE)

# function to load annotated cis-analysis results and sceptre trans-testing results for one dataset 
load_sceptre_dataset <- function(input_dir, dataset) {
  cis_file   <- file.path(input_dir, "results", dataset, "annotated_cis_results.tsv.gz")
  trans_file <- file.path(input_dir, "results", dataset, "results_run_discovery_analysis.rds")

  # NOTE: annotated cis file is comma-separated despite the .tsv.gz extension
  cis   <- read_csv(cis_file, show_col_types = FALSE)
  trans <- as.data.frame(readRDS(trans_file))

  # cis significance is precomputed (20% FDR); derive the p-value threshold and
  # apply it to trans, exactly as calculate_positive_hit_rates.R does
  pval_threshold_cis <- max(cis$p_value[cis$significant == TRUE], na.rm = TRUE)

  trans <- trans %>%
    mutate(significant        = p_value <= pval_threshold_cis,
           regulated_negative = significant & log_2_fold_change < 0,
           regulated_positive = significant & log_2_fold_change >= 0)
  cis <- cis %>%
    mutate(regulated_negative = significant & log_2_fold_change < 0,
           regulated_positive = significant & log_2_fold_change >= 0)

  # valid distal elements (from cis) define which trans perturbations are valid
  valid_elements <- unique(cis$grna_target[cis$valid_element == TRUE])
  trans <- mutate(trans, valid_connection = grna_target %in% valid_elements)

  # filter cis and trans results based on whether they are valid connection for the analysis
  cis <- filter(cis, valid_connection == TRUE, pass_qc == TRUE,
                response_id %in% gene_univ$Ensembl_ID)
  trans <- filter(trans, valid_connection == TRUE, pass_qc == TRUE,
                  response_id %in% gene_univ$Ensembl_ID)

  # create tables with relevant data for bootstrapping analysis
  cis_norm <- tibble(dataset = dataset, perturbation = cis$grna_target, arm = "cis",
                     significant = cis$significant,
                     regulated_positive = cis$regulated_positive,
                     regulated_negative = cis$regulated_negative,
                     dist_to_tss = abs(cis$dist_to_tss))
  trans_norm <- tibble(dataset = dataset, perturbation = trans$grna_target, arm = "trans",
                       significant = trans$significant,
                       regulated_positive = trans$regulated_positive,
                       regulated_negative = trans$regulated_negative,
                       dist_to_tss = NA_real_)
  output <- bind_rows(cis_norm, trans_norm)
  
  return(output)
  
}

# function to load cis and trans results for Gasperini et al., 2019 dataset
load_gasperini_dataset <- function(g) {
  cis    <- read_tsv(g$cis, show_col_types = FALSE)
  trans  <- read_tsv(g$trans, show_col_types = FALSE)  # already contains re-calculated significance
  encode <- read_tsv(g$encode, show_col_types = FALSE)  # contains ValidConnection column

  # get gene symbols for all genes
  genes <- encode %>% 
    select(gene = measuredEnsemblID, gene_symbol = measuredGeneSymbol) %>%
    distinct()
  
  # add gene symbols and ValidConnection lookups from the ENCODE-format table to cis results
  cis <- cis %>%
    left_join(genes, by = "gene") %>%
    mutate(name = paste0(gene_symbol, "|", pert_chr, ":", pert_start, "-", pert_end, ":.")) %>%
    left_join(distinct(select(encode, name, ValidConnection)), by = "name")

  # extract list of valid candidate enhancers not overlapping TSSs based on ValidConnection column
  enh_filter <- c("overlaps potential promoter", "TSS targeting guide(s)")
  valid_enh  <- encode %>%
    filter(!ValidConnection %in% enh_filter) %>%
    pull(PerturbationTargetID)

  # label pairs with invalid candidate enhancers
  trans <- trans %>%
    mutate(pert_id = paste0(pert_chr, ":", pert_start, "-", pert_end, ":."),
           ValidConnection = if_else(pert_id %in% valid_enh, "TRUE", "Invalid enhancer"))

  # cis significance = adjusted p < 0.05 (MAST pipeline); drop NA effect rows
  cis <- cis %>%
    filter(ValidConnection == "TRUE") %>%
    mutate(significant        = pval_adj < 0.05,
           regulated_negative = significant & logFC < 0,
           regulated_positive = significant & logFC >= 0) %>%
    filter(!is.na(regulated_negative), !is.na(regulated_positive))
  trans <- filter(trans, ValidConnection == "TRUE")

  # create tables with relevant data for bootstrapping analysis
  cis_norm <- tibble(dataset = g$dataset, perturbation = cis$perturbation, arm = "cis",
                     significant = cis$significant,
                     regulated_positive = cis$regulated_positive,
                     regulated_negative = cis$regulated_negative,
                     dist_to_tss = abs(cis$dist_to_tss))
  trans_norm <- tibble(dataset = g$dataset, perturbation = trans$perturbation, arm = "trans",
                       significant = trans$significant,
                       regulated_positive = trans$regulated_positive,
                       regulated_negative = trans$regulated_negative,
                       dist_to_tss = NA_real_)
  output <- bind_rows(cis_norm, trans_norm)
  
  return(output)
}

message("Loading and pre-processing DE results for all datasets ...")
norm_list <- c(
  lapply(SCEPTRE_DATASETS, load_sceptre_dataset, input_dir = SRC),
  list(load_gasperini_dataset(GAS))
)
norm <- bind_rows(norm_list)
ALL_DATASETS <- c(SCEPTRE_DATASETS, GAS$dataset)

## Export cis support (pairs + perturbations per distance) for the report --------------------------
## `norm` already holds, per cis pair: dataset, perturbation, distance, significance and effect
## direction. This is the observational support behind the direct-effect rate. Since the bootstrap
## resamples the `perturbation` clusters, both the pair and perturbation counts are informative.
cis_support <- norm %>%
  filter(arm == "cis") %>%
  transmute(dataset, perturbation, dist_to_tss,
            significant, regulated_positive, regulated_negative)
write_tsv(cis_support, file.path(OUTDIR, "cis_support_by_distance.tsv"))

# regenerate just the support table (no bootstrap) with: SUPPORT_ONLY=1 Rscript bootstrap_directness_CIs.R
if (nzchar(Sys.getenv("SUPPORT_ONLY"))) {
  message("SUPPORT_ONLY set -- wrote cis_support_by_distance.tsv and exiting before the bootstrap.")
  quit(save = "no")
}

## Section 2: precompute per-dataset bootstrap structures ------------------------------------------
## For each dataset: perturbation index (shared by cis+trans), fixed distance
## bins (from full cis data), and the logical hit vectors. The bootstrap then
## only reweights rows by per-perturbation draw counts (cluster bootstrap).

# coerce a hit indicator to numeric with NA -> 0. This reproduces the upstream
# rate scripts, which count hits with sum(..., na.rm = TRUE) while keeping the
# full pair count n() in the denominator (NA = not a hit, still a tested pair).
na0 <- function(x) { x <- as.numeric(x); x[is.na(x)] <- 0; x }

build_dataset_struct <- function(d) {
  sub   <- filter(norm, dataset == d)
  cis   <- filter(sub, arm == "cis")
  trans <- filter(sub, arm == "trans")

  perts <- unique(c(cis$perturbation, trans$perturbation))  # should be the same for cis and trans

  # fixed distance bins for this dataset (from full cis data), reused every rep
  max_dist <- ceiling(max(cis$dist_to_tss, na.rm = TRUE) / 1e6) * 1e6
  breaks   <- seq(0, max_dist, by = BIN_SIZE)
  cis_bin  <- cut(cis$dist_to_tss, breaks = breaks, include.lowest = TRUE)

  # construct object containing relevant information from the given dataset
  output <- list(
    dataset = d,
    n_pert  = length(perts),
    cis = list(pert_idx = match(cis$perturbation, perts),
               bin      = cis_bin,
               dist     = cis$dist_to_tss,
               sig      = na0(cis$significant),
               pos      = na0(cis$regulated_positive),
               neg      = na0(cis$regulated_negative)),
    trans = list(pert_idx = match(trans$perturbation, perts),
                 sig = na0(trans$significant),
                 pos = na0(trans$regulated_positive),
                 neg = na0(trans$regulated_negative))
  )
  
  return(output)
  
}

ds_structs <- lapply(ALL_DATASETS, build_dataset_struct)
names(ds_structs) <- ALL_DATASETS

## Section 2b: rate + power-law functions ----------------------------------------------------------

# per-dataset direct rates (per bin) + trans/indirect rates, given perturbation
# draw counts (weights). weights = rep(1, n_pert) reproduces the point estimate.
dataset_rates <- function(ds, weights) {
  
  # calculate how many times each E-G pair occurs in the bootstrap sample based on weights
  wc <- weights[ds$cis$pert_idx]
  wt <- weights[ds$trans$pert_idx]

  # cis: weighted sums per distance bin
  M <- rowsum(cbind(w = wc, sig = wc * ds$cis$sig, pos = wc * ds$cis$pos,
                    neg = wc * ds$cis$neg, wd = wc * ds$cis$dist),
              group = ds$cis$bin, reorder = FALSE)
  tot <- M[, "w"]
  keep <- tot > 0
  M <- M[keep, , drop = FALSE]; tot <- tot[keep]

  # cis hit rate across distance bins per rate type
  cis_rate <- data.frame(
    dataset   = ds$dataset,
    dist_bin  = rownames(M),
    mean_dist = M[, "wd"] / tot,
    cis_all   = M[, "sig"] / tot,
    cis_pos   = M[, "pos"] / tot,
    cis_neg   = M[, "neg"] / tot,
    stringsAsFactors = FALSE, row.names = NULL
  )

  # distance independent trans hit rate per rate type
  tt <- sum(wt)
  ind <- c(all = sum(wt * ds$trans$sig) / tt,
           positive = sum(wt * ds$trans$pos) / tt,
           negative = sum(wt * ds$trans$neg) / tt)

  # direct rate = max(cis - indirect, 0), per type
  cis_rate$direct_all <- pmax(cis_rate$cis_all - ind["all"], 0)
  cis_rate$direct_pos <- pmax(cis_rate$cis_pos - ind["positive"], 0)
  cis_rate$direct_neg <- pmax(cis_rate$cis_neg - ind["negative"], 0)

  # create output containing boootstrapped direct and indirect hit rates
  output <- list(cis = cis_rate, indirect = ind)
  
  return(output)
  
}

# fit the three power-law models from all datasets' per-bin direct rates;
# returns intercept (= log a) and slope (= b) for each type
fit_models <- function(rates_list) {
  
  # get cis hit rates within given maximum distance
  cis <- do.call(rbind, lapply(rates_list, `[[`, "cis"))
  cis <- cis[cis$mean_dist <= MAX_DIST, , drop = FALSE]

  # average across datasets within each distance bin
  avg <- cis %>%
    group_by(dist_bin) %>%
    summarize(dist_to_tss = mean(mean_dist),
              direct_all  = mean(direct_all),
              direct_pos  = mean(direct_pos),
              direct_neg  = mean(direct_neg),
              .groups = "drop")

  # fit power law model to bootstrapped direct effects rates
  fit_one <- function(y) {
    df <- data.frame(direct_rate = avg[[y]], dist_to_tss = avg$dist_to_tss)
    df <- df[df$direct_rate > 0, , drop = FALSE]
    co <- coef(lm(log(direct_rate) ~ log(dist_to_tss), data = df))
    c(intercept = unname(co[1]), slope = unname(co[2]))
  }
  
  # create output
  output <- list(all      = fit_one("direct_all"),
                 positive = fit_one("direct_pos"),
                 negative = fit_one("direct_neg"))
  
  return(output)
  
}

# one bootstrap replicate: resample each dataset's perturbations with replacement
one_replicate <- function(seed_i) {
  set.seed(seed_i)
  rates_list <- lapply(ds_structs, function(ds) {
    counts <- tabulate(sample.int(ds$n_pert, ds$n_pert, replace = TRUE), nbins = ds$n_pert)
    dataset_rates(ds, counts)
  })
  models   <- fit_models(rates_list)
  indirect <- sapply(rates_list, `[[`, "indirect")           # types x datasets
  colnames(indirect) <- ALL_DATASETS
  list(models = models, indirect = indirect)
}

## Section 3: point estimate + validation ----------------------------------------------------------

point_rates    <- lapply(ds_structs, function(ds) dataset_rates(ds, rep(1, ds$n_pert)))
point_models   <- fit_models(point_rates)
point_indirect <- sapply(point_rates, `[[`, "indirect"); colnames(point_indirect) <- ALL_DATASETS

# validate power-law coefficients against the published model
pub_model <- readRDS(PUB_MODEL_FILE)
pub_coef  <- sapply(c(all = "direct_rate_all", positive = "direct_rate_positive",
                      negative = "direct_rate_negative"),
                    function(nm) coef(pub_model[[nm]]))
message("\n== Validation: recomputed vs published power-law coefficients ==")
for (ty in RATE_TYPES) {
  rc <- point_models[[ty]]
  pc <- pub_coef[, ty]
  message(sprintf("  %-9s intercept %.5f (pub %.5f)  slope %.5f (pub %.5f)",
                  ty, rc["intercept"], pc[1], rc["slope"], pc[2]))
  stopifnot(abs(rc["intercept"] - pc[1]) < 1e-3, abs(rc["slope"] - pc[2]) < 1e-3)
}

# validate indirect (trans) rates against the published table
pub_trans <- read_tsv(PUB_TRANS_FILE, show_col_types = FALSE)
message("== Validation: recomputed vs published indirect rates (significant) ==")
for (d in ALL_DATASETS) {
  pubv <- pub_trans$positive_rate_significant[pub_trans$dataset == d]
  recv <- point_indirect["all", d]
  message(sprintf("  %-16s recomputed %.6f  published %.6f", d, recv, pubv))
  stopifnot(abs(recv - pubv) < 1e-4)
}
message("Validation passed: re-implementation reproduces the published point estimate.\n")

## Section 4: bootstrap loop -----------------------------------------------------------------------

message("Running ", B, " bootstrap replicates ...")
seeds <- SEED + seq_len(B)
boot <- if (CORES > 1) {
  parallel::mclapply(seeds, one_replicate, mc.cores = CORES)
} else {
  lapply(seeds, one_replicate)
}

# tidy the bootstrap draws
coef_rows <- lapply(seq_along(boot), function(i) {
  m <- boot[[i]]$models
  do.call(rbind, lapply(RATE_TYPES, function(ty) {
    data.frame(replicate = i, type = ty,
               intercept = m[[ty]]["intercept"], slope = m[[ty]]["slope"],
               a = exp(m[[ty]]["intercept"]), b = m[[ty]]["slope"],
               row.names = NULL)
    }))
})
coef_df <- do.call(rbind, coef_rows)

indirect_rows <- lapply(seq_along(boot), function(i) {
  ind <- boot[[i]]$indirect
  do.call(rbind, lapply(ALL_DATASETS, function(d)
    data.frame(replicate = i, dataset = d, type = RATE_TYPES,
               indirect_rate = ind[RATE_TYPES, d], row.names = NULL)))
})
indirect_df <- do.call(rbind, indirect_rows)

# add the point estimate as replicate 0 for reference
point_coef <- do.call(rbind, lapply(RATE_TYPES, function(ty)
  data.frame(replicate = 0L, type = ty,
             intercept = point_models[[ty]]["intercept"], slope = point_models[[ty]]["slope"],
             a = exp(point_models[[ty]]["intercept"]), b = point_models[[ty]]["slope"],
             row.names = NULL)))
coef_df <- rbind(point_coef, coef_df)

write_tsv(coef_df, file.path(OUTDIR, "bootstrap_powerlaw_coefficients.tsv"))
write_tsv(indirect_df, file.path(OUTDIR, "bootstrap_indirect_rates.tsv"))

## Section 5: propagate to directness CIs ----------------------------------------------------------

# helper: directness probability across all B replicates for a vector of
# distances, one cell type, one effect type -> returns a B x length(dist) matrix
directness_matrix <- function(dist, type, dataset) {
  inter <- coef_df$intercept[coef_df$replicate >= 1 & coef_df$type == type]
  slope <- coef_df$slope[coef_df$replicate >= 1 & coef_df$type == type]
  ind   <- indirect_df$indirect_rate[indirect_df$dataset == dataset & indirect_df$type == type]
  logd  <- log(dist)
  # B x n: direct rate (capped at 1), then probability
  direct <- pmin(exp(outer(inter, rep(1, length(dist))) + outer(slope, logd)), 1)
  direct / (direct + ind)          # ind (length B) recycles down each column
}

# point directness for a distance vector (published coefficients + published indirect)
point_directness <- function(dist, type, dataset) {
  co <- point_models[[type]]
  ind <- point_indirect[type, dataset]
  direct <- pmin(exp(co["intercept"] + co["slope"] * log(dist)), 1)
  unname(direct / (direct + ind))
}

## 5a. smooth CI curves over a distance grid (for the ribbon plot)
grid <- 10 ^ seq(log10(1e3), log10(MAX_DIST), length.out = 200)
curve_rows <- list(); k <- 1
for (ct in names(CELLTYPE_TO_DATASET)) {
  d <- CELLTYPE_TO_DATASET[[ct]]
  for (ty in DCTAP_TYPES) {
    mat <- directness_matrix(grid, ty, d)
    qs  <- apply(mat, 2, quantile, probs = ci_probs, na.rm = TRUE)
    curve_rows[[k]] <- data.frame(
      dist_to_tss = grid, cell_type = ct, type = ty,
      point  = point_directness(grid, ty, d),
      median = qs[2, ], lower = qs[1, ], upper = qs[3, ], row.names = NULL)
    k <- k + 1
  }
}
curve_df <- do.call(rbind, curve_rows)
write_tsv(curve_df, file.path(OUTDIR, "directness_probability_CI_by_distance.tsv"))

## 5b. per-pair CIs for the DC-TAP results file
message("Computing per-pair directness CIs for the DC-TAP results ...")
dctap <- read_tsv(DCTAP_FILE, show_col_types = FALSE)
dctap$.dist <- pmax(abs(dctap[[DIST_COL]]), 1)      # guard log(0)/negatives
dctap$.dataset <- unname(CELLTYPE_TO_DATASET[dctap$cell_type])

for (ty in DCTAP_TYPES) {
  lo <- md <- up <- rep(NA_real_, nrow(dctap))
  for (ct in names(CELLTYPE_TO_DATASET)) {
    idx <- which(dctap$cell_type == ct)
    if (!length(idx)) next
    mat <- directness_matrix(dctap$.dist[idx], ty, CELLTYPE_TO_DATASET[[ct]])
    qs  <- apply(mat, 2, quantile, probs = ci_probs, na.rm = TRUE)
    lo[idx] <- qs[1, ]; md[idx] <- qs[2, ]; up[idx] <- qs[3, ]
  }
  dctap[[paste0("direct_vs_indirect_", ty, "_median")]] <- md
  dctap[[paste0("direct_vs_indirect_", ty, "_lower")]]  <- lo
  dctap[[paste0("direct_vs_indirect_", ty, "_upper")]]  <- up
}
dctap <- select(dctap, -c(.dist, .dataset))

out_dctap <- file.path(OUTDIR, sub("\\.tsv$", "_with_CIs.tsv", basename(DCTAP_FILE)))
write_tsv(dctap, out_dctap)

## Save a compact object for the report + session info ---------------------------------------------
saveRDS(list(coef = coef_df, indirect = indirect_df, curves = curve_df,
             point_models = point_models, point_indirect = point_indirect,
             pub_coef = pub_coef, B = B, seed = SEED, ci_level = CI_LEVEL,
             celltype_to_dataset = CELLTYPE_TO_DATASET, dctap_types = DCTAP_TYPES,
             dctap_ci_file = out_dctap),
        file.path(OUTDIR, "bootstrap_summary.rds"))

message("\nDone. Outputs written to ", OUTDIR)
message("  - bootstrap_powerlaw_coefficients.tsv")
message("  - bootstrap_indirect_rates.tsv")
message("  - directness_probability_CI_by_distance.tsv")
message("  - ", basename(out_dctap))
message("  - bootstrap_summary.rds")
print(sessionInfo())
