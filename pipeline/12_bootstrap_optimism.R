#!/usr/bin/env Rscript
# =============================================================================
# STEP 12 (Signature v2) — Bootstrap Optimism Correction
# -----------------------------------------------------------------------------
# Bootstrap estimate of optimism for the v2 signature
#
# Features:
#   - the test set is not used
#   - the full pipeline is rebuilt in every replicate (screening + scaling + CV + glmnet)
#   - scaling parameters come from the bootstrap sample only
#   - lambda is chosen by inner CV
#   - C_corrected = C_apparent − mean(optimism)
#
# Outputs are saved in signature_v2/.
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(data.table)
  library(jsonlite)
})

output_dir <- WORK_DIR
sig_dir    <- file.path(output_dir, "signature_v2")
dir.create(sig_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

# ======================================================================
# 0) Parameters
# ======================================================================
B                <- 1000
ALPHA            <- 0.5
UNIV_P_STRICT    <- 0.01   # as in step 7 v2
UNIV_P_RELAX     <- 0.05   # fallback
MIN_SCREEN_N     <- 20     # threshold for the fallback
N_CV_FOLDS       <- 5
SEED             <- 2024
MIN_EVENTS       <- 5

normalize_tcga_id <- function(ids) {
  ids <- toupper(as.character(ids))
  ids <- gsub("\\.", "-", ids)
  m <- regmatches(ids, regexpr("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", ids))
  out <- rep(NA_character_, length(ids))
  ok <- nchar(m) > 0
  out[ok] <- m[ok]
  out
}

# ======================================================================
# 1) Rebuild the preprocessed train (as in step 7 v2)
# ======================================================================
log_msg("============================================================")
log_msg("STEP 12 (Signature v2) — Bootstrap Optimism Correction")
log_msg("============================================================")
log_msg("rebuilding the train data ...")

isoform_data  <- read.csv(file.path(output_dir, "mirna_isoform_rpm_matrix.csv"),
                          row.names = 1, check.names = FALSE)
clinical_data <- readRDS(file.path(output_dir, "raw_clinical_cdr_integrated.rds"))
train_df      <- read.csv(file.path(output_dir, "train_set.csv"),
                          stringsAsFactors = FALSE)
train_ids     <- normalize_tcga_id(train_df$patient_barcode)

# primary tumor filter
mirna_raw <- readRDS(file.path(output_dir, "raw_mirna.rds"))
if (!is.data.table(mirna_raw)) mirna_raw <- as.data.table(mirna_raw)
mirna_meta <- unique(mirna_raw[, .(sample_barcode, is_primary_tumor)])
mirna_meta[, key := normalize_tcga_id(sample_barcode)]
primary <- unique(mirna_meta[is_primary_tumor == TRUE, key])
primary <- primary[!is.na(primary)]

iso_ids  <- normalize_tcga_id(rownames(isoform_data))
clin_ids <- normalize_tcga_id(as.character(clinical_data$patient_barcode))

common <- Reduce(intersect, list(train_ids, iso_ids, clin_ids, primary))
log_msg(sprintf("   shared train samples: %d", length(common)))

idx_iso  <- match(common, iso_ids)
idx_clin <- match(common, clin_ids)

iso_train <- isoform_data[idx_iso, , drop = FALSE]
rownames(iso_train) <- common
clin <- as.data.frame(clinical_data[idx_clin, , drop = FALSE])
rownames(clin) <- common

clin$status        <- clin$status_final
clin$survival_time <- clin$survival_time_final

keep <- !is.na(clin$status) & !is.na(clin$survival_time) &
        clin$survival_time >= 0
clin <- clin[keep, , drop = FALSE]
iso_train <- iso_train[rownames(clin), , drop = FALSE]

X_raw <- log2(as.matrix(iso_train) + 1)

keep_rows <- complete.cases(X_raw)
X_raw  <- X_raw[keep_rows, , drop = FALSE]
clin   <- clin[keep_rows, , drop = FALSE]

# expression pre-filter (same as step 7 v2)
n_samples <- nrow(X_raw)
keep_expr <- colSums(X_raw > 1, na.rm = TRUE) >= 0.5 * n_samples
X_log     <- X_raw[, keep_expr, drop = FALSE]
vars      <- apply(X_log, 2, var, na.rm = TRUE)
keep_var  <- vars > quantile(vars, 0.2, na.rm = TRUE)
keep_final <- colnames(X_log)[keep_var]

X_full <- X_raw[, keep_final, drop = FALSE]

# Age
age_raw <- suppressWarnings(as.numeric(clin$age_at_diagnosis))
if (all(is.na(age_raw))) age_raw <- suppressWarnings(as.numeric(clin$age))
age_years <- if (median(age_raw, na.rm = TRUE) > 200) age_raw / 365.25 else age_raw
age_years[is.na(age_years)] <- median(age_years, na.rm = TRUE)

clin_cov_matrix <- matrix(age_years, ncol = 1, dimnames = list(NULL, "Age"))
y_full <- Surv(clin$survival_time, clin$status)

n <- nrow(X_full)
log_msg(sprintf("   n = %d | events = %d | miRNAs = %d | clinical = %d",
                n, sum(clin$status == 1), ncol(X_full), ncol(clin_cov_matrix)))

use_method <- "survConcordance"

# ======================================================================
# 2) Helper functions
# ======================================================================
compute_cindex <- function(y_surv, lp, method = "survConcordance") {
  lp <- as.numeric(lp)
  if (length(lp) == 0 || all(!is.finite(lp)) ||
      length(unique(lp[is.finite(lp)])) <= 1) {
    return(0.5)
  }
  if (method == "survConcordance") {
    return(survival::survConcordance(y_surv ~ lp)$concordance)
  } else {
    return(survival::concordance(y_surv ~ lp, reverse = TRUE)$concordance)
  }
}

fit_scaler <- function(X) {
  center <- colMeans(X, na.rm = TRUE)
  scale_value <- apply(X, 2, sd, na.rm = TRUE)
  scale_value[!is.finite(scale_value) | scale_value <= 0] <- 1
  list(center = center, scale = scale_value)
}

apply_scaler <- function(X, scaler) {
  X <- sweep(X, 2, scaler$center, "-")
  X <- sweep(X, 2, scaler$scale, "/")
  X
}

univariate_screen <- function(X, y) {
  p_values <- rep(NA_real_, ncol(X))
  names(p_values) <- colnames(X)

  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    if (all(!is.finite(x)) || sd(x, na.rm = TRUE) == 0) next
    fit <- tryCatch(survival::coxph(y ~ x), error = function(e) NULL)
    if (!is.null(fit)) {
      ct <- summary(fit)$coefficients
      if (nrow(ct) >= 1 && ncol(ct) >= 5) p_values[j] <- ct[1, 5]
    }
  }

  # two-stage, as in step 7 v2
  selected <- which(!is.na(p_values) & p_values < UNIV_P_STRICT)
  if (length(selected) < MIN_SCREEN_N) {
    selected <- which(!is.na(p_values) & p_values < UNIV_P_RELAX)
  }

  list(p_values = p_values, selected = selected)
}

# ======================================================================
# 3) Fit the full pipeline on one sample
# ======================================================================
fit_bootstrap_pipeline <- function(X, clin, y, alpha = 0.5, nfolds = 5) {
  n_events <- sum(y[, "status"] == 1)
  if (n_events < MIN_EVENTS) {
    return(list(success = FALSE, reason = "Too few events"))
  }

  # 1) univariate screening
  screening <- univariate_screen(X, y)
  selected_idx <- screening$selected
  if (length(selected_idx) < 1) {
    return(list(success = FALSE, reason = "No miRNA passed screening"))
  }

  X_mir <- X[, selected_idx, drop = FALSE]

  # 2) scaling of the miRNAs from the bootstrap
  mir_scaler   <- fit_scaler(X_mir)
  X_mir_scaled <- apply_scaler(X_mir, mir_scaler)

  valid_mir <- apply(X_mir_scaled, 2, function(z) all(is.finite(z)))
  X_mir_scaled <- X_mir_scaled[, valid_mir, drop = FALSE]
  if (ncol(X_mir_scaled) < 1) {
    return(list(success = FALSE, reason = "No valid miRNA after scaling"))
  }

  # 3) scaling of the clinical variable
  clin_scaler   <- fit_scaler(clin)
  X_clin_scaled <- apply_scaler(clin, clin_scaler)

  X_combined <- cbind(X_mir_scaled, X_clin_scaled)
  if (ncol(X_combined) < 2) {
    return(list(success = FALSE, reason = "Insufficient predictors"))
  }

  # 4) event-balanced folds
  folds <- min(nfolds, n_events)
  if (folds < 3) {
    return(list(success = FALSE, reason = "Too few events for CV"))
  }

  event_idx <- which(y[, "status"] == 1)
  cens_idx  <- which(y[, "status"] == 0)
  foldid    <- integer(length(y))
  foldid[event_idx] <- sample(rep(seq_len(folds), length.out = length(event_idx)))
  foldid[cens_idx]  <- sample(rep(seq_len(folds), length.out = length(cens_idx)))

  cv_fit <- tryCatch(
    cv.glmnet(x = X_combined, y = y, family = "cox",
              alpha = alpha, nfolds = folds, foldid = foldid,
              standardize = FALSE),
    error = function(e) NULL
  )
  if (is.null(cv_fit)) {
    return(list(success = FALSE, reason = "cv.glmnet failed"))
  }

  final_fit <- tryCatch(
    glmnet(x = X_combined, y = y, family = "cox", alpha = alpha,
           lambda = cv_fit$lambda.min, standardize = FALSE),
    error = function(e) NULL
  )
  if (is.null(final_fit)) {
    return(list(success = FALSE, reason = "glmnet failed"))
  }

  coef_matrix <- as.matrix(coef(final_fit))
  nonzero <- which(coef_matrix[, 1] != 0)
  nonzero_names <- rownames(coef_matrix)[nonzero]
  nonzero_mirna <- setdiff(nonzero_names, colnames(clin))

  list(
    success = TRUE,
    reason = NA_character_,
    model = final_fit,
    mirna_names = colnames(X_mir_scaled),
    selected_screening_names = colnames(X_mir),
    nonzero_names = nonzero_names,
    nonzero_mirna = nonzero_mirna,
    mir_scaler = mir_scaler,
    clin_scaler = clin_scaler,
    lambda = cv_fit$lambda.min,
    n_screened = ncol(X_mir),
    n_nonzero = length(nonzero),
    X_boot = X_combined
  )
}

# ======================================================================
# 4) Apparent model
# ======================================================================
log_msg("")
log_msg("------------------------------------------------------------")
log_msg("Fitting apparent model")
log_msg("------------------------------------------------------------")

set.seed(SEED)
apparent_fit <- fit_bootstrap_pipeline(X_full, clin_cov_matrix, y_full,
                                        alpha = ALPHA, nfolds = N_CV_FOLDS)

if (!isTRUE(apparent_fit$success)) {
  stop("Apparent fit failed: ", apparent_fit$reason)
}

lp_apparent <- predict(apparent_fit$model,
                       newx = apparent_fit$X_boot, type = "link")
C_apparent <- compute_cindex(y_full, lp_apparent, use_method)

log_msg(sprintf("C-index apparent = %.4f", C_apparent))
log_msg(sprintf("miRNA after screening = %d", apparent_fit$n_screened))
log_msg(sprintf("Non-zero coefficients = %d", apparent_fit$n_nonzero))

# ======================================================================
# 5) Bootstrap (serial — for cross-system compatibility, instead of parallel)
# ======================================================================
log_msg("")
log_msg(sprintf("%d bootstrap replicates (sequential) ...", B))
log_msg("   (estimate: ~%.0f minutes)", B * 0.15)

set.seed(SEED)

results_list <- vector("list", B)
t0 <- Sys.time()

for (b in seq_len(B)) {
  boot_idx <- sample(seq_len(n), size = n, replace = TRUE)
  X_boot   <- X_full[boot_idx, , drop = FALSE]
  clin_boot <- clin_cov_matrix[boot_idx, , drop = FALSE]
  y_boot   <- y_full[boot_idx]

  n_events_boot <- sum(y_boot[, "status"] == 1)

  if (n_events_boot < MIN_EVENTS) {
    results_list[[b]] <- data.frame(
      rep = b, success = FALSE,
      reason = "Too few bootstrap events",
      n_events_boot = n_events_boot,
      n_screened = NA, n_nonzero = NA, lambda = NA,
      c_boot = NA, c_original = NA, optimism = NA,
      stringsAsFactors = FALSE)
    next
  }

  fit_b <- tryCatch(
    fit_bootstrap_pipeline(X_boot, clin_boot, y_boot,
                           alpha = ALPHA, nfolds = N_CV_FOLDS),
    error = function(e) list(success = FALSE,
                             reason = paste0("Pipeline: ", e$message))
  )

  if (!isTRUE(fit_b$success)) {
    results_list[[b]] <- data.frame(
      rep = b, success = FALSE, reason = fit_b$reason,
      n_events_boot = n_events_boot,
      n_screened = NA, n_nonzero = NA, lambda = NA,
      c_boot = NA, c_original = NA, optimism = NA,
      stringsAsFactors = FALSE)
    next
  }

  # C-index on the bootstrap
  lp_boot <- predict(fit_b$model, newx = fit_b$X_boot, type = "link")
  c_boot <- tryCatch(compute_cindex(y_boot, lp_boot, use_method),
                     error = function(e) NA_real_)

  # rebuild the original train with the same bootstrap features and scaling
  selected_names <- intersect(fit_b$selected_screening_names,
                              colnames(X_full))
  if (length(selected_names) < 1) {
    results_list[[b]] <- data.frame(
      rep = b, success = FALSE,
      reason = "Selected vars absent from original",
      n_events_boot = n_events_boot,
      n_screened = fit_b$n_screened, n_nonzero = fit_b$n_nonzero,
      lambda = fit_b$lambda, c_boot = c_boot, c_original = NA,
      optimism = NA, stringsAsFactors = FALSE)
    next
  }

  X_orig_mir <- X_full[, selected_names, drop = FALSE]
  orig_mir_scaled <- apply_scaler(X_orig_mir, fit_b$mir_scaler)

  model_mir_names <- intersect(fit_b$mirna_names, colnames(orig_mir_scaled))
  if (length(model_mir_names) < 1) {
    results_list[[b]] <- data.frame(
      rep = b, success = FALSE,
      reason = "No model miRNA in original",
      n_events_boot = n_events_boot,
      n_screened = fit_b$n_screened, n_nonzero = fit_b$n_nonzero,
      lambda = fit_b$lambda, c_boot = c_boot, c_original = NA,
      optimism = NA, stringsAsFactors = FALSE)
    next
  }

  orig_mir_scaled <- orig_mir_scaled[, model_mir_names, drop = FALSE]
  orig_clin_scaled <- apply_scaler(clin_cov_matrix, fit_b$clin_scaler)

  X_original <- cbind(orig_mir_scaled, orig_clin_scaled)
  model_columns <- colnames(fit_b$X_boot)

  if (!all(model_columns %in% colnames(X_original))) {
    results_list[[b]] <- data.frame(
      rep = b, success = FALSE,
      reason = "Column mismatch",
      n_events_boot = n_events_boot,
      n_screened = fit_b$n_screened, n_nonzero = fit_b$n_nonzero,
      lambda = fit_b$lambda, c_boot = c_boot, c_original = NA,
      optimism = NA, stringsAsFactors = FALSE)
    next
  }

  X_original <- X_original[, model_columns, drop = FALSE]

  lp_original <- tryCatch(
    predict(fit_b$model, newx = X_original, type = "link"),
    error = function(e) NULL
  )
  if (is.null(lp_original)) {
    results_list[[b]] <- data.frame(
      rep = b, success = FALSE,
      reason = "Prediction on original failed",
      n_events_boot = n_events_boot,
      n_screened = fit_b$n_screened, n_nonzero = fit_b$n_nonzero,
      lambda = fit_b$lambda, c_boot = c_boot, c_original = NA,
      optimism = NA, stringsAsFactors = FALSE)
    next
  }

  c_original <- tryCatch(compute_cindex(y_full, lp_original, use_method),
                         error = function(e) NA_real_)

  results_list[[b]] <- data.frame(
    rep = b, success = TRUE, reason = NA_character_,
    n_events_boot = n_events_boot,
    n_screened = fit_b$n_screened,
    n_nonzero = fit_b$n_nonzero,
    lambda = fit_b$lambda,
    c_boot = c_boot, c_original = c_original,
    optimism = c_boot - c_original,
    stringsAsFactors = FALSE)

  if (b %% 50 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    rate <- b / max(elapsed, 0.01)
    eta  <- (B - b) / max(rate, 0.01)
    log_msg(sprintf("   %d/%d | %.1f/min | ETA %.0f min", b, B, rate, eta))
  }
}

results <- do.call(rbind, results_list)
log_msg("")
log_msg("Bootstrap finished.")

# ======================================================================
# 6) Success/failure summary
# ======================================================================
results$success <- as.logical(results$success)
n_success <- sum(results$success & is.finite(results$optimism))
n_failed  <- B - n_success

log_msg("============================================================")
log_msg("BOOTSTRAP QUALITY CONTROL")
log_msg("============================================================")
log_msg(sprintf("Requested  : %d", B))
log_msg(sprintf("Successful : %d", n_success))
log_msg(sprintf("Failed     : %d", n_failed))
log_msg(sprintf("Success %%  : %.1f%%", 100 * n_success / B))

if (n_failed > 0) {
  failure_table <- sort(table(results$reason[!results$success]),
                        decreasing = TRUE)
  log_msg("")
  log_msg("Failure reasons:")
  print(failure_table)
}

if (n_success < 500)
  warning(sprintf("Only %d successful replicates — optimism estimate may be unstable.",
                  n_success))
if (n_success == 0) stop("No successful bootstrap replicates.")

# ======================================================================
# 7) Compute the optimism
# ======================================================================
valid <- results[results$success & is.finite(results$optimism), , drop = FALSE]
optimism_vec   <- valid$optimism
mean_optimism  <- mean(optimism_vec, na.rm = TRUE)
sd_optimism    <- sd(optimism_vec, na.rm = TRUE)
se_optimism    <- sd_optimism / sqrt(n_success)
optimism_q     <- quantile(optimism_vec, probs = c(0.025, 0.50, 0.975),
                           na.rm = TRUE)

C_corrected <- max(min(C_apparent - mean_optimism, 1), 0)

mean_c_boot     <- mean(valid$c_boot, na.rm = TRUE)
mean_c_original <- mean(valid$c_original, na.rm = TRUE)
mean_screened   <- mean(valid$n_screened, na.rm = TRUE)
median_screened <- median(valid$n_screened, na.rm = TRUE)
mean_nonzero    <- mean(valid$n_nonzero, na.rm = TRUE)
median_nonzero  <- median(valid$n_nonzero, na.rm = TRUE)

# ======================================================================
# 8) Print the results
# ======================================================================
log_msg("")
log_msg("============================================================")
log_msg("FINAL BOOTSTRAP RESULTS (Signature v2)")
log_msg("============================================================")
log_msg(sprintf("C-index apparent             : %.4f", C_apparent))
log_msg(sprintf("Mean bootstrap C-index       : %.4f", mean_c_boot))
log_msg(sprintf("Mean original-train C-index  : %.4f", mean_c_original))
log_msg(sprintf("Mean optimism                : %.4f", mean_optimism))
log_msg(sprintf("SE optimism                  : %.4f", se_optimism))
log_msg(sprintf("Optimism median              : %.4f", optimism_q[2]))
log_msg(sprintf("Optimism 2.5%%–97.5%%          : %.4f – %.4f",
                optimism_q[1], optimism_q[3]))
log_msg(sprintf("C-index optimism-corrected   : %.4f", C_corrected))
log_msg("")
log_msg(sprintf("Mean screened miRNAs         : %.2f", mean_screened))
log_msg(sprintf("Median screened miRNAs       : %.0f", median_screened))
log_msg(sprintf("Mean non-zero coefficients   : %.2f", mean_nonzero))
log_msg(sprintf("Median non-zero coefficients : %.0f", median_nonzero))

# ======================================================================
# 9) Save
# ======================================================================
write.csv(results,
          file.path(sig_dir, "signature_v2_bootstrap_results.csv"),
          row.names = FALSE)

if (n_failed > 0) {
  failure_df <- data.frame(
    reason = names(failure_table),
    n = as.integer(failure_table),
    stringsAsFactors = FALSE)
  write.csv(failure_df,
            file.path(sig_dir, "signature_v2_bootstrap_failure_reasons.csv"),
            row.names = FALSE)
}

summary_lines <- c(
  "============================================================",
  "STEP 12 (Signature v2) — BOOTSTRAP OPTIMISM CORRECTION",
  "============================================================",
  paste("Execution date:", as.character(Sys.time())),
  "",
  paste("Training patients   :", n),
  paste("Training events     :", sum(y_full[, "status"])),
  paste("Candidate miRNAs    :", ncol(X_full)),
  paste("Clinical variables  :", ncol(clin_cov_matrix)),
  "",
  paste("Requested replicates:", B),
  paste("Successful          :", n_success),
  paste("Failed              :", n_failed),
  paste("Success rate        :", sprintf("%.2f%%", 100 * n_success / B)),
  "",
  paste("Univariate threshold strict:", UNIV_P_STRICT),
  paste("Univariate threshold relax :", UNIV_P_RELAX),
  paste("Elastic Net alpha          :", ALPHA),
  paste("CV folds                   :", N_CV_FOLDS),
  "",
  paste("Apparent C-index           :", round(C_apparent, 4)),
  paste("Mean bootstrap C-index     :", round(mean_c_boot, 4)),
  paste("Mean original-train C-index:", round(mean_c_original, 4)),
  paste("Mean optimism              :", round(mean_optimism, 4)),
  paste("SE optimism                :", round(se_optimism, 4)),
  paste("Optimism 2.5%-97.5%        :",
        round(optimism_q[1], 4), "to", round(optimism_q[3], 4)),
  "",
  paste("Optimism-corrected C-index :", round(C_corrected, 4)),
  "",
  paste("Mean screened miRNAs       :", round(mean_screened, 2)),
  paste("Median screened miRNAs     :", round(median_screened, 0)),
  paste("Mean non-zero coefficients :", round(mean_nonzero, 2)),
  paste("Median non-zero coefficients:", round(median_nonzero, 0)),
  "",
  "IMPORTANT:",
  "  - The 2.5%-97.5% range is the empirical distribution of bootstrap",
  "    optimism, NOT a formal 95% CI for the corrected C-index.",
  "  - Test set was NOT used at any stage of this procedure."
)
writeLines(summary_lines,
           file.path(sig_dir, "signature_v2_step12_summary.txt"))

bootstrap_summary <- list(
  method         = "Bootstrap optimism correction (Signature v2)",
  B_requested    = B,
  B_success      = n_success,
  B_failed       = n_failed,
  C_apparent     = C_apparent,
  mean_optimism  = mean_optimism,
  sd_optimism    = sd_optimism,
  se_optimism    = se_optimism,
  optimism_quantiles = optimism_q,
  C_corrected    = C_corrected,
  mean_c_boot    = mean_c_boot,
  mean_c_original = mean_c_original,
  mean_screened  = mean_screened,
  median_screened = median_screened,
  mean_nonzero   = mean_nonzero,
  median_nonzero = median_nonzero,
  univ_p_strict  = UNIV_P_STRICT,
  univ_p_relax   = UNIV_P_RELAX,
  alpha          = ALPHA,
  n_cv_folds     = N_CV_FOLDS,
  use_method     = use_method,
  seed           = SEED,
  results        = results
)
saveRDS(bootstrap_summary,
        file.path(sig_dir, "signature_v2_bootstrap_optimism_corrected_cindex.rds"))

log_msg("")
log_msg("outputs in: ", sig_dir)
log_msg("STEP 12 (Signature v2) COMPLETED")
