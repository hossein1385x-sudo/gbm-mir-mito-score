#!/usr/bin/env Rscript
# =============================================================================
# STEP 07 (Signature v2) — LASSO-Cox + Nested CV + Final Model
# -----------------------------------------------------------------------------
# This script is a self-contained pipeline that builds the 6-miRNA + Age signature.
# Parameters:
#   consensus threshold : 0.60
#   EPV min             : 10
#
# Outputs (in the signature_v2/ folder):
#   signature_features.csv
#   signature_v2_model.rds
#   signature_v2_scaling.rds
#   signature_v2_cv_performance.csv
#   signature_v2_selection_frequency.csv
#   step07_signature_v2_manifest.json
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(jsonlite)
  library(data.table)
})

output_dir <- paste0(WORK_DIR, "/")
sig_dir    <- file.path(output_dir, "signature_v2")
dir.create(sig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Parameters ----
FREQ_THRESHOLD <- 0.60
EPV_MIN        <- 10
ALPHA          <- 0.5
N_REPEATS      <- 10
N_OUTER        <- 5
N_INNER        <- 5

log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

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
# 1) Load the data
# ======================================================================
log_msg("============================================================")
log_msg("STEP 07 (Signature v2) — Training")
log_msg("============================================================")
log_msg("loading ...")

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
mirna_meta[, sample_clean := normalize_tcga_id(sample_barcode)]
primary_barcodes <- unique(mirna_meta[is_primary_tumor == TRUE, sample_clean])
primary_barcodes <- primary_barcodes[!is.na(primary_barcodes)]

# ======================================================================
# 2) Synchronisation
# ======================================================================
iso_ids  <- normalize_tcga_id(rownames(isoform_data))
clin_ids <- normalize_tcga_id(as.character(clinical_data$patient_barcode))

common <- Reduce(intersect, list(iso_ids, clin_ids, train_ids, primary_barcodes))
log_msg(sprintf("   shared train samples: %d", length(common)))

idx_iso  <- match(common, iso_ids)
isoform_data <- isoform_data[idx_iso, , drop = FALSE]
rownames(isoform_data) <- common

idx_clin <- match(common, clin_ids)
clin <- as.data.frame(clinical_data[idx_clin, , drop = FALSE])
rownames(clin) <- common

# ======================================================================
# 3) log2 + survival
# ======================================================================
X_full <- log2(as.matrix(isoform_data) + 1)

clin$status        <- clin$status_final
clin$survival_time <- clin$survival_time_final

keep <- !is.na(clin$status) & !is.na(clin$survival_time) & clin$survival_time >= 0
clin   <- clin[keep, , drop = FALSE]
X_full <- X_full[rownames(clin), , drop = FALSE]

keep_rows <- complete.cases(X_full)
X_full <- X_full[keep_rows, , drop = FALSE]
clin   <- clin[keep_rows, , drop = FALSE]

y_full <- Surv(clin$survival_time, clin$status)
log_msg(sprintf("train: %d patients | %d events",
                nrow(clin), sum(clin$status == 1)))

# ======================================================================
# 4) Expression pre-filter (train only)
# ======================================================================
log_msg("")
log_msg("expression pre-filter ...")

n_samples <- nrow(X_full)
X_rpm_raw <- as.matrix(isoform_data[rownames(clin), , drop = FALSE])

# Filter 1: RPM > 1 in ≥50% of samples
keep_expr <- colSums(X_rpm_raw > 1, na.rm = TRUE) >= 0.5 * n_samples
log_msg(sprintf("  expression filter: %d → %d", ncol(X_rpm_raw), sum(keep_expr)))

# Filter 2: log2 variance above the 20th percentile
X_log <- X_full[, keep_expr, drop = FALSE]
vars <- apply(X_log, 2, var, na.rm = TRUE)
keep_var <- vars > quantile(vars, 0.2, na.rm = TRUE)
log_msg(sprintf("  variance filter: %d → %d", sum(keep_expr), sum(keep_var)))

keep_final <- colnames(X_log)[keep_var]
X_full <- X_full[, keep_final, drop = FALSE]
log_msg(sprintf("  final isoforms: %d", ncol(X_full)))

# ======================================================================
# 5) Age
# ======================================================================
age_raw <- suppressWarnings(as.numeric(clin$age_at_diagnosis))
if (all(is.na(age_raw))) age_raw <- suppressWarnings(as.numeric(clin$age))
if (median(age_raw, na.rm = TRUE) > 200) age_years <- age_raw / 365.25 else age_years <- age_raw
age_years[is.na(age_years)] <- median(age_years, na.rm = TRUE)

clin_cov <- data.frame(Age = age_years)
rownames(clin_cov) <- rownames(clin)

# ======================================================================
# 6) Nested CV (repeated)
# ======================================================================
X_matrix <- as.matrix(X_full)
n <- nrow(X_matrix)

log_msg("")
log_msg(sprintf("%d repeats × %d outer folds ...", N_REPEATS, N_OUTER))

run_one_repeat <- function(seed) {
  set.seed(seed)

  foldid <- integer(n)
  for (ev in c(0, 1)) {
    idx <- which(y_full[, "status"] == ev)
    foldid[idx] <- sample(rep(1:N_OUTER, length.out = length(idx)))
  }

  c_index_outer <- c()
  selected_mirna_names <- list()

  for (fold in 1:N_OUTER) {
    test_idx  <- which(foldid == fold)
    train_idx <- setdiff(1:n, test_idx)
    y_train   <- y_full[train_idx]
    y_test    <- y_full[test_idx]

    univ_p <- apply(X_matrix[train_idx, , drop = FALSE], 2, function(col) {
      if (sd(col, na.rm = TRUE) == 0) return(NA)
      fit <- tryCatch(coxph(y_train ~ col), error = function(e) NULL)
      if (!is.null(fit)) summary(fit)$coefficients[5] else NA
    })

    keep_uni <- which(!is.na(univ_p) & univ_p < 0.01)
    if (length(keep_uni) < 20) {
      keep_uni <- which(!is.na(univ_p) & univ_p < 0.05)
    }
    if (length(keep_uni) < 2) next

    sub_tr <- X_matrix[train_idx, keep_uni, drop = FALSE]
    tr_mean <- colMeans(sub_tr)
    tr_sd   <- apply(sub_tr, 2, sd, na.rm = TRUE)
    valid   <- !is.na(tr_sd) & tr_sd > 0
    if (sum(valid) < 2) next
    keep_final_fold <- keep_uni[valid]

    X_tr_mir <- scale(X_matrix[train_idx, keep_final_fold, drop = FALSE],
                      center = tr_mean[valid], scale = tr_sd[valid])
    X_te_mir <- scale(X_matrix[test_idx, keep_final_fold, drop = FALSE],
                      center = tr_mean[valid], scale = tr_sd[valid])

    X_tr_clin <- scale(clin_cov[train_idx, , drop = FALSE],
                       center = TRUE, scale = TRUE)
    X_te_clin <- scale(clin_cov[test_idx, , drop = FALSE],
                       center = attr(X_tr_clin, "scaled:center"),
                       scale  = attr(X_tr_clin, "scaled:scale"))

    X_tr <- cbind(X_tr_mir, X_tr_clin)
    X_te <- cbind(X_te_mir, X_te_clin)

    n_ev_train <- sum(y_train[, "status"] == 1)
    if (n_ev_train < N_INNER) next
    if (sum(y_test[, "status"] == 1) == 0) next

    cv_fit <- tryCatch(
      cv.glmnet(X_tr, y_train, family = "cox", alpha = ALPHA,
                nfolds = min(N_INNER, n_ev_train)),
      error = function(e) NULL)
    if (is.null(cv_fit)) next

    fold_model <- glmnet(X_tr, y_train, family = "cox", alpha = ALPHA,
                         lambda = cv_fit$lambda.min)
    lp_test <- as.numeric(predict(fold_model, newx = X_te, type = "link"))

    c_index <- if (length(unique(lp_test)) > 1) {
      survConcordance(y_test ~ lp_test)$concordance
    } else 0.5

    c_index_outer <- c(c_index_outer, c_index)

    coefs <- as.matrix(coef(fold_model))
    nz    <- rownames(coefs)[coefs[, 1] != 0]
    mirna_nz <- setdiff(nz, colnames(clin_cov))
    selected_mirna_names[[length(selected_mirna_names) + 1]] <- mirna_nz
  }

  list(c_index = c_index_outer, selected = selected_mirna_names)
}

all_results <- lapply(1:N_REPEATS, function(r) run_one_repeat(seed = 100 + r))
all_cindex  <- unlist(lapply(all_results, `[[`, "c_index"))
log_msg(sprintf("Nested CV C-index: %.3f ± %.3f (n=%d)",
                mean(all_cindex), sd(all_cindex), length(all_cindex)))

# ======================================================================
# 7) miRNA selection frequency
# ======================================================================
all_selected_lists <- unlist(lapply(all_results, `[[`, "selected"), recursive = FALSE)
n_valid <- length(all_selected_lists)

freq_tbl <- table(unlist(all_selected_lists))
sel_freq <- data.frame(
  miRNA      = names(freq_tbl),
  n_selected = as.integer(freq_tbl),
  frequency  = as.numeric(freq_tbl) / n_valid,
  stringsAsFactors = FALSE
)
sel_freq <- sel_freq[order(-sel_freq$frequency), ]

write.csv(sel_freq,
          file.path(sig_dir, "signature_v2_selection_frequency.csv"),
          row.names = FALSE)

log_msg("")
log_msg("most frequently selected miRNAs:")
print(head(sel_freq, 15))

# ======================================================================
# 8) Final model on the full train set
# ======================================================================
n_events    <- sum(clin$status == 1)
n_clin_vars <- ncol(clin_cov)
max_mirna_by_epv <- max(1, floor(n_events / EPV_MIN) - n_clin_vars)

log_msg("")
log_msg(sprintf("EPV_MIN=%d | events=%d | max_mirna=%d",
                EPV_MIN, n_events, max_mirna_by_epv))

consensus <- sel_freq$miRNA[sel_freq$frequency >= FREQ_THRESHOLD]

if (length(consensus) > max_mirna_by_epv) {
  consensus <- head(consensus, max_mirna_by_epv)
  log_msg(sprintf("   EPV cap applied → %d miRNAs", length(consensus)))
} else if (length(consensus) < 2) {
  n_take <- min(max_mirna_by_epv, max(2, nrow(sel_freq)))
  consensus <- head(sel_freq$miRNA, n_take)
}

log_msg(sprintf("signature: %d miRNAs", length(consensus)))
print(consensus)

# build the final model
X_mir_scaled <- scale(X_matrix[, consensus, drop = FALSE],
                      center = TRUE, scale = TRUE)
mir_center <- attr(X_mir_scaled, "scaled:center")
mir_scale  <- attr(X_mir_scaled, "scaled:scale")

na_cols <- apply(X_mir_scaled, 2, function(x) any(is.na(x)))
if (any(na_cols)) {
  X_mir_scaled <- X_mir_scaled[, !na_cols, drop = FALSE]
  mir_center   <- mir_center[!na_cols]
  mir_scale    <- mir_scale[!na_cols]
}

X_clin_scaled <- scale(clin_cov, center = TRUE, scale = TRUE)
clin_center <- attr(X_clin_scaled, "scaled:center")
clin_scale  <- attr(X_clin_scaled, "scaled:scale")

X_full_scaled <- cbind(X_mir_scaled, X_clin_scaled)

set.seed(42)
final_cv <- cv.glmnet(X_full_scaled, y_full, family = "cox",
                      alpha = ALPHA, nfolds = min(10, n_events))
final_model <- glmnet(X_full_scaled, y_full, family = "cox",
                      alpha = ALPHA, lambda = final_cv$lambda.min)

coef_final <- as.matrix(coef(final_model))
selected <- coef_final[coef_final[, 1] != 0, , drop = FALSE]
selected_names <- rownames(selected)

lp_full <- as.numeric(predict(final_model, newx = X_full_scaled, type = "link"))
c_train <- survConcordance(y_full ~ lp_full)$concordance

log_msg("")
log_msg(sprintf("final features: %d | C-index (train): %.3f",
                nrow(selected), c_train))
print(selected)

# ======================================================================
# 9) miRNA-only analysis (on train)
# ======================================================================
mirna_in_sig    <- intersect(selected_names, colnames(X_mir_scaled))
clinical_in_sig <- intersect(selected_names, colnames(clin_cov))

log_msg("")
log_msg(sprintf("miRNA-only score (%d miRNA) ...", length(mirna_in_sig)))

mir_coef <- selected[mirna_in_sig, 1]
score_mir_only <- as.numeric(X_mir_scaled[, mirna_in_sig, drop = FALSE] %*% mir_coef)

conc_mir  <- concordance(y_full ~ score_mir_only, reverse = TRUE)
c_mir_only <- as.numeric(conc_mir$concordance)

cox_mir <- coxph(y_full ~ scale(score_mir_only))
s_mir <- summary(cox_mir)
beta_m <- s_mir$coefficients[1, "coef"]
se_m   <- s_mir$coefficients[1, "se(coef)"]
p_m    <- s_mir$coefficients[1, "Pr(>|z|)"]

log_msg(sprintf("   C-index (miRNA-only): %.4f", c_mir_only))
log_msg(sprintf("   HR per SD: %.3f (95%% CI: %.3f–%.3f, p=%.4g)",
                exp(beta_m),
                exp(beta_m - 1.96 * se_m),
                exp(beta_m + 1.96 * se_m),
                p_m))

# ======================================================================
# 10) Save the model + scaling params
# ======================================================================
write.csv(data.frame(Feature = selected_names, Coefficient = selected[, 1]),
          file.path(sig_dir, "signature_features.csv"),
          row.names = FALSE)

write.csv(data.frame(c_index = all_cindex),
          file.path(sig_dir, "signature_v2_cv_performance.csv"),
          row.names = FALSE)

saveRDS(final_model,
        file.path(sig_dir, "signature_v2_model.rds"))

scaling_params <- list(
  mirna = list(
    center   = mir_center,
    scale    = mir_scale,
    features = colnames(X_mir_scaled)
  ),
  clinical = list(
    center   = clin_center,
    scale    = clin_scale,
    features = colnames(clin_cov)
  ),
  alpha = ALPHA,
  mirna_features    = mirna_in_sig,
  clinical_features = clinical_in_sig,
  selected_features = selected_names
)
saveRDS(scaling_params,
        file.path(sig_dir, "signature_v2_scaling.rds"))

manifest <- list(
  step              = "step07_signature_v2",
  freq_threshold    = FREQ_THRESHOLD,
  epv_min           = EPV_MIN,
  alpha             = ALPHA,
  n_repeats         = N_REPEATS,
  n_train           = nrow(clin),
  n_events          = n_events,
  n_mirna_input     = ncol(X_matrix),
  max_mirna_by_epv  = max_mirna_by_epv,
  consensus_mirnas  = consensus,
  final_features    = selected_names,
  mirna_features    = mirna_in_sig,
  clinical_features = clinical_in_sig,
  cv_cindex_mean    = mean(all_cindex),
  cv_cindex_sd      = sd(all_cindex),
  cv_n_folds        = length(all_cindex),
  train_cindex      = c_train,
  mirna_only_cindex = c_mir_only,
  mirna_only_hr     = exp(beta_m),
  mirna_only_hr_ci  = c(exp(beta_m - 1.96*se_m), exp(beta_m + 1.96*se_m)),
  mirna_only_p      = p_m,
  timestamp         = as.character(Sys.time()),
  r_version         = R.version.string
)
write_json(manifest,
           file.path(sig_dir, "step07_signature_v2_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE)

# ======================================================================
# 11) Summary
# ======================================================================
log_msg("")
log_msg("============================================================")
log_msg("Summary of Signature v2")
log_msg("============================================================")
log_msg(sprintf("consensus miRNAs: %d", length(consensus)))
log_msg(sprintf("final features: %d", nrow(selected)))
log_msg(sprintf("C-index (train full): %.3f", c_train))
log_msg(sprintf("C-index (train miRNA-only): %.4f", c_mir_only))
log_msg(sprintf("HR miRNA-only (train): %.3f (p=%.4g)", exp(beta_m), p_m))
log_msg(sprintf("Nested CV: %.3f ± %.3f", mean(all_cindex), sd(all_cindex)))
log_msg("")
log_msg("outputs in: ", sig_dir)
log_msg("STEP 07 (Signature v2) COMPLETED")
