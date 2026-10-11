#!/usr/bin/env Rscript
# =============================================================================
# STEP 08 (Signature v2) — Signature Lock
# -----------------------------------------------------------------------------
# Freeze the v2 signature (6 miRNA + Age) with a SHA-256 hash
# Outputs are saved in signature_v2/ - no interference with the original version
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(digest)
  library(jsonlite)
})

output_dir <- paste0(WORK_DIR, "/")
sig_dir    <- file.path(output_dir, "signature_v2")

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

# ======================================================================
# 1) Load the Step 07 (v2) outputs
# ======================================================================
log_msg("Loading Step 07 (v2) outputs ...")

model_path   <- file.path(sig_dir, "signature_v2_model.rds")
feat_path    <- file.path(sig_dir, "signature_features.csv")
perf_path    <- file.path(sig_dir, "signature_v2_cv_performance.csv")
scaling_path <- file.path(sig_dir, "signature_v2_scaling.rds")
manifest7_path <- file.path(sig_dir, "step07_signature_v2_manifest.json")

for (p in c(model_path, feat_path, perf_path, scaling_path)) {
  if (!file.exists(p)) stop("not found: ", p)
}

final_model       <- readRDS(model_path)
selected_features <- read.csv(feat_path, stringsAsFactors = FALSE)
cv_perf           <- read.csv(perf_path)
scaling_params    <- readRDS(scaling_path)

manifest7 <- if (file.exists(manifest7_path)) {
  fromJSON(manifest7_path, simplifyVector = FALSE)
} else NULL

log_msg(sprintf("model, %d features, %d CV folds loaded.",
                nrow(selected_features), nrow(cv_perf)))

# ======================================================================
# 2) Separate miRNA from clinical features — using scaling_params
# ======================================================================
clinical_vars  <- scaling_params$clinical$features
mirna_in_sig   <- intersect(selected_features$Feature,
                            scaling_params$mirna$features)
clinical_in_sig <- intersect(selected_features$Feature, clinical_vars)

all_scaling_features <- c(scaling_params$mirna$features,
                          scaling_params$clinical$features)
missing_in_scaling <- setdiff(selected_features$Feature, all_scaling_features)
if (length(missing_in_scaling) > 0) {
  stop("features in the model that are missing from scaling_params: ",
       paste(missing_in_scaling, collapse = ", "))
}

log_msg(sprintf("   - miRNA: %d | clinical: %d (%s)",
                length(mirna_in_sig),
                length(clinical_in_sig),
                ifelse(length(clinical_in_sig) == 0, "none",
                       paste(clinical_in_sig, collapse = ", "))))

if (length(mirna_in_sig) == 0) {
  warning("no miRNA in the signature.")
}

# ======================================================================
# 3) Structural checks of scaling_params
# ======================================================================
if (length(scaling_params$mirna$center) != length(scaling_params$mirna$features))
  stop("length of mirna center does not match features.")
if (length(scaling_params$mirna$scale) != length(scaling_params$mirna$features))
  stop("length of mirna scale does not match features.")
if (length(scaling_params$clinical$center) != length(scaling_params$clinical$features))
  stop("length of clinical center does not match features.")
if (length(scaling_params$clinical$scale) != length(scaling_params$clinical$features))
  stop("length of clinical scale does not match features.")

if (any(scaling_params$mirna$scale == 0, na.rm = TRUE))
  stop("zero scale in miRNA.")
if (any(scaling_params$clinical$scale == 0, na.rm = TRUE))
  stop("zero scale in clinical.")

# check that the model matches the scaling
coef_mat <- as.matrix(coef(final_model))
model_features_nz <- rownames(coef_mat)[coef_mat[, 1] != 0]

missing_from_model <- setdiff(model_features_nz, all_scaling_features)
if (length(missing_from_model) > 0) {
  stop("the model has features that are missing from scaling_params: ",
       paste(missing_from_model, collapse = ", "))
}

if (length(model_features_nz) != nrow(selected_features)) {
  warning(sprintf("number of non-zero model features (%d) differs from selected_features (%d).",
                  length(model_features_nz), nrow(selected_features)))
}

log_msg("model and scaling params are consistent.")

# ======================================================================
# 4) Build the signature object
# ======================================================================
coef_final <- as.matrix(coef(final_model))
coef_df <- data.frame(
  Feature     = rownames(coef_final),
  Coefficient = coef_final[, 1],
  Type        = ifelse(rownames(coef_final) %in% clinical_vars,
                       "clinical", "miRNA_isoform"),
  stringsAsFactors = FALSE
)
coef_df <- coef_df[coef_df$Coefficient != 0, ]
rownames(coef_df) <- NULL

cv_summary <- list(
  mean_cindex     = mean(cv_perf$c_index),
  sd_cindex       = sd(cv_perf$c_index),
  median_cindex   = median(cv_perf$c_index),
  min_cindex      = min(cv_perf$c_index),
  max_cindex      = max(cv_perf$c_index),
  n_valid_folds   = nrow(cv_perf)
)

cohort_summary <- if (!is.null(manifest7)) {
  list(
    n_train             = manifest7$n_train,
    n_events            = manifest7$n_events,
    n_mirna_input       = manifest7$n_mirna_input,
    n_mirna_prefiltered = manifest7$n_mirna_input,
    clinical_vars       = manifest7$clinical_features,
    epv_min             = manifest7$epv_min,
    freq_threshold      = manifest7$freq_threshold
  )
} else list()

signature_object <- list(
  meta = list(
    pipeline         = "miR-Mito Score - TCGA-GBM (Signature v2)",
    step             = "08_signature_lock_v2",
    locked_on        = as.character(Sys.time()),
    r_version        = R.version.string,
    package_versions = list(
      glmnet   = as.character(packageVersion("glmnet")),
      survival = as.character(packageVersion("survival")),
      digest   = as.character(packageVersion("digest"))
    )
  ),
  model                  = final_model,
  coefficients           = coef_df,
  mirna_features         = mirna_in_sig,
  clinical_features      = clinical_in_sig,
  scaling                = scaling_params,
  cv_performance_summary = cv_summary,
  cohort_summary         = cohort_summary,
  alpha                  = scaling_params$alpha,
  survival_vars          = list(
    status_col = "status_final",
    time_col   = "survival_time_final"
  )
)

# ======================================================================
# 5) SHA-256 hash
# ======================================================================
hash_payload <- list(
  coefficients      = signature_object$coefficients,
  mirna_features    = signature_object$mirna_features,
  clinical_features = signature_object$clinical_features,
  scaling           = signature_object$scaling,
  alpha             = signature_object$alpha
)
sig_hash       <- digest(hash_payload, algo = "sha256")
sig_hash_short <- substr(sig_hash, 1, 12)
locked_filename <- paste0("signature_v2_locked_", sig_hash_short, ".rds")
locked_path     <- file.path(sig_dir, locked_filename)

log_msg(sprintf("SHA-256 hash: %s", sig_hash))
log_msg(sprintf("   lock file: %s", locked_filename))

# ======================================================================
# 6) Save + lock
# ======================================================================
if (file.exists(locked_path)) {
  log_msg("A file with the same hash already exists — not saving again.")
} else {
  saveRDS(signature_object, locked_path)
  tryCatch(Sys.chmod(locked_path, mode = "0444"),
           error = function(e) NULL)
  log_msg(sprintf("locked: %s (read-only)", locked_filename))
}

# check for a previous v2 version
latest_path <- file.path(sig_dir, "signature_v2_locked_latest.rds")
if (file.exists(latest_path)) {
  old_sig <- tryCatch(readRDS(latest_path), error = function(e) NULL)
  if (!is.null(old_sig)) {
    old_hash <- tryCatch({
      digest(list(
        coefficients      = old_sig$coefficients,
        mirna_features    = old_sig$mirna_features,
        clinical_features = old_sig$clinical_features,
        scaling           = old_sig$scaling,
        alpha             = old_sig$alpha
      ), algo = "sha256")
    }, error = function(e) NA_character_)

    if (!is.na(old_hash) && old_hash != sig_hash) {
      log_msg("signature_v2_locked_latest.rds differs from this signature!")
      log_msg(sprintf("   old hash: %s", substr(old_hash, 1, 12)))
      log_msg(sprintf("   new hash: %s", sig_hash_short))
      backup_path <- file.path(sig_dir,
        paste0("signature_v2_locked_latest.BACKUP_",
               substr(old_hash, 1, 8), ".rds"))
      file.copy(latest_path, backup_path, overwrite = TRUE)
      log_msg(sprintf("   old version backed up: %s", basename(backup_path)))
    }
  }
}

saveRDS(signature_object, latest_path)
log_msg("signature_v2_locked_latest.rds saved.")

# ======================================================================
# 7) signature_v2_formula.txt
# ======================================================================
coef_lines <- sprintf("  %-28s  beta = %+9.5f   (%s)",
                      coef_df$Feature, coef_df$Coefficient, coef_df$Type)

formula_lines <- c(
  "==================== Signature v2 formula (miR-Mito Score, TCGA-GBM) ====================",
  paste("Lock date:    ", signature_object$meta$locked_on),
  paste("SHA-256 hash: ", sig_hash),
  paste("Lock file:    ", locked_filename),
  "",
  sprintf("Number of features: %d  (miRNA: %d | clinical: %d)",
          nrow(coef_df), length(mirna_in_sig), length(clinical_in_sig)),
  "",
  sprintf("Design parameters: freq_threshold=%.2f | EPV_MIN=%d",
          manifest7$freq_threshold %||% NA,
          manifest7$epv_min %||% NA),
  "",
  "Cox model coefficients (on the z-score scale):",
  "",
  coef_lines,
  "",
  "--- Performance summary (nested CV on train) ---",
  sprintf("  Mean C-index: %.3f ± %.3f",
          cv_summary$mean_cindex, cv_summary$sd_cindex),
  sprintf("  Range: [%.3f, %.3f]   |   valid folds: %d",
          cv_summary$min_cindex, cv_summary$max_cindex, cv_summary$n_valid_folds),
  "",
  "--- Cohort summary ---",
  if (!is.null(manifest7))
    sprintf("  n_train = %d | n_events = %d",
            manifest7$n_train, manifest7$n_events)
  else "  (Step 07 manifest not available)",
  "",
  "This is a nested-CV number, not performance on an independent test set.",
  "   Official test performance is reported only in Step 09.",
  "",
  "Reproducibility note:",
  "  - miRNAs and Age must be normalised with the center/scale stored in 'scaling'.",
  "  - This locked file must not be read back before Step 09.",
  "  - Any change to the pipeline requires re-running Steps 07 and 08."
)

formula_path <- file.path(sig_dir, "signature_v2_formula.txt")
writeLines(formula_lines, formula_path)
log_msg("signature_v2_formula.txt saved.")

# ======================================================================
# 8) Step 08 manifest (v2)
# ======================================================================
manifest8 <- list(
  step                 = "08_signature_lock_v2",
  locked_on            = as.character(Sys.time()),
  signature_hash       = sig_hash,
  signature_hash_short = sig_hash_short,
  locked_file          = locked_filename,
  latest_file          = "signature_v2_locked_latest.rds",
  formula_file         = "signature_v2_formula.txt",
  model_file           = "signature_v2_model.rds",
  scaling_file         = "signature_v2_scaling.rds",

  n_features_total     = nrow(coef_df),
  n_mirna_features     = length(mirna_in_sig),
  n_clinical_features  = length(clinical_in_sig),
  mirna_features       = mirna_in_sig,
  clinical_features    = clinical_in_sig,
  coefficients         = setNames(coef_df$Coefficient, coef_df$Feature),

  cv_cindex_mean       = cv_summary$mean_cindex,
  cv_cindex_sd         = cv_summary$sd_cindex,
  cv_n_folds           = cv_summary$n_valid_folds,

  cohort_n_train       = if (!is.null(manifest7)) manifest7$n_train else NA,
  cohort_n_events      = if (!is.null(manifest7)) manifest7$n_events else NA,
  freq_threshold       = if (!is.null(manifest7)) manifest7$freq_threshold else NA,
  epv_min              = if (!is.null(manifest7)) manifest7$epv_min else NA,

  r_version            = R.version.string
)

manifest8_path <- file.path(sig_dir, "step08_signature_v2_manifest.json")
write_json(manifest8, manifest8_path, pretty = TRUE, auto_unbox = TRUE)
log_msg("step08_signature_v2_manifest.json saved.")

# ======================================================================
# 9) Final summary
# ======================================================================
log_msg("")
log_msg("==================== Summary of the locked v2 signature ====================")
log_msg(sprintf("hash: %s", sig_hash_short))
log_msg(sprintf("lock file: %s", locked_filename))
log_msg(sprintf("features (%d total):", nrow(coef_df)))
for (i in seq_len(nrow(coef_df))) {
  log_msg(sprintf("   %-28s %+9.5f  (%s)",
                  coef_df$Feature[i], coef_df$Coefficient[i], coef_df$Type[i]))
}
log_msg(sprintf("Nested CV C-index: %.3f ± %.3f",
                cv_summary$mean_cindex, cv_summary$sd_cindex))
log_msg("")
log_msg("Step 08 (v2) complete. Signature v2 locked.")
