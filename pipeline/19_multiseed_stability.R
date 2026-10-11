#!/usr/bin/env Rscript
# =============================================================================
# STEP 18 (Signature v2) — Multi-Seed Stability Analysis
# -----------------------------------------------------------------------------
# Stability analysis of the v2 signature with 50 repeated 80/20 splits
#
# Inputs:
#   train_set.csv, mirna_isoform_rpm_matrix.csv, raw_clinical_cdr_integrated.rds
#   signature_v2/signature_v2_locked_latest.rds
#
# Outputs in signature_v2/18_multi_seed/
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(data.table)
})

set.seed(2026)

# ======================================================================
# 0) Settings
# ======================================================================
output_dir <- paste0(WORK_DIR, "/")
sig_dir    <- file.path(output_dir, "signature_v2")
out_dir    <- file.path(sig_dir, "18_multi_seed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_seeds        <- 50
seeds          <- 1001:(1000 + n_seeds)
train_fraction <- 0.80
alpha          <- 0.5
screen_p       <- 0.10
max_miRNAs     <- 10

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

norm_id <- function(x) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  substr(x, 1, 12)
}

# ======================================================================
# 1) Main paper signature (10-miRNA)
# ======================================================================
original_signature <- c(
  "hsa-miR-181a-2-3p",
  "hsa-miR-196a-5p",
  "hsa-miR-4497-5p",
  "hsa-miR-4766-5p",
  "hsa-miR-4420-3p",
  "hsa-miR-4773-3p",
  "hsa-miR-4668-5p",
  "hsa-miR-219a-1-3p",
  "hsa-miR-450a-2-3p",
  "hsa-miR-4786-5p"
)

log_msg("============================================================")
log_msg("STEP 18 (Signature v2) — Multi-Seed Stability Analysis")
log_msg("============================================================")
log_msg("")
log_msg("Original paper signature (10 miRNAs):")
print(original_signature)

# ======================================================================
# 2) v2 signature (6-miRNA + Age) as the second reference
# ======================================================================
sig_v2 <- readRDS(file.path(sig_dir, "signature_v2_locked_latest.rds"))
signature_v2_mirnas <- sig_v2$mirna_features

log_msg("")
log_msg(sprintf("Signature v2 (%d miRNAs):", length(signature_v2_mirnas)))
print(signature_v2_mirnas)

# ======================================================================
# 3) Rebuild train X + clinical + y
# ======================================================================
log_msg("")
log_msg("rebuilding the train data ...")

isoform  <- read.csv(file.path(output_dir, "mirna_isoform_rpm_matrix.csv"),
                     row.names = 1, check.names = FALSE)
clin     <- readRDS(file.path(output_dir, "raw_clinical_cdr_integrated.rds"))
train_df <- fread(file.path(output_dir, "train_set.csv"))

# primary tumor filter
mirna_raw <- readRDS(file.path(output_dir, "raw_mirna.rds"))
if (!is.data.table(mirna_raw)) mirna_raw <- as.data.table(mirna_raw)
primary <- unique(norm_id(mirna_raw[is_primary_tumor == TRUE, sample_barcode]))
primary <- primary[!is.na(primary)]

# synchronisation
train_keys <- norm_id(train_df$patient_barcode)
iso_keys   <- norm_id(rownames(isoform))
clin_keys  <- norm_id(as.character(clin$patient_barcode))

common <- Reduce(intersect, list(train_keys, iso_keys, clin_keys, primary))
log_msg(sprintf("   shared samples: %d", length(common)))

idx_iso  <- match(common, iso_keys)
idx_clin <- match(common, clin_keys)

X_raw <- log2(as.matrix(isoform[idx_iso, , drop = FALSE]) + 1)
rownames(X_raw) <- common

clin_tr <- as.data.frame(clin[idx_clin, , drop = FALSE])
rownames(clin_tr) <- common

clin_tr$status        <- clin_tr$status_final
clin_tr$survival_time <- clin_tr$survival_time_final

# survival filter
keep <- !is.na(clin_tr$status) & !is.na(clin_tr$survival_time) &
        clin_tr$survival_time >= 0
clin_tr <- clin_tr[keep, , drop = FALSE]
X_raw   <- X_raw[rownames(clin_tr), , drop = FALSE]

# full filter
keep_rows <- complete.cases(X_raw)
X_raw   <- X_raw[keep_rows, , drop = FALSE]
clin_tr <- clin_tr[keep_rows, , drop = FALSE]

# ======================================================================
# 4) Expression pre-filter (same as Step 07 v2)
# ======================================================================
log_msg("")
log_msg("expression pre-filter ...")

n_samples <- nrow(X_raw)
X_rpm_raw <- as.matrix(isoform[rownames(clin_tr), , drop = FALSE])

# Filter 1: RPM > 1 in ≥50% of samples
keep_expr <- colSums(X_rpm_raw > 1, na.rm = TRUE) >= 0.5 * n_samples
X_log     <- X_raw[, keep_expr, drop = FALSE]

# Filter 2: log2 variance above the 20th percentile
vars      <- apply(X_log, 2, var, na.rm = TRUE)
keep_var  <- vars > quantile(vars, 0.2, na.rm = TRUE)
keep_final <- colnames(X_log)[keep_var]

X_full <- X_raw[, keep_final, drop = FALSE]
log_msg(sprintf("   final isoforms: %d", ncol(X_full)))

# ======================================================================
# 5) Age (the only clinical variable)
# ======================================================================
age_raw <- suppressWarnings(as.numeric(clin_tr$age_at_diagnosis))
if (all(is.na(age_raw))) age_raw <- suppressWarnings(as.numeric(clin_tr$age))
age_years <- if (median(age_raw, na.rm = TRUE) > 200) age_raw / 365.25 else age_raw
age_years[is.na(age_years)] <- median(age_years, na.rm = TRUE)

clin_cov_matrix <- matrix(age_years, ncol = 1,
                          dimnames = list(NULL, "Age"))

y_full      <- Surv(clin_tr$survival_time, clin_tr$status)
time_full   <- as.numeric(clin_tr$survival_time)
status_full <- as.numeric(clin_tr$status)

log_msg(sprintf("   n = %d | events = %d | miRNAs = %d | clinical = 1",
                nrow(X_full), sum(status_full), ncol(X_full)))

# ======================================================================
# 6) Univariate Cox screening function
# ======================================================================
safe_univ_cox <- function(x, time, status) {
  ok <- is.finite(x) & is.finite(time) & is.finite(status)
  x <- x[ok]; time <- time[ok]; status <- status[ok]

  if (length(x) < 20) return(NA_real_)
  if (length(unique(x)) < 2) return(NA_real_)
  if (sum(status) < 5) return(NA_real_)

  fit <- tryCatch(
    coxph(Surv(time, status) ~ x, ties = "breslow"),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      tryCatch(coxph(Surv(time, status) ~ x, ties = "breslow"),
               error = function(e) NULL))
  )
  if (is.null(fit)) return(NA_real_)

  sm <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(sm) || nrow(sm$coefficients) < 1) return(NA_real_)

  p <- sm$coefficients[1, "Pr(>|z|)"]
  if (!is.finite(p)) NA_real_ else p
}

# ======================================================================
# 7) Scale train/test
# ======================================================================
scale_train_test <- function(train_x, test_x) {
  center    <- apply(train_x, 2, mean, na.rm = TRUE)
  scale_val <- apply(train_x, 2, sd,   na.rm = TRUE)
  scale_val[!is.finite(scale_val) | scale_val == 0] <- 1

  train_scaled <- sweep(train_x, 2, center, "-")
  train_scaled <- sweep(train_scaled, 2, scale_val, "/")

  test_scaled <- sweep(test_x, 2, center, "-")
  test_scaled <- sweep(test_scaled, 2, scale_val, "/")

  list(train = train_scaled, test = test_scaled,
       center = center, scale = scale_val)
}

# ======================================================================
# 8) One split
# ======================================================================
run_one_seed <- function(seed) {
  tryCatch({
    set.seed(seed)

    n <- nrow(X_full)
    train_n   <- floor(n * train_fraction)
    train_idx <- sample(seq_len(n), size = train_n, replace = FALSE)
    test_idx  <- setdiff(seq_len(n), train_idx)

    X_train <- X_full[train_idx, , drop = FALSE]
    X_test  <- X_full[test_idx,  , drop = FALSE]

    C_train <- clin_cov_matrix[train_idx, , drop = FALSE]
    C_test  <- clin_cov_matrix[test_idx,  , drop = FALSE]

    time_train   <- time_full[train_idx]
    status_train <- status_full[train_idx]
    time_test    <- time_full[test_idx]
    status_test  <- status_full[test_idx]

    # --- Univariate screening on train only ---
    pvals <- sapply(seq_len(ncol(X_train)), function(j)
      safe_univ_cox(X_train[, j], time_train, status_train))
    names(pvals) <- colnames(X_train)

    screened <- names(pvals)[is.finite(pvals) & pvals < screen_p]

    if (length(screened) == 0) {
      return(list(status = "failed", seed = seed,
                  message = "No miRNA passed univariate screening"))
    }

    Xs_train <- X_train[, screened, drop = FALSE]
    Xs_test  <- X_test[,  screened, drop = FALSE]

    # --- Scaling from train ---
    sc_mirna <- scale_train_test(Xs_train, Xs_test)
    sc_clin  <- scale_train_test(C_train, C_test)

    Z_train <- cbind(sc_mirna$train, sc_clin$train)
    Z_test  <- cbind(sc_mirna$test,  sc_clin$test)

    # --- Elastic Net Cox ---
    y_train <- Surv(time_train, status_train)

    cvfit <- tryCatch(
      cv.glmnet(x = as.matrix(Z_train), y = y_train, family = "cox",
                alpha = alpha, nfolds = 5, type.measure = "deviance",
                standardize = FALSE),
      error = function(e) NULL)
    if (is.null(cvfit)) {
      return(list(status = "failed", seed = seed,
                  message = "cv.glmnet failed"))
    }

    coef_mat <- as.matrix(coef(cvfit, s = "lambda.min"))
    coef_vec <- as.numeric(coef_mat[, 1])
    names(coef_vec) <- rownames(coef_mat)

    nonzero <- names(coef_vec)[is.finite(coef_vec) &
                                 abs(coef_vec) > 1e-10]

    selected_mirnas   <- intersect(nonzero, screened)
    selected_clinical <- intersect(nonzero, colnames(clin_cov_matrix))

    # --- rank miRNAs by absolute coefficient, capped at 10 ---
    if (length(selected_mirnas) > 0) {
      mirna_coef <- coef_vec[selected_mirnas]
      mirna_coef <- mirna_coef[order(abs(mirna_coef), decreasing = TRUE)]
      selected_mirnas <- names(mirna_coef)[
        seq_len(min(max_miRNAs, length(mirna_coef)))]
    }

    # --- Final refit on 10 miRNAs + clinical ---
    final_features <- unique(intersect(c(selected_mirnas,
                                          colnames(clin_cov_matrix)),
                                        colnames(Z_train)))
    Z_train_final <- Z_train[, final_features, drop = FALSE]
    Z_test_final  <- Z_test[,  final_features, drop = FALSE]

    final_fit <- glmnet(x = as.matrix(Z_train_final), y = y_train,
                        family = "cox", alpha = alpha,
                        lambda = cvfit$lambda.min,
                        standardize = FALSE)

    beta_final <- as.numeric(coef(final_fit)[, 1])
    names(beta_final) <- rownames(coef(final_fit))

    risk_train <- as.numeric(predict(final_fit,
                                      newx = as.matrix(Z_train_final),
                                      type = "link"))
    risk_test  <- as.numeric(predict(final_fit,
                                      newx = as.matrix(Z_test_final),
                                      type = "link"))

    c_train <- tryCatch(
      concordance(Surv(time_train, status_train) ~ risk_train,
                  reverse = TRUE)$concordance,
      error = function(e) NA_real_)
    c_test <- tryCatch(
      concordance(Surv(time_test, status_test) ~ risk_test,
                  reverse = TRUE)$concordance,
      error = function(e) NA_real_)

    # --- Comparison with the paper signature ---
    overlap_paper    <- intersect(selected_mirnas, original_signature)
    missing_paper    <- setdiff(original_signature, selected_mirnas)
    novel_paper      <- setdiff(selected_mirnas, original_signature)
    union_paper      <- length(union(selected_mirnas, original_signature))
    jaccard_paper    <- if (union_paper > 0) length(overlap_paper) / union_paper else NA

    # --- Comparison with the v2 signature ---
    overlap_v2       <- intersect(selected_mirnas, signature_v2_mirnas)
    missing_v2       <- setdiff(signature_v2_mirnas, selected_mirnas)
    novel_v2         <- setdiff(selected_mirnas, signature_v2_mirnas)
    union_v2         <- length(union(selected_mirnas, signature_v2_mirnas))
    jaccard_v2       <- if (union_v2 > 0) length(overlap_v2) / union_v2 else NA

    list(
      status = "success",
      seed = seed,
      n_train = length(train_idx),
      events_train = sum(status_train),
      n_validation = length(test_idx),
      events_validation = sum(status_test),
      n_screened = length(screened),
      n_selected = length(selected_mirnas),
      selected_mirnas = selected_mirnas,
      selected_clinical = selected_clinical,

      overlap_paper = overlap_paper,
      missing_paper = missing_paper,
      novel_paper   = novel_paper,
      overlap_n_paper = length(overlap_paper),
      jaccard_paper   = jaccard_paper,

      overlap_v2 = overlap_v2,
      missing_v2 = missing_v2,
      novel_v2   = novel_v2,
      overlap_n_v2 = length(overlap_v2),
      jaccard_v2   = jaccard_v2,

      c_train = c_train,
      c_validation = c_test,
      train_validation_gap = c_train - c_test,
      coefficients = beta_final
    )

  }, error = function(e) {
    list(status = "failed", seed = seed,
         message = conditionMessage(e))
  })
}

# ======================================================================
# 9) Run 50 seeds
# ======================================================================
log_msg("")
log_msg("============================================================")
log_msg(sprintf("starting %d seeds | max miRNAs = %d", n_seeds, max_miRNAs))
log_msg("============================================================")

results <- vector("list", length(seeds))

for (i in seq_along(seeds)) {
  seed <- seeds[i]
  cat(sprintf("[%02d/%02d] Seed %d ... ", i, length(seeds), seed))
  flush.console()

  res <- run_one_seed(seed)
  results[[i]] <- res

  if (isTRUE(res$status == "success")) {
    cat(sprintf("OK | C_val=%.3f | selected=%d | ov_paper=%d | ov_v2=%d\n",
                res$c_validation, res$n_selected,
                res$overlap_n_paper, res$overlap_n_v2))
  } else {
    cat("FAILED: ", res$message, "\n")
  }

  saveRDS(results, file.path(out_dir, "results_progress.rds"))
}

# ======================================================================
# 10) Successful runs
# ======================================================================
successful <- Filter(function(x) isTRUE(x$status == "success"), results)

log_msg("")
log_msg("============================================================")
log_msg(sprintf("successful: %d / %d", length(successful), length(results)))
log_msg("============================================================")

if (length(successful) == 0) stop("no seed was successful.")

# ======================================================================
# 11) Summary table
# ======================================================================
summary_df <- do.call(rbind, lapply(successful, function(x) {
  data.frame(
    seed              = x$seed,
    n_train           = x$n_train,
    events_train      = x$events_train,
    n_validation      = x$n_validation,
    events_validation = x$events_validation,
    n_screened        = x$n_screened,
    n_selected        = x$n_selected,
    overlap_paper     = x$overlap_n_paper,
    jaccard_paper     = x$jaccard_paper,
    overlap_v2        = x$overlap_n_v2,
    jaccard_v2        = x$jaccard_v2,
    c_train           = x$c_train,
    c_validation      = x$c_validation,
    train_validation_gap = x$train_validation_gap,
    stringsAsFactors  = FALSE)
}))

write.csv(summary_df,
          file.path(out_dir, "multi_seed_summary.csv"),
          row.names = FALSE)

# ======================================================================
# 12) miRNA list per seed
# ======================================================================
seed_mirna_df <- do.call(rbind, lapply(successful, function(x) {
  data.frame(seed = x$seed,
             rank = seq_along(x$selected_mirnas),
             miRNA = x$selected_mirnas,
             stringsAsFactors = FALSE)
}))

write.csv(seed_mirna_df,
          file.path(out_dir, "selected_miRNAs_every_seed.csv"),
          row.names = FALSE)

# ======================================================================
# 13) Wide table
# ======================================================================
seed_wide_df <- do.call(rbind, lapply(successful, function(x) {
  vals <- rep(NA_character_, max_miRNAs)
  if (length(x$selected_mirnas) > 0) {
    n_fill <- min(max_miRNAs, length(x$selected_mirnas))
    vals[seq_len(n_fill)] <- x$selected_mirnas[seq_len(n_fill)]
  }
  out <- as.data.frame(as.list(vals), stringsAsFactors = FALSE)
  colnames(out) <- paste0("miRNA_", seq_len(max_miRNAs))
  out$seed <- x$seed
  out[, c("seed", paste0("miRNA_", seq_len(max_miRNAs)))]
}))

write.csv(seed_wide_df,
          file.path(out_dir, "selected_miRNAs_wide_by_seed.csv"),
          row.names = FALSE)

# ======================================================================
# 14) Selection frequency
# ======================================================================
all_selected <- unlist(lapply(successful, function(x) x$selected_mirnas))
freq_tbl <- sort(table(all_selected), decreasing = TRUE)

frequency_df <- data.frame(
  miRNA     = names(freq_tbl),
  n_seeds   = as.integer(freq_tbl),
  frequency = as.numeric(freq_tbl) / length(successful),
  stringsAsFactors = FALSE)

write.csv(frequency_df,
          file.path(out_dir, "miRNA_selection_frequency.csv"),
          row.names = FALSE)

# ======================================================================
# 15) Comparison with the paper signature
# ======================================================================
comparison_df <- do.call(rbind, lapply(successful, function(x) {
  data.frame(
    seed              = x$seed,
    overlap_paper_n   = x$overlap_n_paper,
    jaccard_paper     = x$jaccard_paper,
    overlap_paper     = paste(x$overlap_paper, collapse = "; "),
    missing_paper     = paste(x$missing_paper, collapse = "; "),
    novel_paper       = paste(x$novel_paper, collapse = "; "),
    overlap_v2_n      = x$overlap_n_v2,
    jaccard_v2        = x$jaccard_v2,
    overlap_v2        = paste(x$overlap_v2, collapse = "; "),
    missing_v2        = paste(x$missing_v2, collapse = "; "),
    novel_v2          = paste(x$novel_v2, collapse = "; "),
    stringsAsFactors  = FALSE)
}))

write.csv(comparison_df,
          file.path(out_dir, "comparison_with_signatures.csv"),
          row.names = FALSE)

# ======================================================================
# 16) Frequency of the paper-signature miRNAs
# ======================================================================
original_frequency_paper <- data.frame(
  miRNA = original_signature,
  n_seeds = sapply(original_signature, function(m) {
    sum(sapply(successful, function(x) m %in% x$selected_mirnas))
  }),
  stringsAsFactors = FALSE)
original_frequency_paper$frequency <-
  original_frequency_paper$n_seeds / length(successful)
original_frequency_paper <- original_frequency_paper[
  order(original_frequency_paper$frequency, decreasing = TRUE), ]

write.csv(original_frequency_paper,
          file.path(out_dir, "original_signature_frequency.csv"),
          row.names = FALSE)

# ======================================================================
# 17) Frequency of the v2 miRNAs
# ======================================================================
original_frequency_v2 <- data.frame(
  miRNA = signature_v2_mirnas,
  n_seeds = sapply(signature_v2_mirnas, function(m) {
    sum(sapply(successful, function(x) m %in% x$selected_mirnas))
  }),
  stringsAsFactors = FALSE)
original_frequency_v2$frequency <-
  original_frequency_v2$n_seeds / length(successful)
original_frequency_v2 <- original_frequency_v2[
  order(original_frequency_v2$frequency, decreasing = TRUE), ]

write.csv(original_frequency_v2,
          file.path(out_dir, "signature_v2_frequency.csv"),
          row.names = FALSE)

# ======================================================================
# 18) Overall summary
# ======================================================================
overall_summary <- data.frame(
  successful_seeds      = length(successful),
  failed_seeds          = length(results) - length(successful),
  mean_validation_C     = mean(summary_df$c_validation, na.rm = TRUE),
  sd_validation_C       = sd(summary_df$c_validation,   na.rm = TRUE),
  median_validation_C   = median(summary_df$c_validation, na.rm = TRUE),
  min_validation_C      = min(summary_df$c_validation, na.rm = TRUE),
  max_validation_C      = max(summary_df$c_validation, na.rm = TRUE),
  mean_train_C          = mean(summary_df$c_train, na.rm = TRUE),
  mean_train_validation_gap = mean(summary_df$train_validation_gap, na.rm = TRUE),
  mean_selected_miRNAs  = mean(summary_df$n_selected, na.rm = TRUE),
  mean_overlap_paper    = mean(summary_df$overlap_paper, na.rm = TRUE),
  mean_Jaccard_paper    = mean(summary_df$jaccard_paper, na.rm = TRUE),
  mean_overlap_v2       = mean(summary_df$overlap_v2, na.rm = TRUE),
  mean_Jaccard_v2       = mean(summary_df$jaccard_v2, na.rm = TRUE)
)

write.csv(overall_summary,
          file.path(out_dir, "OVERALL_SUMMARY.csv"),
          row.names = FALSE)

# ======================================================================
# 19) Print the results
# ======================================================================
log_msg("")
log_msg("============================================================")
log_msg("final results")
log_msg("============================================================")
print(overall_summary)

log_msg("")
log_msg("Top 30 miRNAs by frequency:")
print(head(frequency_df, 30))

log_msg("")
log_msg("Frequency of the paper-signature miRNAs:")
print(original_frequency_paper)

log_msg("")
log_msg("Frequency of the v2 signature miRNAs:")
print(original_frequency_v2)

# ======================================================================
# 20) Save the full object
# ======================================================================
saveRDS(list(
  results             = results,
  successful          = successful,
  summary             = summary_df,
  seed_mirna          = seed_mirna_df,
  seed_wide           = seed_wide_df,
  selection_frequency = frequency_df,
  original_frequency_paper = original_frequency_paper,
  original_frequency_v2    = original_frequency_v2,
  comparison          = comparison_df,
  overall_summary     = overall_summary,
  original_signature  = original_signature,
  signature_v2_mirnas = signature_v2_mirnas
), file.path(out_dir, "MULTI_SEED_COMPLETE.rds"))

log_msg("")
log_msg("outputs in: ", out_dir)
log_msg("STEP 18 (Signature v2) COMPLETED")
