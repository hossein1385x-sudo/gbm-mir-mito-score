#!/usr/bin/env Rscript
# =========================================================
# lambda_min_vs_1se.R
# Compare lambda.min and lambda.1se on the test set
# Final version - with C-index correction
# =========================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
})

# ---------------------------------------------------------
# 1. Paths
# ---------------------------------------------------------
output_dir <- WORK_DIR
sig_dir    <- file.path(output_dir, "signature_v2")
dir.create(sig_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# 2. Helper functions
# ---------------------------------------------------------
normalize_tcga_id <- function(ids) {
  ids <- toupper(as.character(ids))
  ids <- gsub("\\.", "-", ids)
  m <- regmatches(ids, regexpr("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", ids))
  out <- rep(NA_character_, length(ids))
  out[nchar(m) > 0] <- m[nchar(m) > 0]
  out
}

scale_fit <- function(X) {
  list(center = colMeans(X, na.rm = TRUE),
       scale  = apply(X, 2, sd, na.rm = TRUE))
}

apply_scale <- function(X, s) {
  X <- sweep(X, 2, s$center, "-")
  sweep(X, 2, s$scale, "/")
}

safe_to_years <- function(v) {
  m <- median(v, na.rm = TRUE)
  if (is.na(m)) return(v)
  if (m > 200) v / 365.25 else v
}

# ---------------------------------------------------------
# 3. Evaluation function — with reverse=TRUE for the C-index
# ---------------------------------------------------------
evaluate_model <- function(lp_in, y_time_in, y_status_in, label = "") {
  ok <- is.finite(lp_in) & !is.na(y_time_in) & !is.na(y_status_in) & y_time_in >= 0
  cat("  [", label, "] n =", sum(ok), "/", length(ok), "\n")

  if (sum(ok) < 5) {
    return(list(C_index = NA, HR = NA, CI_low = NA,
                CI_high = NA, P = NA, n_used = sum(ok)))
  }

  t_use  <- as.numeric(y_time_in[ok])
  s_use  <- as.numeric(y_status_in[ok])
  lp_use <- as.numeric(lp_in[ok])

  if (sd(lp_use) == 0) {
    cat("   sd(lp) = 0\n")
    return(list(C_index = NA, HR = NA, CI_low = NA,
                CI_high = NA, P = NA, n_used = sum(ok)))
  }

  # C-index with reverse=TRUE
  y_surv <- survival::Surv(t_use, s_use)
  c_idx  <- survival::concordance(y_surv ~ lp_use, reverse = TRUE)$concordance

  # HR — with coxph.fit
  lp_scaled <- as.numeric(scale(lp_use))
  n_use     <- length(t_use)
  rn        <- as.character(seq_len(n_use))

  x_mat <- matrix(lp_scaled, ncol = 1)
  colnames(x_mat) <- "lpvar"
  rownames(x_mat) <- rn

  y_mat <- cbind(t_use, s_use)
  colnames(y_mat) <- c("time", "status")
  rownames(y_mat) <- rn

  fit <- survival::coxph.fit(
    x        = x_mat,
    y        = y_mat,
    strata   = NULL,
    offset   = NULL,
    init     = NULL,
    control  = survival::coxph.control(),
    weights  = NULL,
    method   = "efron",
    rownames = rn
  )

  coef_val <- fit$coefficients[1]
  se_val   <- sqrt(fit$var[1, 1])

  hr      <- exp(coef_val)
  ci_low  <- exp(coef_val - 1.96 * se_val)
  ci_high <- exp(coef_val + 1.96 * se_val)
  z       <- coef_val / se_val
  p       <- 2 * pnorm(-abs(z))

  list(C_index = c_idx, HR = hr, CI_low = ci_low,
       CI_high = ci_high, P = p, n_used = sum(ok))
}

# ---------------------------------------------------------
# 4. Load the data
# ---------------------------------------------------------
cat("=== loading the data ===\n")

isoform_data  <- read.csv(file.path(output_dir, "mirna_isoform_rpm_matrix.csv"),
                          row.names = 1, check.names = FALSE)
clinical_data <- readRDS(file.path(output_dir, "raw_clinical_cdr_integrated.rds"))
train_df      <- read.csv(file.path(output_dir, "train_set.csv"),
                          stringsAsFactors = FALSE)
test_df       <- read.csv(file.path(output_dir, "test_set_LOCKED.csv"),
                          stringsAsFactors = FALSE)
mirna_raw     <- readRDS(file.path(output_dir, "raw_mirna.rds"))

test_df <- test_df[!is.na(test_df$vital_status), ]

# ---------------------------------------------------------
# 5. Normalise the barcodes
# ---------------------------------------------------------
train_ids <- normalize_tcga_id(train_df$patient_barcode)
test_ids  <- normalize_tcga_id(test_df$patient_barcode)
iso_ids   <- normalize_tcga_id(rownames(isoform_data))
clin_ids  <- normalize_tcga_id(as.character(clinical_data$patient_barcode))

mirna_meta <- unique(as.data.frame(mirna_raw)[, c("sample_barcode","is_primary_tumor")])
mirna_meta$key <- normalize_tcga_id(mirna_meta$sample_barcode)
primary <- unique(mirna_meta$key[mirna_meta$is_primary_tumor == TRUE])
primary <- primary[!is.na(primary)]

cat("train_ids:", length(train_ids), "\n")
cat("test_ids :", length(test_ids),  "\n")
cat("iso_ids  :", length(iso_ids),   "\n")
cat("clin_ids :", length(clin_ids),  "\n")
cat("primary  :", length(primary),   "\n\n")

# ---------------------------------------------------------
# 6. build_set
# ---------------------------------------------------------
build_set <- function(ids, label = "") {
  common <- Reduce(intersect, list(ids, iso_ids, clin_ids, primary))

  X    <- log2(as.matrix(isoform_data[match(common, iso_ids), , drop = FALSE]) + 1)
  clin <- as.data.frame(clinical_data[match(common, clin_ids), , drop = FALSE])

  rownames(X)    <- common
  rownames(clin) <- common

  clin$status        <- clin$status_final
  clin$survival_time <- clin$survival_time_final

  keep <- !is.na(clin$status) & !is.na(clin$survival_time) & clin$survival_time >= 0
  X    <- X[keep, , drop = FALSE]
  clin <- clin[keep, , drop = FALSE]

  X    <- X[complete.cases(X), , drop = FALSE]
  clin <- clin[rownames(X), , drop = FALSE]

  cat("[", label, "] n =", nrow(X),
      "| events =", sum(clin$status == 1), "\n")

  list(X = X, clin = clin)
}

train <- build_set(train_ids, "train")
test  <- build_set(test_ids,  "test")
cat("\n")

# ---------------------------------------------------------
# 7. Pre-filter
# ---------------------------------------------------------
n_train   <- nrow(train$X)
keep_expr <- colSums(train$X > 1, na.rm = TRUE) >= 0.5 * n_train
X_tr      <- train$X[, keep_expr, drop = FALSE]

vars     <- apply(X_tr, 2, var, na.rm = TRUE)
keep_var <- vars > quantile(vars, 0.2, na.rm = TRUE)
X_tr     <- X_tr[, keep_var, drop = FALSE]
X_te     <- test$X[, colnames(X_tr), drop = FALSE]

cat("miRNA after prefilter:", ncol(X_tr), "\n\n")

# ---------------------------------------------------------
# 8. Age
# ---------------------------------------------------------
age_all    <- safe_to_years(suppressWarnings(as.numeric(clinical_data$age_at_diagnosis)))
median_age <- median(age_all, na.rm = TRUE)

age_tr <- safe_to_years(suppressWarnings(as.numeric(train$clin$age_at_diagnosis)))
age_tr[is.na(age_tr)] <- median_age

age_te <- safe_to_years(suppressWarnings(as.numeric(test$clin$age_at_diagnosis)))
age_te[is.na(age_te)] <- median_age

cat("Median age:", round(median_age, 1), "years\n")
cat("age_tr range:", round(range(age_tr), 1), "\n")
cat("age_te range:", round(range(age_te), 1), "\n\n")

clin_tr <- matrix(age_tr, ncol = 1, dimnames = list(NULL, "Age"))
clin_te <- matrix(age_te, ncol = 1, dimnames = list(NULL, "Age"))

y_tr <- Surv(train$clin$survival_time, train$clin$status)
y_te <- Surv(test$clin$survival_time,  test$clin$status)

cat("y_tr — events:", sum(y_tr[, "status"]), "/", nrow(y_tr), "\n")
cat("y_te — events:", sum(y_te[, "status"]), "/", nrow(y_te), "\n\n")

# ---------------------------------------------------------
# 9. Univariate screen
# ---------------------------------------------------------
p_vals <- apply(X_tr, 2, function(x) {
  fit <- tryCatch(coxph(y_tr ~ x), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  summary(fit)$coefficients[1, 5]
})

sel <- which(!is.na(p_vals) & p_vals < 0.01)
if (length(sel) < 20) sel <- which(!is.na(p_vals) & p_vals < 0.05)
X_sel <- X_tr[, sel, drop = FALSE]

cat("Selected miRNAs after screening:", ncol(X_sel), "\n\n")

# ---------------------------------------------------------
# 10. Scaling
# ---------------------------------------------------------
mir_scaler  <- scale_fit(X_sel)
clin_scaler <- scale_fit(clin_tr)

X_tr_s <- cbind(apply_scale(X_sel, mir_scaler),
                apply_scale(clin_tr, clin_scaler))

X_te_sel <- X_te[, colnames(X_sel), drop = FALSE]

for (j in seq_len(ncol(X_te_sel))) {
  na_idx <- is.na(X_te_sel[, j])
  if (any(na_idx)) X_te_sel[na_idx, j] <- mir_scaler$center[j]
}
for (j in seq_len(ncol(clin_te))) {
  na_idx <- is.na(clin_te[, j])
  if (any(na_idx)) clin_te[na_idx, j] <- clin_scaler$center[j]
}

X_te_s <- cbind(apply_scale(X_te_sel, mir_scaler),
                apply_scale(clin_te, clin_scaler))
X_te_s <- X_te_s[, colnames(X_tr_s), drop = FALSE]

cat("colnames match:", identical(colnames(X_tr_s), colnames(X_te_s)), "\n")
cat("NA in X_tr_s:", sum(is.na(X_tr_s)), "\n")
cat("NA in X_te_s:", sum(is.na(X_te_s)), "\n\n")

# ---------------------------------------------------------
# 11. cv.glmnet
# ---------------------------------------------------------
set.seed(2024)
cvfit <- cv.glmnet(x = X_tr_s, y = y_tr,
                   family = "cox", alpha = 0.5,
                   standardize = FALSE)

cat("lambda.min =", round(cvfit$lambda.min, 5), "\n")
cat("lambda.1se =", round(cvfit$lambda.1se, 5), "\n\n")

model_min <- glmnet(x = X_tr_s, y = y_tr, family = "cox",
                    alpha = 0.5, lambda = cvfit$lambda.min,
                    standardize = FALSE)
model_1se <- glmnet(x = X_tr_s, y = y_tr, family = "cox",
                    alpha = 0.5, lambda = cvfit$lambda.1se,
                    standardize = FALSE)

nz_min <- sum(coef(model_min)[, 1] != 0)
nz_1se <- sum(coef(model_1se)[, 1] != 0)
cat("Non-zero coef — lambda.min:", nz_min, "\n")
cat("Non-zero coef — lambda.1se:", nz_1se, "\n\n")

# ---------------------------------------------------------
# 12. predict
# ---------------------------------------------------------
lp_min <- as.numeric(predict(model_min, newx = X_te_s, type = "link"))
lp_1se <- as.numeric(predict(model_1se, newx = X_te_s, type = "link"))

cat("=== lp check ===\n")
cat("lp_min — NA:", sum(is.na(lp_min)),
    "| Inf:", sum(is.infinite(lp_min)),
    "| sd:", round(sd(lp_min), 4), "\n")
cat("lp_1se — NA:", sum(is.na(lp_1se)),
    "| Inf:", sum(is.infinite(lp_1se)),
    "| sd:", round(sd(lp_1se), 4), "\n\n")

# ---------------------------------------------------------
# 13. Evaluation
# ---------------------------------------------------------
y_time_te   <- as.numeric(y_te[, "time"])
y_status_te <- as.numeric(y_te[, "status"])

cat("--- lambda.min ---\n")
res_min <- evaluate_model(lp_min, y_time_te, y_status_te, "lambda.min")

cat("--- lambda.1se ---\n")
res_1se <- evaluate_model(lp_1se, y_time_te, y_status_te, "lambda.1se")

# ---------------------------------------------------------
# 14. Comparison table
# ---------------------------------------------------------
comparison <- data.frame(
  Lambda     = c("lambda.min", "lambda.1se"),
  Value      = round(c(cvfit$lambda.min, cvfit$lambda.1se), 5),
  NonZero    = c(nz_min, nz_1se),
  N_used     = c(res_min$n_used, res_1se$n_used),
  C_index    = round(c(res_min$C_index, res_1se$C_index), 4),
  HR_per_SD  = round(c(res_min$HR, res_1se$HR), 3),
  HR_CI_low  = round(c(res_min$CI_low, res_1se$CI_low), 3),
  HR_CI_high = round(c(res_min$CI_high, res_1se$CI_high), 3),
  P_value    = signif(c(res_min$P, res_1se$P), 3)
)

cat("\n============================================================\n")
cat("COMPARISON: lambda.min vs lambda.1se\n")
cat("============================================================\n")
print(comparison)

# ---------------------------------------------------------
# 15. Save
# ---------------------------------------------------------
write.csv(comparison,
          file.path(sig_dir, "lambda_min_vs_1se_test.csv"),
          row.names = FALSE)

saveRDS(list(cvfit = cvfit,
             model_min = model_min,
             model_1se = model_1se,
             comparison = comparison,
             res_min = res_min,
             res_1se = res_1se),
        file.path(sig_dir, "lambda_comparison.rds"))

cat("\nsaved in:", sig_dir, "\n")
