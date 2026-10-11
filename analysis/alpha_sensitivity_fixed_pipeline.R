#!/usr/bin/env Rscript
# ============================================================
#  Alpha sensitivity — fixed pipeline
# ============================================================
rm(list = ls())

source(here::here("config.R"))  # defines WORK_DIR (see config.R)

ROOT <- WORK_DIR
OUT  <- file.path(ROOT, "ALPHA_SENSITIVITY")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(glmnet); library(survival)
})

set.seed(42)

## ---------- 1) Load ----------
pp <- readRDS(file.path(ROOT, "train_full_preprocessed.rds"))
X  <- as.matrix(pp$X_full)                     # 152 × 2196
Y  <- pp$y_full
age_raw <- as.numeric(pp$clin_cov_matrix[, "Age"])

cat("Original X dim:", dim(X), " events:", sum(Y[, "status"]), "\n")

## ---------- 2) FIX 1 — remove zero-variance columns ----------
sd_vec <- apply(X, 2, sd, na.rm = TRUE)
keep_var <- sd_vec > 0
X <- X[, keep_var, drop = FALSE]
cat("After removing SD=0:", ncol(X), "features\n")

## ---------- 3) FIX 2 — Age in years, not days ----------
age_years <- age_raw / 365.25
age_years <- as.numeric(scale(age_years))   # z-score

## ---------- 4) FIX 3 — proper prefiltering ----------
# Use variance percentile (top 60% by variance) as a fast prefilter
var_X <- apply(X, 2, var, na.rm = TRUE)
thr   <- quantile(var_X, probs = 0.60)
keep_pre <- var_X >= thr
X_pre <- X[, keep_pre, drop = FALSE]
cat("After variance prefilter:", ncol(X_pre), "features\n")

## ---------- 5) FIX 4 — univariate Cox screen ----------
# Faster than sapply loop; also let coxph handle ties properly
uni_p <- numeric(ncol(X_pre))
for (j in seq_len(ncol(X_pre))) {
  fit <- tryCatch(
    coxph(Y ~ X_pre[, j], ties = "efron"),
    error = function(e) NULL
  )
  uni_p[j] <- if (is.null(fit)) NA_real_ else summary(fit)$coefficients[1, "Pr(>|z|)"]
}
names(uni_p) <- colnames(X_pre)

# strict 0.01, fallback 0.05
keep_strict <- names(uni_p)[!is.na(uni_p) & uni_p < 0.01]
if (length(keep_strict) >= 20) {
  keep <- keep_strict
} else {
  keep <- names(uni_p)[!is.na(uni_p) & uni_p < 0.05]
}
cat("After univariate screen:", length(keep), "features\n")

# Cap to top 40 by p-value to keep glmnet stable
if (length(keep) > 40) {
  keep <- names(sort(uni_p[keep]))[1:40]
  cat("Capped to top 40 by p-value\n")
}

X_final <- cbind(X_pre[, keep, drop = FALSE], Age = age_years)
X_final <- scale(X_final)                  # ensure everything on same scale
cat("Final feature space:", ncol(X_final), "features\n")

## ---------- 6) Alpha loop ----------
alphas <- c(0.25, 0.50, 0.75, 1.00)
res    <- list()

for (a in alphas) {
  cat(sprintf("\n--- alpha = %.2f ---\n", a))
  set.seed(42)

  cv <- cv.glmnet(
    x            = X_final,
    y            = Y,
    family       = "cox",
    alpha        = a,
    nfolds       = 5,
    type.measure = "C",          # C-index is valid for Cox
    standardize  = FALSE,
    cox.ties     = "efron"       # silence the warning
  )

  # --- use lambda.min ---
  beta_min <- as.numeric(coef(cv, s = "lambda.min"))
  nz_min   <- sum(beta_min != 0)

  # --- use lambda.1se ---
  beta_1se <- as.numeric(coef(cv, s = "lambda.1se"))
  nz_1se   <- sum(beta_1se != 0)

  cat(sprintf("  lambda.min  -> %d features\n", nz_min))
  cat(sprintf("  lambda.1se  -> %d features\n", nz_1se))

  # choose lambda.min as primary
  idx_cv <- which.min(abs(cv$lambda - cv$lambda.min))
  c_cv   <- cv$cvm[idx_cv]
  c_se   <- cv$cvsd[idx_cv]

  chosen <- if (nz_min > 0) beta_min else beta_1se
  nz <- sum(chosen != 0)
  nz_names <- colnames(X_final)[chosen != 0]

  res[[length(res) + 1]] <- data.frame(
    alpha          = a,
    n_features_min = nz_min,
    n_features_1se = nz_1se,
    features       = paste(nz_names, collapse = "; "),
    C_index_CV     = round(c_cv, 4),
    C_index_SE     = round(c_se, 4),
    lambda_min     = signif(cv$lambda.min, 4)
  )

  rm(cv); gc()
}

res_df <- do.call(rbind, res)
cat("\n=========================================================\n")
cat("SUMMARY\n")
cat("=========================================================\n")
print(res_df[, c("alpha","n_features_min","n_features_1se",
                 "C_index_CV","C_index_SE")],
      row.names = FALSE)

write.csv(res_df,
          file.path(OUT, "alpha_sensitivity_summary.csv"),
          row.names = FALSE)
cat("\nSaved:", file.path(OUT, "alpha_sensitivity_summary.csv"), "\n")
