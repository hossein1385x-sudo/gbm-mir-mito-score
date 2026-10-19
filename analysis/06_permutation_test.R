#!/usr/bin/env Rscript
# =============================================================================
# analysis/06_permutation_test.R
# Permutation test for the miR-Mito Score
#
# Two modes:
#   MODE = "score_vs_survival_test"
#     H0: the locked score is unrelated to survival in the test set
#     each permutation: test survival is shuffled, C-index recomputed
#     (requires: DataS2 only)
#
#   MODE = "full_pipeline"
#     H0: no association exists in the population
#     each permutation: y_train is shuffled, LASSO-Cox refitted with CV
#              -> applied to the real test set
#     (requires: test miRNA matrix - TEST_X_PATH)
#
# Output: figures/Fig_permutation_null.{pdf,png}
#         supplementary/Table_S_permutation.csv
# =============================================================================

suppressPackageStartupMessages({
  library(survival); library(glmnet)
  library(dplyr);   library(readr); library(ggplot2)
})

# ---------------------------------------------------------------- 0) Configuration
MODE           <- "score_vs_survival_test"   # or "full_pipeline"
N_PERM         <- 200                        # use 1000 on a server
USE_PARALLEL   <- TRUE                       # parallelizes full_pipeline
N_CORES        <- 4

# Only for MODE = "full_pipeline"
TEST_X_PATH    <- NULL   # e.g. file.path(WORK_DIR, "test_mirna_matrix.rds")
TEST_AGE_PATH  <- NULL   # the model contains Age, so test Age must be supplied

source(here::here("config.R"))   # defines WORK_DIR (see README)
OUT <- WORK_DIR
FIG  <- here::here("figures")
SUPP <- here::here("supplementary")
dir.create(FIG,  showWarnings = FALSE, recursive = TRUE)
dir.create(SUPP, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- 1) Locked signature
sig  <- readRDS(file.path(OUT, "signature_locked_latest.rds"))
scal <- sig$scaling
ALPHA <- sig$alpha
beta <- setNames(sig$coefficients$Coefficient, sig$coefficients$Feature)

# ---------------------------------------------------------------- 2) Test data
d2 <- readr::read_csv(file.path(OUT, "supplementary_final", "DataS2_risk_scores.csv"),
                      show_col_types = FALSE)
test_df <- d2 %>% filter(cohort == "TCGA_Test") %>%
  transmute(sample_id = patient_id,
            score     = as.numeric(risk_score),
            OS_time   = as.numeric(survival_time),
            OS_event  = as.integer(event)) %>%
  filter(!is.na(score), !is.na(OS_time), !is.na(OS_event), OS_time > 0)
cat("Test n =", nrow(test_df), " events =", sum(test_df$OS_event), "\n")

# ---------------------------------------------------------------- 3) C-index helper
cindex_from_score <- function(time, event, score) {
  ok <- is.finite(time) & is.finite(event) & is.finite(score) & time > 0
  time <- time[ok]; event <- event[ok]; score <- score[ok]
  if (length(time) < 10 || length(unique(event)) < 2) return(NA_real_)
  df <- data.frame(t = time, e = event, s = score)
  # Fixed direction (higher = more risk); refitting coxph on permuted data would re-estimate beta
  # and make C = max(C, 1-C), inflating the null mean (~0.53 instead of ~0.5)
  as.numeric(survival::concordance(survival::Surv(t, e) ~ s, data = df, reverse = TRUE)$concordance)
}

# =============================================================================
# Mode A: score_vs_survival_test
# =============================================================================
run_mode_A <- function() {
  c_real <- cindex_from_score(test_df$OS_time, test_df$OS_event, test_df$score)
  cat(sprintf("Real test C-index = %.4f\n\n", c_real))

  cat(sprintf("Running %d permutations (shuffling survival in test)...\n", N_PERM))
  c_perm <- numeric(N_PERM)
  t0 <- Sys.time()
  for (i in seq_len(N_PERM)) {
    set.seed(10000 + i)
    idx <- sample(seq_len(nrow(test_df)))
    c_perm[i] <- cindex_from_score(
      test_df$OS_time[idx], test_df$OS_event[idx], test_df$score)
    if (i %% 25 == 0)
      cat(sprintf("  perm %3d / %d  |  elapsed %.1f min\n",
                  i, N_PERM, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  c_perm <- c_perm[is.finite(c_perm)]

  list(c_real = c_real, c_perm = c_perm, mode = "score_vs_survival_test")
}

# =============================================================================
# Mode B: full_pipeline
# =============================================================================
run_mode_B <- function() {
  if (is.null(TEST_X_PATH) || !file.exists(TEST_X_PATH))
    stop("MODE='full_pipeline' requires TEST_X_PATH (miRNA matrix for test).")
  if (!is.null(TEST_AGE_PATH) && !file.exists(TEST_AGE_PATH))
    stop("TEST_AGE_PATH given but not found.")

  tr <- readRDS(file.path(OUT, "train_full_preprocessed.rds"))
  X  <- as.matrix(tr$X_full); y <- tr$y_full
  ts <- readr::read_csv(file.path(OUT, "train_set.csv"), show_col_types = FALSE)
  ts$Age_years <- ts$age_at_diagnosis / 365.25
  idx <- match(rownames(X), ts$patient_barcode)
  Age_train <- ts$Age_years[idx]
  if (any(is.na(Age_train))) Age_train[is.na(Age_train)] <- median(ts$Age_years, na.rm = TRUE)

  lock_mirnas <- sig$mirna_features
  Xte <- readRDS(TEST_X_PATH)               # matrix or data.frame with rownames = patient_id
  Xte <- as.matrix(Xte)
  # Align Xte rows with test_df (time/event) by patient_id instead of assuming the same order
  mm <- match(test_df$sample_id, rownames(Xte))
  if (anyNA(mm)) stop("Xte rownames do not cover all test patients (", sum(is.na(mm)), " missing)")
  Xte <- Xte[mm, , drop = FALSE]
  # Refit on ALL shared candidate miRNAs (not only the 2 locked ones), so that feature selection is also permuted;
  # otherwise the null is biased upwards because the 2 miRNAs were selected on the real y.
  feat_mirna <- intersect(colnames(X), colnames(Xte))
  if (!all(lock_mirnas %in% feat_mirna))
    stop("Locked miRNAs missing from train/test matrices: ",
         paste(setdiff(lock_mirnas, feat_mirna), collapse = ", "))
  feat_all   <- c(feat_mirna, "Age")

  if (!is.null(TEST_AGE_PATH)) {
    age_te_df <- readr::read_csv(TEST_AGE_PATH, show_col_types = FALSE)
    # Expected columns: patient_id and Age_years
    age_te <- age_te_df$Age_years[match(rownames(Xte), age_te_df$patient_id)]
    age_te[is.na(age_te)] <- median(Age_train, na.rm = TRUE)
  } else {
    age_te <- rep(median(Age_train, na.rm = TRUE), nrow(Xte))
  }

  # --- Real score on the test set ---
  zte1 <- (Xte[, lock_mirnas[1]] - scal$mirna$center[lock_mirnas[1]]) /
    scal$mirna$scale [lock_mirnas[1]]
  zte2 <- (Xte[, lock_mirnas[2]] - scal$mirna$center[lock_mirnas[2]]) /
    scal$mirna$scale [lock_mirnas[2]]
  zteA <- (age_te           - scal$clinical$center["Age"])      /
    scal$clinical$scale ["Age"]
  lp_test_real <- beta[lock_mirnas[1]]*zte1 + beta[lock_mirnas[2]]*zte2 + beta["Age"]*zteA

  c_real <- cindex_from_score(test_df$OS_time, test_df$OS_event, lp_test_real)
  cat(sprintf("Real test C-index = %.4f\n\n", c_real))

  # --- One permutation: shuffle y_train, refit LASSO-Cox on X_train, apply to test ---
  one_perm <- function(i) {
    set.seed(20000 + i)
    # shuffle train survival
    perm_idx <- sample(seq_len(nrow(X)))
    y_perm <- y[perm_idx, ]

    Xtr <- cbind(X[, feat_mirna, drop = FALSE],
                 Age = (Age_train - scal$clinical$center["Age"]) /
                   scal$clinical$scale ["Age"])

    fit <- tryCatch(
      glmnet::cv.glmnet(x = Xtr, y = y_perm,
                        family = "cox", alpha = ALPHA, nfolds = 10,
                        standardize = TRUE),
      error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)

    coefs <- as.numeric(coef(fit, s = "lambda.min"))
    names(coefs) <- c(colnames(Xtr), "dummy")[seq_along(coefs)]
    coefs <- coefs[colnames(Xtr)]   # same order as Xtr

    # apply to test
    Xte_all <- cbind(Xte[, feat_mirna, drop = FALSE], Age = zteA)
    lp <- as.numeric(as.matrix(Xte_all) %*% coefs)
    cindex_from_score(test_df$OS_time, test_df$OS_event, lp)
  }

  cat(sprintf("Running %d full-pipeline permutations%s...\n",
              N_PERM, if (USE_PARALLEL) sprintf(" (parallel, %d cores)", N_CORES) else ""))

  t0 <- Sys.time()
  if (USE_PARALLEL) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      cat("[warn] future.apply not installed — falling back to sequential.\n")
      USE_PARALLEL <- FALSE
    }
  }
  if (USE_PARALLEL) {
    future::plan(future::multisession, workers = N_CORES)
    c_perm <- future.apply::future_sapply(seq_len(N_PERM), one_perm, future.seed = TRUE)
    future::plan(future::sequential)
  } else {
    c_perm <- sapply(seq_len(N_PERM), one_perm)
  }
  c_perm <- as.numeric(c_perm); c_perm <- c_perm[is.finite(c_perm)]
  cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  list(c_real = c_real, c_perm = c_perm, mode = "full_pipeline")
}

# ---------------------------------------------------------------- 4) Run
res <- if (MODE == "score_vs_survival_test") run_mode_A() else run_mode_B()

# ---------------------------------------------------------------- 5) p-value
c_real <- res$c_real
c_perm <- res$c_perm
n      <- length(c_perm)
p_emp  <- (1 + sum(c_perm >= c_real)) / (1 + n)   # Phipson & Smyth corrected

cat(sprintf("\n=== Permutation test result ===\n"))
cat(sprintf("Mode            : %s\n", res$mode))
cat(sprintf("Real C-index    : %.4f\n", c_real))
cat(sprintf("N permutations  : %d\n", n))
cat(sprintf("Null mean       : %.4f\n", mean(c_perm)))
cat(sprintf("Null sd         : %.4f\n", sd(c_perm)))
cat(sprintf("Null quantiles  : 2.5%%=%.4f  median=%.4f  97.5%%=%.4f\n",
            quantile(c_perm, 0.025), median(c_perm), quantile(c_perm, 0.975)))
cat(sprintf("Empirical p     : %.4f  (one-sided, Phipson-corrected)\n", p_emp))

# ---------------------------------------------------------------- 6) Output
tab <- data.frame(
  metric = c("real_C", "n_perm", "null_mean", "null_sd",
             "null_q2.5", "null_median", "null_q97.5", "p_one_sided",
             "p_min_possible", "mode_full_pipeline(1)_or_score_only(0)"),   # record the mode and the minimum attainable p
  value  = round(c(c_real, n, mean(c_perm), sd(c_perm),
                   quantile(c_perm, 0.025), median(c_perm),
                   quantile(c_perm, 0.975), p_emp,
                   1 / (n + 1), as.numeric(res$mode == "full_pipeline")), 4))
readr::write_csv(tab, file.path(SUPP, "Table_S_permutation.csv"))

perm_df <- data.frame(C = c_perm)
p <- ggplot(perm_df, aes(x = C)) +
  geom_histogram(bins = 30, fill = "#B0BEC5", color = "white") +
  geom_vline(xintercept = c_real, color = "red", linewidth = 1) +
  annotate("text", x = c_real, y = Inf,
           label = sprintf("Observed C = %.3f\np = %.4f", c_real, p_emp),
           hjust = -0.05, vjust = 1.3, color = "red", size = 3.5) +
  labs(title = "Permutation null distribution of test C-index",
       subtitle = sprintf("%s | n_perm = %d | one-sided p", res$mode, n),
       x = "Test C-index (null)", y = "Count") +
  theme_bw()
ggsave(file.path(FIG, "Fig_permutation_null.pdf"), p, width = 6.5, height = 4.2)
ggsave(file.path(FIG, "Fig_permutation_null.png"), p, width = 6.5, height = 4.2, dpi = 300)

cat("\nDone.\n")
cat("  Table : supplementary/Table_S_permutation.csv\n")
cat("  Figure: figures/Fig_permutation_null.pdf / .png\n")
