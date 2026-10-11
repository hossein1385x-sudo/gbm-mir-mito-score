#!/usr/bin/env Rscript
# =============================================================================
# STEP 10 (Signature v2) — Performance Evaluation
# -----------------------------------------------------------------------------
# Performance evaluation of the v2 signature (6 miRNA + Age) on the test set
#
# Includes:
#   - C-index (full + miRNA-only) with bootstrap CI
#   - HR univariate (full + miRNA-only)
#   - HR adjusted for Age (age-independent analysis)
#   - Time-dependent AUC (6mo / 1yr / 2yr)
#   - Kaplan-Meier (median + quartile)
#
# Outputs are saved in signature_v2/.
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
# 1) Load the v2 model and rebuild it on train
# ======================================================================
log_msg("============================================================")
log_msg("STEP 10 (Signature v2) — Performance Evaluation")
log_msg("============================================================")
log_msg("loading ...")

isoform_data  <- read.csv(file.path(output_dir, "mirna_isoform_rpm_matrix.csv"),
                          row.names = 1, check.names = FALSE)
clinical_data <- readRDS(file.path(output_dir, "raw_clinical_cdr_integrated.rds"))

# locked v2 signature
sig_path <- file.path(sig_dir, "signature_v2_locked_latest.rds")
if (!file.exists(sig_path)) stop("signature_v2_locked_latest.rds is missing")
signature_object <- readRDS(sig_path)

# list of features
sel_feats <- read.csv(file.path(sig_dir, "signature_features.csv"),
                      stringsAsFactors = FALSE)
mir_feats <- sel_feats$Feature[grepl("^hsa-", sel_feats$Feature)]
log_msg(sprintf("   miRNAs in the signature: %d | clinical: %s",
                length(mir_feats),
                paste(setdiff(sel_feats$Feature, mir_feats), collapse = ", ")))

# train + test ids
train_df  <- read.csv(file.path(output_dir, "train_set.csv"),
                      stringsAsFactors = FALSE)
test_df   <- read.csv(file.path(output_dir, "test_set_LOCKED.csv"),
                      stringsAsFactors = FALSE)
train_ids <- normalize_tcga_id(train_df$patient_barcode)
test_ids  <- normalize_tcga_id(test_df$patient_barcode)

# primary filter
mirna_raw <- readRDS(file.path(output_dir, "raw_mirna.rds"))
if (!is.data.table(mirna_raw)) mirna_raw <- as.data.table(mirna_raw)
mirna_meta <- unique(mirna_raw[, .(sample_barcode, is_primary_tumor)])
mirna_meta[, key := normalize_tcga_id(sample_barcode)]
primary <- unique(mirna_meta[is_primary_tumor == TRUE, key])
primary <- primary[!is.na(primary)]

iso_ids  <- normalize_tcga_id(rownames(isoform_data))
clin_ids <- normalize_tcga_id(as.character(clinical_data$patient_barcode))

prepare_cohort <- function(ids) {
  common <- Reduce(intersect, list(ids, iso_ids, clin_ids, primary))
  idx_iso  <- match(common, iso_ids)
  idx_clin <- match(common, clin_ids)

  iso <- isoform_data[idx_iso, , drop = FALSE]
  rownames(iso) <- common

  clin <- as.data.frame(clinical_data[idx_clin, , drop = FALSE])
  rownames(clin) <- common
  clin$status        <- clin$status_final
  clin$survival_time <- clin$survival_time_final

  keep <- !is.na(clin$status) & !is.na(clin$survival_time) &
          clin$survival_time >= 0
  clin <- clin[keep, , drop = FALSE]
  iso  <- iso[rownames(clin), , drop = FALSE]

  X <- log2(as.matrix(iso) + 1)
  keep_rows <- complete.cases(X)
  X <- X[keep_rows, , drop = FALSE]
  clin <- clin[keep_rows, , drop = FALSE]

  list(X = X, clin = clin)
}

train_data <- prepare_cohort(train_ids)
test_data  <- prepare_cohort(test_ids)
log_msg(sprintf("   train: %d | test: %d",
                nrow(train_data$clin), nrow(test_data$clin)))

y_train <- Surv(train_data$clin$survival_time, train_data$clin$status)
y_test  <- Surv(test_data$clin$survival_time,  test_data$clin$status)

# ======================================================================
# 2) Rebuild the v2 model on train (without using test)
# ======================================================================
X_train_all <- train_data$X
X_test_all  <- test_data$X

# --- miRNA scaling from train ---
mir_train <- X_train_all[, mir_feats, drop = FALSE]
mir_test  <- X_test_all[,  mir_feats, drop = FALSE]

mir_center <- colMeans(mir_train)
mir_scale  <- apply(mir_train, 2, sd)
mir_scale[mir_scale == 0] <- 1

X_tr_mir <- sweep(mir_train, 2, mir_center, "-")
X_tr_mir <- sweep(X_tr_mir, 2, mir_scale, "/")
X_te_mir <- sweep(mir_test,  2, mir_center, "-")
X_te_mir <- sweep(X_te_mir,  2, mir_scale, "/")

# --- Age ---
extract_age <- function(clin) {
  age_raw <- suppressWarnings(as.numeric(clin$age_at_diagnosis))
  if (all(is.na(age_raw))) age_raw <- suppressWarnings(as.numeric(clin$age))
  if (median(age_raw, na.rm = TRUE) > 200) age_raw / 365.25 else age_raw
}
age_train <- extract_age(train_data$clin)
age_test  <- extract_age(test_data$clin)

age_center <- mean(age_train, na.rm = TRUE)
age_scale  <- sd(age_train,   na.rm = TRUE)

age_train[is.na(age_train)] <- age_center
age_test[is.na(age_test)]   <- age_center

X_tr_age <- (age_train - age_center) / age_scale
X_te_age <- (age_test  - age_center) / age_scale

# --- combination ---
X_tr_full <- cbind(X_tr_mir, Age = X_tr_age)
X_te_full <- cbind(X_te_mir, Age = X_te_age)

set.seed(42)
cv_fit <- cv.glmnet(X_tr_full, y_train, family = "cox", alpha = 0.5, nfolds = 10)
final_model <- glmnet(X_tr_full, y_train, family = "cox", alpha = 0.5,
                      lambda = cv_fit$lambda.min)

coef_final <- as.matrix(coef(final_model))
log_msg("")
log_msg("v2 model coefficients (train):")
print(coef_final[coef_final[, 1] != 0, , drop = FALSE])

# ======================================================================
# 3) Apply to test
# ======================================================================
lp_full  <- as.numeric(predict(final_model, newx = X_te_full, type = "link"))

coefs_all <- coef_final[, 1]
mir_coefs <- coefs_all[mir_feats]
lp_mirna  <- as.numeric(X_te_mir %*% mir_coefs)

log_msg("")
log_msg(sprintf("lp_test full  : [%.3f, %.3f] | SD=%.3f",
                min(lp_full), max(lp_full), sd(lp_full)))
log_msg(sprintf("lp_test mirna : [%.3f, %.3f] | SD=%.3f",
                min(lp_mirna), max(lp_mirna), sd(lp_mirna)))

# ======================================================================
# 4) C-index (with bootstrap CI)
# ======================================================================
c_full <- as.numeric(concordance(y_test ~ lp_full,  reverse = TRUE)$concordance)
c_mir  <- as.numeric(concordance(y_test ~ lp_mirna, reverse = TRUE)$concordance)

set.seed(42)
boot_c <- function(score, B = 1000) {
  n <- length(score)
  vals <- numeric(B)
  for (b in 1:B) {
    idx <- sample(1:n, n, replace = TRUE)
    y_b <- Surv(test_data$clin$survival_time[idx], test_data$clin$status[idx])
    tryCatch({
      vals[b] <- as.numeric(
        concordance(y_b ~ score[idx], reverse = TRUE)$concordance)
    }, error = function(e) vals[b] <<- NA)
  }
  vals[!is.na(vals)]
}

boot_full <- boot_c(lp_full)
boot_mir  <- boot_c(lp_mirna)

c_full_ci <- quantile(boot_full, c(0.025, 0.975))
c_mir_ci  <- quantile(boot_mir,  c(0.025, 0.975))

log_msg("")
log_msg("============================================================")
log_msg("TEST C-INDEX")
log_msg("============================================================")
log_msg(sprintf("Full score:  C = %.4f (95%% CI: %.4f–%.4f)",
                c_full, c_full_ci[1], c_full_ci[2]))
log_msg(sprintf("miRNA-only:  C = %.4f (95%% CI: %.4f–%.4f)",
                c_mir, c_mir_ci[1], c_mir_ci[2]))

# ======================================================================
# 5) HR — univariate + age-adjusted
# ======================================================================
risk_z_full <- as.numeric(scale(lp_full))
risk_z_mir  <- as.numeric(scale(lp_mirna))

extract_hr <- function(cox_fit, coef_name) {
  s <- summary(cox_fit)
  b  <- s$coefficients[coef_name, "coef"]
  se <- s$coefficients[coef_name, "se(coef)"]
  p  <- s$coefficients[coef_name, "Pr(>|z|)"]
  c(hr = exp(b),
    lo = exp(b - 1.96 * se),
    hi = exp(b + 1.96 * se),
    p  = p)
}

# --- Full score univariate ---
cox_full <- coxph(y_test ~ risk_z_full)
hr_full  <- extract_hr(cox_full, "risk_z_full")

# --- miRNA-only univariate ---
cox_mir <- coxph(y_test ~ risk_z_mir)
hr_mir  <- extract_hr(cox_mir, "risk_z_mir")

# --- miRNA-only + Age (age-independent analysis) ---
cox_mir_adj <- coxph(y_test ~ risk_z_mir + age_test)
hr_mir_adj  <- extract_hr(cox_mir_adj, "risk_z_mir")

# --- Age alone ---
cox_age <- coxph(y_test ~ age_test)
hr_age  <- extract_hr(cox_age, "age_test")

log_msg("")
log_msg("============================================================")
log_msg("TEST HR per 1-SD")
log_msg("============================================================")
log_msg(sprintf("Full score (uni)           : HR=%.3f (95%% CI: %.3f–%.3f, p=%.4g)",
                hr_full["hr"], hr_full["lo"], hr_full["hi"], hr_full["p"]))
log_msg(sprintf("miRNA-only (uni)           : HR=%.3f (95%% CI: %.3f–%.3f, p=%.4g)",
                hr_mir["hr"], hr_mir["lo"], hr_mir["hi"], hr_mir["p"]))
log_msg(sprintf("miRNA-only (adj for Age)   : HR=%.3f (95%% CI: %.3f–%.3f, p=%.4g)",
                hr_mir_adj["hr"], hr_mir_adj["lo"], hr_mir_adj["hi"], hr_mir_adj["p"]))
log_msg(sprintf("Age alone (per year)       : HR=%.3f (95%% CI: %.3f–%.3f, p=%.4g)",
                hr_age["hr"], hr_age["lo"], hr_age["hi"], hr_age["p"]))

# ======================================================================
# 6) Time-dependent AUC
# ======================================================================
log_msg("")
log_msg("Time-dependent AUC ...")

if (!requireNamespace("timeROC", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  tryCatch(BiocManager::install("timeROC", update = FALSE, ask = FALSE),
           error = function(e) NULL)
}

times_days  <- c(180, 365, 730)
times_label <- c("6-month", "1-year", "2-year")

auc_full <- data.frame()
auc_mir  <- data.frame()

if (requireNamespace("timeROC", quietly = TRUE)) {
  suppressPackageStartupMessages(library(timeROC))

  roc_full <- timeROC(T = test_data$clin$survival_time,
                      delta = test_data$clin$status,
                      marker = lp_full, cause = 1,
                      times = times_days, iid = FALSE)

  auc_full <- data.frame(
    score_type = "full",
    horizon    = times_label,
    time_days  = times_days,
    AUC        = as.numeric(roc_full$AUC),
    n_at_risk  = sapply(times_days, function(t)
                        sum(test_data$clin$survival_time >= t)),
    n_events   = sapply(times_days, function(t)
                        sum(test_data$clin$survival_time <= t &
                            test_data$clin$status == 1))
  )

  roc_mir <- timeROC(T = test_data$clin$survival_time,
                     delta = test_data$clin$status,
                     marker = lp_mirna, cause = 1,
                     times = times_days, iid = FALSE)

  auc_mir <- data.frame(
    score_type = "mirna_only",
    horizon    = times_label,
    time_days  = times_days,
    AUC        = as.numeric(roc_mir$AUC),
    n_at_risk  = sapply(times_days, function(t)
                        sum(test_data$clin$survival_time >= t)),
    n_events   = sapply(times_days, function(t)
                        sum(test_data$clin$survival_time <= t &
                            test_data$clin$status == 1))
  )

  log_msg("Full score:")
  for (i in 1:3) log_msg(sprintf("   %-8s AUC = %.4f",
                                  times_label[i], auc_full$AUC[i]))
  log_msg("miRNA-only:")
  for (i in 1:3) log_msg(sprintf("   %-8s AUC = %.4f",
                                  times_label[i], auc_mir$AUC[i]))
} else {
  log_msg("timeROC could not be installed.")
}

# ======================================================================
# 7) Kaplan-Meier
# ======================================================================
log_msg("")
log_msg("Kaplan-Meier ...")

risk_group <- factor(ifelse(lp_full >= median(lp_full), "High", "Low"),
                     levels = c("Low", "High"))

km_diff <- survdiff(y_test ~ risk_group)
logrank_p <- 1 - pchisq(km_diff$chisq, df = 1)

risk_quartile <- cut(lp_full,
                     breaks = quantile(lp_full, probs = seq(0, 1, 0.25)),
                     labels = c("Q1", "Q2", "Q3", "Q4"),
                     include.lowest = TRUE)
km_quart <- survdiff(y_test ~ risk_quartile)
logrank_quart <- 1 - pchisq(km_quart$chisq, df = 3)

log_msg(sprintf("   Median split    log-rank p = %.4g", logrank_p))
log_msg(sprintf("   Quartile split  log-rank p = %.4g", logrank_quart))

# ======================================================================
# 8) Saving (everything in signature_v2/)
# ======================================================================
# Predictions
write.csv(data.frame(
  patient_id        = rownames(test_data$clin),
  risk_score_full   = lp_full,
  risk_score_mirna  = lp_mirna,
  age_years         = age_test,
  survival_time     = test_data$clin$survival_time,
  event             = test_data$clin$status,
  risk_group        = as.character(risk_group),
  risk_quartile     = as.character(risk_quartile)
), file.path(sig_dir, "signature_v2_test_predictions.csv"), row.names = FALSE)

# Performance table
perf_table <- data.frame(
  metric = c(
    "n_test", "n_events",
    "C_full", "C_full_CI_low", "C_full_CI_high",
    "C_mirna", "C_mirna_CI_low", "C_mirna_CI_high",
    "HR_full", "HR_full_CI_low", "HR_full_CI_high", "HR_full_p",
    "HR_mirna", "HR_mirna_CI_low", "HR_mirna_CI_high", "HR_mirna_p",
    "HR_mirna_adjAge", "HR_mirna_adjAge_CI_low",
    "HR_mirna_adjAge_CI_high", "HR_mirna_adjAge_p",
    "LogRank_median", "LogRank_quartile"
  ),
  value = c(
    nrow(test_data$clin), sum(test_data$clin$status == 1),
    c_full, c_full_ci[1], c_full_ci[2],
    c_mir, c_mir_ci[1], c_mir_ci[2],
    hr_full["hr"], hr_full["lo"], hr_full["hi"], hr_full["p"],
    hr_mir["hr"],  hr_mir["lo"],  hr_mir["hi"],  hr_mir["p"],
    hr_mir_adj["hr"], hr_mir_adj["lo"], hr_mir_adj["hi"], hr_mir_adj["p"],
    logrank_p, logrank_quart
  )
)
write.csv(perf_table,
          file.path(sig_dir, "signature_v2_test_performance.csv"),
          row.names = FALSE)

# AUC table
if (nrow(auc_full) > 0) {
  write.csv(rbind(auc_full, auc_mir),
            file.path(sig_dir, "signature_v2_time_dependent_auc.csv"),
            row.names = FALSE)
}

# Manifest
manifest <- list(
  step              = "10_signature_v2_evaluation",
  model             = "signature_v2_6mirna_plus_age",
  signature_hash    = signature_object$meta$locked_on,
  n_test            = nrow(test_data$clin),
  n_events_test     = sum(test_data$clin$status == 1),
  c_index_full      = c_full,
  c_index_full_ci   = c(c_full_ci[1], c_full_ci[2]),
  c_index_mirna     = c_mir,
  c_index_mirna_ci  = c(c_mir_ci[1], c_mir_ci[2]),
  hr_full           = unname(hr_full),
  hr_mirna          = unname(hr_mir),
  hr_mirna_adjAge   = unname(hr_mir_adj),
  hr_age_only       = unname(hr_age),
  logrank_median    = logrank_p,
  logrank_quartile  = logrank_quart,
  timestamp         = as.character(Sys.time()),
  r_version         = R.version.string
)
write_json(manifest,
           file.path(sig_dir, "step10_signature_v2_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE)

# ======================================================================
# 9) Final summary
# ======================================================================
log_msg("")
log_msg("============================================================")
log_msg("Final summary — Signature v2 (Test Set)")
log_msg("============================================================")
log_msg(sprintf("n = %d | events = %d",
                nrow(test_data$clin), sum(test_data$clin$status == 1)))
log_msg("")
log_msg("C-index:")
log_msg(sprintf("   Full        = %.4f (95%% CI: %.4f–%.4f)",
                c_full, c_full_ci[1], c_full_ci[2]))
log_msg(sprintf("   miRNA-only  = %.4f (95%% CI: %.4f–%.4f)",
                c_mir, c_mir_ci[1], c_mir_ci[2]))
log_msg("")
log_msg("HR per 1-SD:")
log_msg(sprintf("   Full               = %.3f (p=%.4g)", hr_full["hr"], hr_full["p"]))
log_msg(sprintf("   miRNA-only         = %.3f (p=%.4g)", hr_mir["hr"],  hr_mir["p"]))
log_msg(sprintf("   miRNA (adj. Age)   = %.3f (p=%.4g)",
                hr_mir_adj["hr"], hr_mir_adj["p"]))
log_msg(sprintf("   Age only (per year)= %.3f (p=%.4g)", hr_age["hr"], hr_age["p"]))
if (nrow(auc_full) > 0) {
  log_msg("")
  log_msg("Time-dependent AUC:")
  for (i in 1:3) {
    log_msg(sprintf("   %-8s Full=%.4f | miRNA=%.4f",
                    times_label[i], auc_full$AUC[i], auc_mir$AUC[i]))
  }
}
log_msg("")
log_msg(sprintf("Log-rank: median p=%.4g | quartile p=%.4g",
                logrank_p, logrank_quart))
log_msg("")
log_msg("saved in: ", sig_dir)
log_msg("STEP 10 (Signature v2) COMPLETED")
