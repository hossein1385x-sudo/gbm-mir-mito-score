#!/usr/bin/env Rscript
# =============================================================================
# Sensitivity Analysis — Elastic Net alpha (0.1 / 0.5 / 0.9)  — FIXED
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(glmnet)
  library(jsonlite)
})

# =============================================================================
# Settings
# =============================================================================
input_dir <- paste0(WORK_DIR, "/")
sig_dir   <- file.path(input_dir, "signature_v2")
supp_dir  <- file.path(input_dir, "supplementary_final")
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

norm_id <- function(x) substr(toupper(gsub("\\.", "-", as.character(x))), 1, 12)
log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

ALPHAS     <- c(0.1, 0.5, 0.9)
N_OUTER    <- 5
N_INNER    <- 5
N_REPEATS  <- 10
CONSENSUS  <- 0.60
EPV_MIN    <- 10
SEED_BASE  <- 2024
REF_ALPHA  <- 0.5

# =============================================================================
# Robust helper functions
# =============================================================================

# --- reliable concordance ---
safe_cindex <- function(y_surv, lp) {
  if (length(lp) < 5 || length(unique(lp)) < 2) return(NA_real_)
  res <- tryCatch(
    survival::concordance(y_surv ~ lp, reverse = TRUE)$concordance,
    error = function(e) NA_real_)
  if (is.na(res)) {
    res <- tryCatch(
      survival::survConcordance(y_surv ~ lp)$concordance,
      error = function(e) NA_real_)
  }
  as.numeric(res)
}

# --- glmnet cox with explicit settings ---
glmnet_cox <- function(x, y, alpha, lambda = NULL) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (is.null(lambda)) {
    glmnet(x, y, family = "cox", alpha = alpha, cox.ties = "breslow")
  } else {
    glmnet(x, y, family = "cox", alpha = alpha, lambda = lambda,
           cox.ties = "breslow")
  }
}

cv_glmnet_cox <- function(x, y, alpha, nfolds) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  cv.glmnet(x, y, family = "cox", alpha = alpha, nfolds = nfolds,
            cox.ties = "breslow")
}

# --- reliable predict ---
safe_predict <- function(model, newx) {
  newx <- as.matrix(newx)
  storage.mode(newx) <- "double"
  p <- tryCatch(
    as.numeric(predict(model, newx = newx, type = "link")),
    error = function(e) NULL)
  p
}

# =============================================================================
# STEP 0 — Load the data
# =============================================================================
log_msg("============================================================")
log_msg("STEP 0 — loading the data")
log_msg("============================================================")

iso <- read.csv(file.path(input_dir, "mirna_isoform_rpm_matrix.csv"),
                row.names = 1, check.names = FALSE)
if (mean(grepl("^TCGA-", rownames(iso))) > 0.5) {
  mirna_mat <- t(as.matrix(iso))
} else {
  mirna_mat <- as.matrix(iso)
}
storage.mode(mirna_mat) <- "numeric"
colnames(mirna_mat) <- norm_id(colnames(mirna_mat))

data_s1 <- fread(file.path(supp_dir, "DataS1_barcodes.csv"))
data_s1[, key := norm_id(patient_barcode_short)]
train_keys <- data_s1[split == "train", key]
test_keys  <- data_s1[split == "test",  key]

clin <- as.data.table(readRDS(file.path(input_dir, "raw_clinical_cdr_integrated.rds")))
clin[, key := norm_id(patient_barcode)]
clin[, status := as.integer(tolower(vital_status) %in%
                              c("dead", "deceased", "1"))]
clin[, time := ifelse(
  status == 1 & !is.na(days_to_death) & days_to_death > 0,
  as.numeric(days_to_death),
  as.numeric(days_to_last_follow_up)
)]

clin_ok <- clin[!is.na(time) & time > 0 & !is.na(status)]

train_common <- intersect(intersect(train_keys, colnames(mirna_mat)), clin_ok$key)
test_common  <- intersect(intersect(test_keys,  colnames(mirna_mat)), clin_ok$key)

log_msg(sprintf("  shared train: %d | shared test: %d",
                length(train_common), length(test_common)))

X_train_full <- t(mirna_mat[, train_common, drop = FALSE])
X_test_full  <- t(mirna_mat[, test_common,  drop = FALSE])

# Surv objects
clin_train <- clin_ok[match(train_common, key)]
clin_test  <- clin_ok[match(test_common, key)]
y_train_full <- Surv(clin_train$time, clin_train$status)
y_test_full  <- Surv(clin_test$time,  clin_test$status)

log_msg(sprintf("  X_train: %d × %d | X_test: %d × %d",
                nrow(X_train_full), ncol(X_train_full),
                nrow(X_test_full),  ncol(X_test_full)))

# =============================================================================
# STEP 1 — prefilter
# =============================================================================
log_msg("")
log_msg("STEP 1 — pre-filter")

keep_expr <- colSums(X_train_full > 1, na.rm = TRUE) >= 0.5 * nrow(X_train_full)
X_train_log <- log2(X_train_full[, keep_expr, drop = FALSE] + 1)

vars <- apply(X_train_log, 2, var, na.rm = TRUE)
keep_var <- vars > quantile(vars, 0.2, na.rm = TRUE)
keep_features <- colnames(X_train_log)[keep_var]
X_train <- X_train_log[, keep_features, drop = FALSE]

log_msg(sprintf("  final features: %d", length(keep_features)))

# =============================================================================
# STEP 2 — nested CV
# =============================================================================
run_nested_cv <- function(X, y, alpha, seed_base = SEED_BASE,
                          n_outer = N_OUTER, n_inner = N_INNER,
                          n_repeats = N_REPEATS) {

  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X)

  all_selected <- list()
  c_index_outer <- numeric(0)

  for (rep_i in seq_len(n_repeats)) {
    set.seed(seed_base + rep_i)

    foldid <- integer(n)
    for (ev in c(0, 1)) {
      idx <- which(y[, "status"] == ev)
      foldid[idx] <- sample(rep(seq_len(n_outer), length.out = length(idx)))
    }

    for (fold in seq_len(n_outer)) {
      test_idx  <- which(foldid == fold)
      train_idx <- setdiff(seq_len(n), test_idx)

      X_tr <- X[train_idx, , drop = FALSE]
      X_te <- X[test_idx,  , drop = FALSE]
      y_tr <- y[train_idx]
      y_te <- y[test_idx]

      # univariate screen
      univ_p <- apply(X_tr, 2, function(col) {
        if (sd(col, na.rm = TRUE) == 0) return(NA_real_)
        fit <- tryCatch(coxph(y_tr ~ col), error = function(e) NULL)
        if (is.null(fit)) return(NA_real_)
        summary(fit)$coefficients[1, "Pr(>|z|)"]
      })

      keep_uni <- which(!is.na(univ_p) & univ_p < 0.01)
      if (length(keep_uni) < 20)
        keep_uni <- which(!is.na(univ_p) & univ_p < 0.05)
      if (length(keep_uni) < 2) next

      sub_tr <- X_tr[, keep_uni, drop = FALSE]
      tr_mean <- colMeans(sub_tr)
      tr_sd   <- apply(sub_tr, 2, sd, na.rm = TRUE)
      valid   <- !is.na(tr_sd) & tr_sd > 0
      if (sum(valid) < 2) next
      keep_final <- keep_uni[valid]

      X_tr_s <- scale(X_tr[, keep_final, drop = FALSE],
                      center = tr_mean[valid], scale = tr_sd[valid])
      X_te_s <- scale(X_te[, keep_final, drop = FALSE],
                      center = tr_mean[valid], scale = tr_sd[valid])

      X_tr_s <- as.matrix(X_tr_s); storage.mode(X_tr_s) <- "double"
      X_te_s <- as.matrix(X_te_s); storage.mode(X_te_s) <- "double"
      colnames(X_tr_s) <- colnames(X_tr)[keep_final]
      colnames(X_te_s) <- colnames(X_te)[keep_final]

      n_ev_tr <- sum(y_tr[, "status"] == 1)
      if (n_ev_tr < n_inner) next

      cv_fit <- tryCatch(
        cv_glmnet_cox(X_tr_s, y_tr, alpha = alpha,
                      nfolds = min(n_inner, n_ev_tr)),
        error = function(e) NULL)
      if (is.null(cv_fit)) next

      fold_model <- tryCatch(
        glmnet_cox(X_tr_s, y_tr, alpha = alpha,
                   lambda = cv_fit$lambda.min),
        error = function(e) NULL)
      if (is.null(fold_model)) next

      coefs <- as.matrix(coef(fold_model))
      nz <- rownames(coefs)[coefs[, 1] != 0]
      if (length(nz) < 1) next

      all_selected[[length(all_selected) + 1]] <- nz

      # predict on test
      nz_valid <- intersect(nz, colnames(X_te_s))
      if (length(nz_valid) < 1) next
      X_te_aligned <- as.matrix(X_te_s[, nz_valid, drop = FALSE])
      storage.mode(X_te_aligned) <- "double"

      lp_te <- safe_predict(fold_model, X_te_aligned)
      if (is.null(lp_te)) next

      ci <- safe_cindex(y_te, lp_te)
      if (!is.na(ci)) c_index_outer <- c(c_index_outer, ci)
    }
  }

  if (length(all_selected) == 0) {
    return(list(features = character(0),
                c_index_outer = c_index_outer,
                n_valid_folds = 0))
  }

  all_feat <- unlist(all_selected)
  freq_tbl <- table(all_feat) / length(all_selected)
  consensus <- names(freq_tbl)[freq_tbl >= CONSENSUS]
  consensus <- consensus[order(-freq_tbl[consensus])]

  list(features = consensus,
       freq_table = freq_tbl,
       c_index_outer = c_index_outer,
       n_valid_folds = length(all_selected))
}

# =============================================================================
# STEP 3 — final fit + evaluation
# =============================================================================
fit_and_evaluate <- function(X_train, y_train, features,
                             X_test, y_test, alpha) {

  if (length(features) < 2) {
    return(list(features = features, c_train = NA_real_, c_test = NA_real_,
                hr_test = NA_real_, hr_low = NA_real_, hr_high = NA_real_,
                p_test = NA_real_, note = "too few features"))
  }

  # --- make sure the features are in X_train and X_test ---
  features <- intersect(features, colnames(X_train))
  features <- intersect(features, colnames(X_test))
  if (length(features) < 2) {
    return(list(features = features, c_train = NA_real_, c_test = NA_real_,
                hr_test = NA_real_, hr_low = NA_real_, hr_high = NA_real_,
                p_test = NA_real_, note = "features not in matrices"))
  }

  sub_tr <- as.matrix(X_train[, features, drop = FALSE])
  storage.mode(sub_tr) <- "double"

  tr_mean <- colMeans(sub_tr)
  tr_sd   <- apply(sub_tr, 2, sd, na.rm = TRUE)
  tr_sd[!is.finite(tr_sd) | tr_sd == 0] <- 1

  X_tr_s <- scale(sub_tr, center = tr_mean, scale = tr_sd)
  X_tr_s <- as.matrix(X_tr_s); storage.mode(X_tr_s) <- "double"
  colnames(X_tr_s) <- features

  n_ev <- sum(y_train[, "status"] == 1)
  cv_fit <- cv_glmnet_cox(X_tr_s, y_train, alpha = alpha,
                          nfolds = min(10, n_ev))
  final_model <- glmnet_cox(X_tr_s, y_train, alpha = alpha,
                            lambda = cv_fit$lambda.min)

  coefs <- as.matrix(coef(final_model))
  nz <- rownames(coefs)[coefs[, 1] != 0]
  if (length(nz) < 1) {
    return(list(features = features, c_train = NA_real_, c_test = NA_real_,
                hr_test = NA_real_, hr_low = NA_real_, hr_high = NA_real_,
                p_test = NA_real_, note = "all coef zero"))
  }

  # C-index train
  X_tr_nz <- as.matrix(X_tr_s[, nz, drop = FALSE])
  lp_tr <- safe_predict(final_model, X_tr_nz)
  c_train <- if (is.null(lp_tr)) NA_real_ else safe_cindex(y_train, lp_tr)

  # apply on test
  nz_valid <- intersect(nz, colnames(X_test))
  if (length(nz_valid) < 1) {
    return(list(features = features, c_train = c_train, c_test = NA_real_,
                hr_test = NA_real_, hr_low = NA_real_, hr_high = NA_real_,
                p_test = NA_real_, note = "features not in X_test"))
  }

  sub_te <- as.matrix(X_test[, nz_valid, drop = FALSE])
  storage.mode(sub_te) <- "double"
  X_te_s <- scale(sub_te, center = tr_mean[nz_valid], scale = tr_sd[nz_valid])
  X_te_s <- as.matrix(X_te_s); storage.mode(X_te_s) <- "double"
  colnames(X_te_s) <- nz_valid

  lp_te <- safe_predict(final_model, X_te_s)
  if (is.null(lp_te)) {
    return(list(features = features, c_train = c_train, c_test = NA_real_,
                hr_test = NA_real_, hr_low = NA_real_, hr_high = NA_real_,
                p_test = NA_real_, note = "predict failed"))
  }

  c_test <- safe_cindex(y_test, lp_te)

  # HR per SD
  hr <- hr_low <- hr_high <- p_val <- NA_real_
  fit_cox <- tryCatch(
    coxph(y_test ~ scale(lp_te)),
    error = function(e) NULL)
  if (!is.null(fit_cox)) {
    s <- summary(fit_cox)
    b  <- s$coefficients[1, "coef"]
    se <- s$coefficients[1, "se(coef)"]
    hr     <- exp(b)
    hr_low <- exp(b - 1.96 * se)
    hr_high<- exp(b + 1.96 * se)
    p_val  <- s$coefficients[1, "Pr(>|z|)"]
  }

  list(features = features, nonzero_features = nz_valid,
       c_train = c_train, c_test = c_test,
       hr_test = hr, hr_low = hr_low, hr_high = hr_high, p_test = p_val,
       note = "OK")
}

# =============================================================================
# STEP 4 — run for each alpha
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("STEP 4 — running the pipeline for each alpha")
log_msg("============================================================")

results <- list()
all_selected_dt <- data.table()

for (a in ALPHAS) {
  log_msg("")
  log_msg(sprintf("▶ alpha = %.1f", a))

  t0 <- Sys.time()
  cv_res <- tryCatch(run_nested_cv(X_train, y_train_full, alpha = a),
                     error = function(e) {
                       log_msg("   error in nested CV: ", conditionMessage(e))
                       NULL })
  if (is.null(cv_res)) next

  t_cv <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_msg(sprintf("   nested CV: %d folds | %d-feature consensus | %.1f minutes",
                  cv_res$n_valid_folds, length(cv_res$features), t_cv))

  if (length(cv_res$c_index_outer) > 0) {
    log_msg(sprintf("   C-index outer: %.3f ± %.3f",
                    mean(cv_res$c_index_outer, na.rm = TRUE),
                    sd(cv_res$c_index_outer, na.rm = TRUE)))
  }

  n_events_train <- sum(y_train_full[, "status"] == 1)
  max_features <- max(1, floor(n_events_train / EPV_MIN) - 1)

  features_final <- cv_res$features
  if (length(features_final) > max_features) {
    features_final <- head(features_final, max_features)
    log_msg(sprintf("   EPV cap → %d feature", length(features_final)))
  }
  if (length(features_final) < 2) {
    log_msg("   not enough features")
    next
  }

  ev <- tryCatch(
    fit_and_evaluate(X_train, y_train_full, features_final,
                     X_test_full, y_test_full, alpha = a),
    error = function(e) {
      log_msg("   error in fit/evaluate: ", conditionMessage(e))
      NULL })
  if (is.null(ev)) next

  log_msg(sprintf("   Final: C_train = %.3f | C_test = %.3f | HR = %.3f | P = %.4g",
                  ev$c_train, ev$c_test, ev$hr_test, ev$p_test))

  results[[as.character(a)]] <- list(
    alpha = a,
    n_features_consensus = length(cv_res$features),
    n_features_final = length(ev$features),
    consensus_features = cv_res$features,
    final_features = ev$features,
    nonzero_features = ev$nonzero_features %||% ev$features,
    cv_cindex_mean = mean(cv_res$c_index_outer, na.rm = TRUE),
    cv_cindex_sd   = sd(cv_res$c_index_outer, na.rm = TRUE),
    cv_n_folds     = length(cv_res$c_index_outer),
    c_train = ev$c_train, c_test = ev$c_test,
    hr_test = ev$hr_test, hr_low = ev$hr_low,
    hr_high = ev$hr_high, p_test = ev$p_test,
    note = ev$note)

  all_selected_dt <- rbind(all_selected_dt, data.table(
    alpha = a,
    feature = ev$features,
    nonzero = ev$features %in% (ev$nonzero_features %||% ev$features),
    rank = seq_along(ev$features)), fill = TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# STEP 5 — summary table
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("STEP 5 — summary table")
log_msg("============================================================")

if (length(results) == 0) {
  stop("no alpha succeeded.")
}

summary_dt <- rbindlist(lapply(results, function(r) {
  data.table(
    alpha = r$alpha,
    n_consensus = r$n_features_consensus,
    n_final     = r$n_features_final,
    cv_cindex   = sprintf("%.3f ± %.3f", r$cv_cindex_mean, r$cv_cindex_sd),
    c_train     = round(r$c_train, 3),
    c_test      = round(r$c_test, 3),
    HR_test     = sprintf("%.3f (%.3f–%.3f)",
                          r$hr_test, r$hr_low, r$hr_high),
    P_test      = signif(r$p_test, 3),
    note        = r$note)
}))

print(summary_dt, row.names = FALSE)
write.csv(summary_dt,
          file.path(supp_dir, "TableS25_alpha_sensitivity.csv"),
          row.names = FALSE)
write.csv(all_selected_dt,
          file.path(supp_dir, "DataS3_alpha_selected_features.csv"),
          row.names = FALSE)

# =============================================================================
# STEP 6 — overlap with alpha = 0.5
# =============================================================================
log_msg("")
log_msg("STEP 6 — overlap with alpha = 0.5")

ref_key <- as.character(REF_ALPHA)
if (!is.null(results[[ref_key]])) {
  ref_features <- results[[ref_key]]$final_features
  overlap_dt <- rbindlist(lapply(setdiff(ALPHAS, REF_ALPHA), function(a) {
    a_key <- as.character(a)
    if (is.null(results[[a_key]])) return(NULL)
    f <- results[[a_key]]$final_features
    inter <- intersect(f, ref_features)
    union <- union(f, ref_features)
    data.table(
      alpha = a,
      n_features = length(f),
      n_overlap  = length(inter),
      jaccard    = round(length(inter) / length(union), 3),
      overlap_names = paste(inter, collapse = "; "),
      missing_names = paste(setdiff(ref_features, f), collapse = "; "))
  }), fill = TRUE)

  if (nrow(overlap_dt) > 0) {
    print(overlap_dt, row.names = FALSE)
    write.csv(overlap_dt,
              file.path(supp_dir, "TableS25b_alpha_overlap.csv"),
              row.names = FALSE)
  }
}

# =============================================================================
# STEP 7 — manifest
# =============================================================================
manifest <- list(
  step = "alpha_sensitivity",
  alphas_tested = ALPHAS, reference_alpha = REF_ALPHA,
  n_outer_folds = N_OUTER, n_inner_folds = N_INNER,
  n_repeats = N_REPEATS, consensus_threshold = CONSENSUS, epv_min = EPV_MIN,
  n_train = length(train_common),
  n_events_train = sum(y_train_full[, "status"] == 1),
  n_test = length(test_common),
  n_events_test = sum(y_test_full[, "status"] == 1),
  n_features_prefiltered = length(keep_features),
  timestamp = as.character(Sys.time()),
  r_version = R.version.string)
write_json(manifest,
           file.path(supp_dir, "alpha_sensitivity_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE)

log_msg("")
log_msg("outputs:")
log_msg("   TableS25_alpha_sensitivity.csv")
log_msg("   TableS25b_alpha_overlap.csv")
log_msg("   DataS3_alpha_selected_features.csv")
log_msg("   alpha_sensitivity_manifest.json")
log_msg("DONE")
