#!/usr/bin/env Rscript
# =============================================================================
# STEP 09 (Signature v2) — Apply Locked Signature to Test
# -----------------------------------------------------------------------------
# Apply the locked v2 signature (6 miRNA + Age) to the test set, without refitting.
# Outputs are saved in signature_v2/ so they do not collide with the original version.
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(jsonlite)
  library(digest)
})

output_dir <- paste0(WORK_DIR, "/")
sig_dir    <- file.path(output_dir, "signature_v2")

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

to_named_numeric <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.data.frame(x)) x <- x[[1]]
  if (is.list(x) && !is.numeric(x)) x <- unlist(x, use.names = TRUE)
  v <- as.numeric(x)
  names(v) <- names(x)
  v
}

# ======================================================================
# 1) Load the v2 signature + verify the hash
# ======================================================================
log_msg("============================================================")
log_msg("STEP 09 (Signature v2) — Test Evaluation")
log_msg("============================================================")
log_msg("loading ...")

sig_path       <- file.path(sig_dir, "signature_v2_locked_latest.rds")
test_path      <- file.path(output_dir, "test_set_LOCKED.csv")
manifest8_path <- file.path(sig_dir, "step08_signature_v2_manifest.json")

if (!file.exists(sig_path))  stop("signature_v2_locked_latest.rds is missing")
if (!file.exists(test_path)) stop("test_set_LOCKED.csv is missing")

signature_object <- readRDS(sig_path)
test_df          <- read.csv(test_path, stringsAsFactors = FALSE)

log_msg(sprintf("signature v2: %d features (miRNA: %d | clinical: %d)",
                nrow(signature_object$coefficients),
                length(signature_object$mirna_features),
                length(signature_object$clinical_features)))
log_msg(sprintf("test: %d patients", nrow(test_df)))

# ---- verify the hash against the Step 08 manifest ----
if (file.exists(manifest8_path)) {
  manifest8 <- fromJSON(manifest8_path, simplifyVector = FALSE)
  sig_hash_now <- substr(digest(list(
    coefficients      = signature_object$coefficients,
    mirna_features    = signature_object$mirna_features,
    clinical_features = signature_object$clinical_features,
    scaling           = signature_object$scaling,
    alpha             = signature_object$alpha
  ), algo = "sha256"), 1, 12)

  if (!is.null(manifest8$signature_hash_short) &&
      manifest8$signature_hash_short != sig_hash_now) {
    stop(sprintf("signature hash does not match the manifest. manifest=%s current=%s",
                 manifest8$signature_hash_short, sig_hash_now))
  }
  log_msg(sprintf("hash verified: %s", sig_hash_now))
}

# ======================================================================
# 2) Load the matrices
# ======================================================================
isoform_data  <- read.csv(
  file.path(output_dir, "mirna_isoform_rpm_matrix.csv"),
  row.names = 1, check.names = FALSE
)
clinical_data <- readRDS(file.path(output_dir, "raw_clinical_cdr_integrated.rds"))

# ======================================================================
# 3) Synchronisation
# ======================================================================
test_ids       <- normalize_tcga_id(test_df$patient_barcode)
iso_ids_clean  <- normalize_tcga_id(rownames(isoform_data))
clin_ids_clean <- normalize_tcga_id(as.character(clinical_data$patient_barcode))

test_present <- Reduce(intersect, list(test_ids, iso_ids_clean, clin_ids_clean))
log_msg(sprintf("   - shared test samples: %d of %d",
                length(test_present), length(test_ids)))

if (length(test_present) == 0) stop("no test sample has molecular data available")

iso_test  <- isoform_data[match(test_present, iso_ids_clean), , drop = FALSE]
rownames(iso_test) <- test_present

clin_test <- as.data.frame(clinical_data[match(test_present, clin_ids_clean), , drop = FALSE])
rownames(clin_test) <- test_present

# ======================================================================
# 4) Survival
# ======================================================================
if (!("status_final" %in% colnames(clin_test)))       stop("column status_final is missing")
if (!("survival_time_final" %in% colnames(clin_test))) stop("column survival_time_final is missing")

clin_test$status        <- clin_test$status_final
clin_test$survival_time <- clin_test$survival_time_final

keep <- !is.na(clin_test$status) &
        !is.na(clin_test$survival_time) &
        clin_test$survival_time >= 0
log_msg(sprintf("   - removed for missing survival: %d", sum(!keep)))

clin_test <- clin_test[keep, , drop = FALSE]
iso_test  <- iso_test[rownames(clin_test), , drop = FALSE]
if (nrow(clin_test) == 0) stop("no sample left")

# ======================================================================
# 5) log2 + clinical variable
# ======================================================================
X_test <- log2(as.matrix(iso_test) + 1)

scaling <- signature_object$scaling

# ---- Age: transformed exactly as in train ----
age_days <- NULL
if ("age_at_diagnosis" %in% colnames(clin_test)) {
  age_days <- suppressWarnings(as.numeric(clin_test$age_at_diagnosis))
} else if ("age" %in% colnames(clin_test)) {
  age_days <- suppressWarnings(as.numeric(clin_test$age))
}

if (is.null(age_days) || all(is.na(age_days))) {
  age_years <- rep(NA_real_, nrow(clin_test))
} else if (median(age_days, na.rm = TRUE) > 200) {
  age_years <- age_days / 365.25
  log_msg(sprintf("   - Age converted from days to years: median %.1f → %.1f",
                  median(age_days, na.rm = TRUE),
                  median(age_years, na.rm = TRUE)))
} else {
  age_years <- age_days
}

# ---- fill NA with the train centre ----
clin_center_all <- to_named_numeric(scaling$clinical$center)
if (any(is.na(age_years))) {
  age_center_train <- clin_center_all["Age"]
  if (is.na(age_center_train) || is.null(age_center_train)) {
    age_center_train <- median(age_years, na.rm = TRUE)
  }
  n_na <- sum(is.na(age_years))
  age_years[is.na(age_years)] <- age_center_train
  log_msg(sprintf("   - %d NA values in Age filled with the train median (%.1f years).",
                  n_na, age_center_train))
}

clin_test_cov <- data.frame(Age = age_years, row.names = rownames(clin_test))

clin_feats <- scaling$clinical$features
if (!identical(clin_feats, "Age")) {
  stop(sprintf("this version supports only 'Age', but scaling contains: %s",
               paste(clin_feats, collapse = ", ")))
}

# ======================================================================
# 6) transform with the locked parameters (no refit)
# ======================================================================
mir_center_all  <- to_named_numeric(scaling$mirna$center)
mir_scale_all   <- to_named_numeric(scaling$mirna$scale)
clin_center_all <- to_named_numeric(scaling$clinical$center)
clin_scale_all  <- to_named_numeric(scaling$clinical$scale)

mir_feats <- scaling$mirna$features

# ---- miRNA ----
miss_mir <- setdiff(mir_feats, colnames(X_test))
if (length(miss_mir) > 0) {
  stop(sprintf("miRNAs missing in test: %s", paste(miss_mir, collapse = ", ")))
}

X_test_mir <- X_test[, mir_feats, drop = FALSE]
mir_center_sub <- mir_center_all[mir_feats]
mir_scale_sub  <- mir_scale_all[mir_feats]

if (any(is.na(mir_center_sub)) || any(is.na(mir_scale_sub))) {
  stop("center/scale is NA for some miRNAs")
}
if (any(mir_scale_sub == 0)) {
  stop("zero scale found")
}

X_test_mir_scaled <- sweep(X_test_mir, 2, mir_center_sub, "-")
X_test_mir_scaled <- sweep(X_test_mir_scaled, 2, mir_scale_sub, "/")

# ---- clinical ----
X_test_clin <- as.matrix(clin_test_cov[, clin_feats, drop = FALSE])
storage.mode(X_test_clin) <- "double"

clin_center_sub <- clin_center_all[clin_feats]
clin_scale_sub  <- clin_scale_all[clin_feats]
if (any(is.na(clin_center_sub)) || any(is.na(clin_scale_sub))) {
  stop("clinical center/scale is NA")
}
clin_scale_sub[clin_scale_sub == 0] <- 1

X_test_clin_scaled <- sweep(X_test_clin, 2, clin_center_sub, "-")
X_test_clin_scaled <- sweep(X_test_clin_scaled, 2, clin_scale_sub, "/")

# ---- final combination ----
all_features <- c(mir_feats, clin_feats)
X_test_final <- cbind(X_test_mir_scaled, X_test_clin_scaled)
X_test_final <- as.matrix(X_test_final[, all_features, drop = FALSE])
storage.mode(X_test_final) <- "double"

log_msg(sprintf("   - dimensions of X_test_final: %d × %d",
                nrow(X_test_final), ncol(X_test_final)))

# ======================================================================
# 7) predict (no refit)
# ======================================================================
lp_pred <- predict(signature_object$model,
                   newx = X_test_final,
                   type = "link")

if (is.list(lp_pred)) lp_pred <- lp_pred[[1]]
if (!is.null(dim(lp_pred))) lp_pred <- lp_pred[, 1]
lp_test <- as.numeric(lp_pred)

if (any(is.na(lp_test))) stop("lp_test contains NA")

log_msg(sprintf("   - lp_test: range [%.3f, %.3f] | SD=%.3f",
                min(lp_test), max(lp_test), sd(lp_test)))

y_test <- Surv(time = clin_test$survival_time, event = clin_test$status)

# ======================================================================
# 8) Evaluation (no refit)
# ======================================================================
c_test <- survival::survConcordance(y_test ~ lp_test)$concordance

cox_uni <- coxph(y_test ~ lp_test)
cox_s   <- summary(cox_uni)

beta    <- coef(cox_uni)[1]
se_beta <- cox_s$coefficients[1, "se(coef)"]
p_uni   <- cox_s$coefficients[1, "Pr(>|z|)"]

# ---- HR per 1 unit ----
hr_per_unit <- exp(beta)
hr_unit_ci  <- exp(c(beta - 1.96 * se_beta, beta + 1.96 * se_beta))

# ---- HR per 1-SD ----
sd_lp          <- sd(lp_test)
hr_per_sd      <- exp(beta * sd_lp)
hr_per_sd_low  <- exp((beta - 1.96 * se_beta) * sd_lp)
hr_per_sd_high <- exp((beta + 1.96 * se_beta) * sd_lp)

log_msg("")
log_msg(strrep("=", 70))
log_msg("Test Performance — Signature v2 (no refit)")
log_msg(strrep("=", 70))
log_msg(sprintf("   n test          : %d", nrow(clin_test)))
log_msg(sprintf("   events          : %d", sum(clin_test$status == 1)))
log_msg(sprintf("   C-index         : %.4f", c_test))
log_msg(sprintf("   HR per 1-SD     : %.3f (95%% CI: %.3f–%.3f)",
                hr_per_sd, hr_per_sd_low, hr_per_sd_high))
log_msg(sprintf("   HR per 1 unit   : %.3f (95%% CI: %.3f–%.3f)",
                hr_per_unit, hr_unit_ci[1], hr_unit_ci[2]))
log_msg(sprintf("   p-value Cox     : %.4g", p_uni))
log_msg(strrep("=", 70))

# ======================================================================
# 9) Save in signature_v2/
# ======================================================================
write.csv(data.frame(
  patient_id    = rownames(clin_test),
  risk_score    = lp_test,
  survival_time = clin_test$survival_time,
  event         = clin_test$status,
  stringsAsFactors = FALSE),
  file.path(sig_dir, "signature_v2_test_predictions.csv"),
  row.names = FALSE)

write.csv(data.frame(
  metric = c("n_test", "n_events", "C_index",
             "HR_per_unit", "HR_per_unit_CI_low", "HR_per_unit_CI_high",
             "HR_per_SD",   "HR_per_SD_CI_low",   "HR_per_SD_CI_high",
             "p_value", "sd_lp"),
  value  = c(nrow(clin_test), sum(clin_test$status == 1), c_test,
             hr_per_unit, hr_unit_ci[1], hr_unit_ci[2],
             hr_per_sd, hr_per_sd_low, hr_per_sd_high,
             p_uni, sd_lp),
  stringsAsFactors = FALSE),
  file.path(sig_dir, "signature_v2_test_performance.csv"),
  row.names = FALSE)

# readable summary
sink(file.path(sig_dir, "signature_v2_step09_summary.txt"))
cat("========== Step 09 (Signature v2) — Test Evaluation ==========\n")
cat("Date:", as.character(Sys.time()), "\n")
cat(sprintf("Locked on     : %s\n", signature_object$meta$locked_on))
cat(sprintf("n test        : %d\n", nrow(clin_test)))
cat(sprintf("n events      : %d\n", sum(clin_test$status == 1)))
cat(sprintf("C-index       : %.4f\n", c_test))
cat(sprintf("HR per 1-SD   : %.4f (95%% CI: %.4f–%.4f, p=%.4g)\n",
            hr_per_sd, hr_per_sd_low, hr_per_sd_high, p_uni))
cat(sprintf("HR per unit   : %.4f (95%% CI: %.4f–%.4f)\n",
            hr_per_unit, hr_unit_ci[1], hr_unit_ci[2]))
sink()

log_msg("output files saved in: ", sig_dir)
log_msg("   - signature_v2_test_predictions.csv")
log_msg("   - signature_v2_test_performance.csv")
log_msg("   - signature_v2_step09_summary.txt")
log_msg("STEP 09 (Signature v2) COMPLETED")
