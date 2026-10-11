#!/usr/bin/env Rscript
# =============================================================================
# STEP 14.4 (Signature v2) — Calibration Plot
# -----------------------------------------------------------------------------
# Calibration of the v2 signature on the test set
#   - Predicted vs observed survival at 6 months / 1 year / 2 years
#   - Grouping by quartiles of predicted survival (4 groups)
#
# Inputs:
#   signature_v2/signature_v2_locked_latest.rds
#   signature_v2/signature_v2_test_predictions.csv
#   train_set.csv, mirna_isoform_rpm_matrix.csv, raw_clinical_cdr_integrated.rds
#
# Outputs in signature_v2/14_4_calibration/
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

set.seed(2024)

output_dir <- paste0(WORK_DIR, "/")
sig_dir    <- file.path(output_dir, "signature_v2")
out_dir    <- file.path(sig_dir, "14_4_calibration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

norm_id <- function(x) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  substr(x, 1, 12)
}

log_msg("============================================================")
log_msg("STEP 14.4 (Signature v2) — Calibration")
log_msg("============================================================")

# ======================================================================
# 1) Load the signature + test predictions
# ======================================================================
log_msg("loading ...")

sig_path  <- file.path(sig_dir, "signature_v2_locked_latest.rds")
test_path <- file.path(sig_dir, "signature_v2_test_predictions.csv")

if (!file.exists(sig_path))  stop("not found: ", sig_path)
if (!file.exists(test_path)) stop("not found: ", test_path)

sig         <- readRDS(sig_path)
test_scores <- as.data.frame(fread(test_path))

log_msg(sprintf("   miRNAs in the signature: %d | clinical: %s",
                length(sig$mirna_features),
                paste(sig$clinical_features, collapse = ", ")))
log_msg(sprintf("   test: %d patients | %d events",
                nrow(test_scores), sum(test_scores$event == 1)))

# ======================================================================
# 2) Rebuild the train X_scaled (to compute the baseline)
# ======================================================================
log_msg("")
log_msg("rebuilding train ...")

isoform  <- read.csv(file.path(output_dir, "mirna_isoform_rpm_matrix.csv"),
                     row.names = 1, check.names = FALSE)
clin     <- readRDS(file.path(output_dir, "raw_clinical_cdr_integrated.rds"))
train_df <- fread(file.path(output_dir, "train_set.csv"))

# primary filter
mirna_raw <- readRDS(file.path(output_dir, "raw_mirna.rds"))
if (!is.data.table(mirna_raw)) mirna_raw <- as.data.table(mirna_raw)
primary <- unique(norm_id(mirna_raw[is_primary_tumor == TRUE, sample_barcode]))

train_keys <- norm_id(train_df$patient_barcode)
iso_keys   <- norm_id(rownames(isoform))
clin_keys  <- norm_id(as.character(clin$patient_barcode))

common <- Reduce(intersect, list(train_keys, iso_keys, clin_keys, primary))
log_msg(sprintf("   shared train samples: %d", length(common)))

idx_iso  <- match(common, iso_keys)
idx_clin <- match(common, clin_keys)

X_tr_all <- log2(as.matrix(isoform[idx_iso, , drop = FALSE]) + 1)
cl_tr    <- as.data.frame(clin[idx_clin, , drop = FALSE])
rownames(X_tr_all) <- common
rownames(cl_tr)    <- common

cl_tr$status        <- cl_tr$status_final
cl_tr$survival_time <- cl_tr$survival_time_final

keep <- !is.na(cl_tr$status) & !is.na(cl_tr$survival_time) &
        cl_tr$survival_time >= 0
cl_tr    <- cl_tr[keep, , drop = FALSE]
X_tr_all <- X_tr_all[rownames(cl_tr), , drop = FALSE]

# --- miRNA scaling with the locked parameters ---
mir_feats <- sig$mirna_features
X_tr_mir  <- X_tr_all[, mir_feats, drop = FALSE]

mir_center <- sig$scaling$mirna$center[mir_feats]
mir_scale  <- sig$scaling$mirna$scale[mir_feats]

X_tr_mir_scaled <- sweep(X_tr_mir, 2, mir_center, "-")
X_tr_mir_scaled <- sweep(X_tr_mir_scaled, 2, mir_scale, "/")

# --- Age ---
age_raw <- suppressWarnings(as.numeric(cl_tr$age_at_diagnosis))
if (all(is.na(age_raw))) age_raw <- suppressWarnings(as.numeric(cl_tr$age))
if (median(age_raw, na.rm = TRUE) > 200) {
  age_years <- age_raw / 365.25
} else {
  age_years <- age_raw
}
age_years[is.na(age_years)] <- sig$scaling$clinical$center["Age"]

age_z <- (age_years - sig$scaling$clinical$center["Age"]) /
         sig$scaling$clinical$scale["Age"]

X_train_scaled <- cbind(X_tr_mir_scaled, Age = as.numeric(age_z))
y_train_surv   <- Surv(cl_tr$survival_time, cl_tr$status)

# ======================================================================
# 3) Compute lp_train and the baseline hazard
# ======================================================================
log_msg("")
log_msg("computing the baseline survival ...")

lp_train <- as.numeric(predict(sig$model, newx = X_train_scaled, type = "link"))

# Fit Cox with an offset (the lp_train coefficient is fixed at 1)
cox_refit <- coxph(y_train_surv ~ offset(lp_train))

# Baseline cumulative hazard at lp = 0
h0 <- basehaz(cox_refit, centered = FALSE)

H0_at <- function(t) {
  idx <- findInterval(t, h0$time)
  if (idx == 0) return(0)
  h0$hazard[idx]
}

log_msg(sprintf("   lp_train range: [%.3f, %.3f] | mean = %.3f",
                min(lp_train), max(lp_train), mean(lp_train)))

# ======================================================================
# 4) Predicted survival on test
# ======================================================================
log_msg("")
log_msg("computing the predicted survival ...")

# use the risk_score_full column from step 9
if (!("risk_score_full" %in% names(test_scores)))
  stop("column risk_score_full is missing from the predictions.")

time_horizons_days <- c(180, 365, 730)
time_labels        <- c("6mo", "1yr", "2yr")

for (i in seq_along(time_horizons_days)) {
  t_i  <- time_horizons_days[i]
  S0_t <- exp(-H0_at(t_i))
  test_scores[[time_labels[i]]] <-
    S0_t ^ exp(test_scores$risk_score_full)
}

for (lab in time_labels) {
  r <- range(test_scores[[lab]], na.rm = TRUE)
  log_msg(sprintf("   %-5s predicted survival: [%.3f, %.3f]", lab, r[1], r[2]))
}

# ======================================================================
# 5) Quartile grouping + compute the observed values
# ======================================================================
log_msg("")
log_msg("computing the calibration groups ...")

n_groups <- 4
calibration_results <- list()

for (i in seq_along(time_labels)) {
  t_lab  <- time_labels[i]
  t_days <- time_horizons_days[i]

  df_t <- test_scores
  df_t$pred_surv <- df_t[[t_lab]]

  qs <- quantile(df_t$pred_surv, probs = c(0, 0.25, 0.5, 0.75, 1),
                 na.rm = TRUE)

  if (length(unique(qs)) < 5) {
    df_t$risk_group <- dplyr::ntile(-df_t$pred_surv, n_groups)
  } else {
    df_t$risk_group <- as.integer(
      cut(df_t$pred_surv, breaks = qs,
          include.lowest = TRUE, labels = FALSE)
    )
  }

  group_summary <- df_t %>%
    group_by(risk_group) %>%
    summarise(
      n_group        = n(),
      n_event_group  = sum(event == 1, na.rm = TRUE),
      mean_predicted = mean(pred_surv, na.rm = TRUE),
      .groups = "drop"
    )

  observed_list <- lapply(sort(unique(df_t$risk_group)), function(g) {
    sub <- df_t[df_t$risk_group == g, , drop = FALSE]
    km  <- survfit(Surv(survival_time, event) ~ 1, data = sub)
    s   <- summary(km, times = t_days, extend = TRUE)
    data.frame(
      risk_group     = g,
      observed_surv  = s$surv,
      observed_lower = s$lower,
      observed_upper = s$upper
    )
  })
  observed_df <- do.call(rbind, observed_list)

  merged <- merge(group_summary, observed_df, by = "risk_group")
  merged$time_horizon <- t_lab
  calibration_results[[t_lab]] <- merged
}

calibration_df <- do.call(rbind, calibration_results)
calibration_df$time_horizon <- factor(calibration_df$time_horizon,
                                       levels = time_labels)

write.csv(calibration_df,
          file.path(out_dir, "calibration_data.csv"),
          row.names = FALSE)

log_msg("")
print(calibration_df)

# ======================================================================
# 6) Draw the calibration plot
# ======================================================================
p_calibration <- ggplot(calibration_df,
                        aes(x = mean_predicted, y = observed_surv)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50") +
  geom_errorbar(aes(ymin = observed_lower, ymax = observed_upper),
                width = 0.02, color = "steelblue") +
  geom_point(aes(size = n_group), color = "steelblue") +
  facet_wrap(~ time_horizon, nrow = 1) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title    = "Calibration — Signature v2 (Test Set)",
    subtitle = sprintf("n = %d patients, %d events | dashed = perfect calibration",
                       nrow(test_scores), sum(test_scores$event == 1)),
    x        = "Predicted survival (group mean)",
    y        = "Observed survival (KM ± 95% CI)",
    size     = "n per group"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "calibration_plot.pdf"),
       p_calibration, width = 10, height = 4)
ggsave(file.path(out_dir, "calibration_plot_600dpi.png"),
       p_calibration, width = 10, height = 4, dpi = 600)

log_msg("")
log_msg("outputs in: ", out_dir)
log_msg("STEP 14.4 (Signature v2) COMPLETED")
