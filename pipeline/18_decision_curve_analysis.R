#!/usr/bin/env Rscript
# =============================================================================
# STEP 17 (Signature v2) — Decision Curve Analysis
# -----------------------------------------------------------------------------
# DCA for the v2 signature on the test set at 6 months, 1 year and 2 years
#
# Inputs:
#   signature_v2/signature_v2_locked_latest.rds
#   signature_v2/signature_v2_test_predictions.csv
#   train_set.csv, mirna_isoform_rpm_matrix.csv, raw_clinical_cdr_integrated.rds
#
# Outputs in signature_v2/17_dca/
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(dcurves)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(data.table)
})

# ======================================================================
# 0) Settings
# ======================================================================
output_dir <- paste0(WORK_DIR, "/")
sig_dir    <- file.path(output_dir, "signature_v2")
out_dir    <- file.path(sig_dir, "17_dca")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

norm_id <- function(x) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  substr(x, 1, 12)
}

log_msg("============================================================")
log_msg("STEP 17 (Signature v2) — Decision Curve Analysis")
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

# standardise the columns for DCA
if ("survival_time" %in% names(test_scores) &&
    !("time" %in% names(test_scores))) {
  test_scores$time <- test_scores$survival_time
}
if (!("time" %in% names(test_scores)))
  stop("column time or survival_time not found.")
if (!("event" %in% names(test_scores)))
  stop("column event not found.")

# risk score
if (!("risk_score_full" %in% names(test_scores)))
  stop("column risk_score_full is missing from the predictions.")

n_patients <- nrow(test_scores)
n_events   <- sum(test_scores$event == 1, na.rm = TRUE)
log_msg(sprintf("   DCA cohort: n = %d | events = %d", n_patients, n_events))

# ======================================================================
# 2) Rebuild train to compute the baseline survival
# ======================================================================
log_msg("")
log_msg("rebuilding train ...")

isoform  <- read.csv(file.path(output_dir, "mirna_isoform_rpm_matrix.csv"),
                     row.names = 1, check.names = FALSE)
clin     <- readRDS(file.path(output_dir, "raw_clinical_cdr_integrated.rds"))
train_df <- fread(file.path(output_dir, "train_set.csv"))

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

# miRNA scaling
mir_feats  <- sig$mirna_features
X_tr_mir   <- X_tr_all[, mir_feats, drop = FALSE]
mir_center <- sig$scaling$mirna$center[mir_feats]
mir_scale  <- sig$scaling$mirna$scale[mir_feats]

X_tr_mir_scaled <- sweep(X_tr_mir, 2, mir_center, "-")
X_tr_mir_scaled <- sweep(X_tr_mir_scaled, 2, mir_scale, "/")

# Age
age_raw <- suppressWarnings(as.numeric(cl_tr$age_at_diagnosis))
if (all(is.na(age_raw))) age_raw <- suppressWarnings(as.numeric(cl_tr$age))
age_years <- if (median(age_raw, na.rm = TRUE) > 200) age_raw / 365.25 else age_raw
age_years[is.na(age_years)] <- sig$scaling$clinical$center["Age"]

age_z <- (age_years - sig$scaling$clinical$center["Age"]) /
         sig$scaling$clinical$scale["Age"]

X_train_scaled <- cbind(X_tr_mir_scaled, Age = as.numeric(age_z))
y_train_surv   <- Surv(cl_tr$survival_time, cl_tr$status)

# ======================================================================
# 3) lp_train + baseline hazard
# ======================================================================
log_msg("")
log_msg("computing the baseline ...")

lp_train  <- as.numeric(predict(sig$model, newx = X_train_scaled, type = "link"))
cox_refit <- coxph(y_train_surv ~ offset(lp_train))
h0        <- basehaz(cox_refit, centered = FALSE)

H0_at <- function(t) {
  idx <- findInterval(t, h0$time)
  if (idx == 0) return(0)
  h0$hazard[idx]
}

# ======================================================================
# 4) Compute the event probability at the time horizons
# ======================================================================
log_msg("")
log_msg("computing the event probability ...")

time_horizons_days  <- c(180, 365, 730)
time_labels         <- c("p6mo", "p1yr", "p2yr")
time_labels_display <- c("6-month", "1-year", "2-year")

for (i in seq_along(time_horizons_days)) {
  t_i  <- time_horizons_days[i]
  S0_t <- exp(-H0_at(t_i))
  test_scores[[time_labels[i]]] <-
    1 - S0_t ^ exp(test_scores$risk_score_full)
}

for (lab in time_labels) {
  r <- range(test_scores[[lab]], na.rm = TRUE)
  log_msg(sprintf("   %-5s event prob: [%.3f, %.3f]", lab, r[1], r[2]))
}

# ======================================================================
# 5) DCA for each horizon
# ======================================================================
log_msg("")
log_msg("DCA per horizon ...")

thresholds       <- seq(0.01, 0.50, by = 0.01)
dca_results_list <- list()

for (i in seq_along(time_labels)) {
  lab_col     <- time_labels[i]
  lab_display <- time_labels_display[i]

  formula_str <- as.formula(paste0("Surv(time, event) ~ `", lab_col, "`"))

  dca_fit <- tryCatch(
    dcurves::dca(
      formula_str,
      data       = test_scores,
      time       = time_horizons_days[i],
      thresholds = thresholds,
      label      = setNames(list("Signature v2"), lab_col)
    ),
    error = function(e) { log_msg("   ", conditionMessage(e)); NULL }
  )

  if (is.null(dca_fit)) next

  dca_df <- as_tibble(dca_fit$dca) %>%
    mutate(time_horizon = lab_display)

  dca_results_list[[lab_col]] <- dca_df
  log_msg(sprintf("   %s", lab_display))
}

if (length(dca_results_list) == 0) stop("no DCA was computed.")
dca_all <- bind_rows(dca_results_list)

# ======================================================================
# 6) Crossing detection
# ======================================================================
log_msg("")
log_msg("detecting crossings ...")

detect_crossings <- function(df) {
  df_wide <- df %>%
    select(threshold, variable, net_benefit) %>%
    pivot_wider(names_from = variable, values_from = net_benefit)

  model_col <- setdiff(colnames(df_wide), c("threshold", "all", "none"))
  if (length(model_col) == 0) {
    return(tibble(time_horizon = unique(df$time_horizon),
                  crossing_found = FALSE))
  }

  model_col <- model_col[1]
  df_wide <- df_wide %>%
    mutate(diff_vs_all = .data[[model_col]] - all) %>%
    arrange(threshold)

  sign_changes <- which(diff(sign(df_wide$diff_vs_all)) != 0)
  if (length(sign_changes) == 0) {
    return(tibble(time_horizon = unique(df$time_horizon),
                  crossing_found = FALSE))
  }

  tibble(
    time_horizon       = unique(df$time_horizon),
    crossing_found     = TRUE,
    crossing_threshold = list(df_wide$threshold[sign_changes])
  )
}

crossing_diagnostics <- dca_all %>%
  group_by(time_horizon) %>%
  group_split() %>%
  map_dfr(detect_crossings)

if (any(crossing_diagnostics$crossing_found)) {
  crossing_with_counts <- crossing_diagnostics %>%
    filter(crossing_found) %>%
    rowwise() %>%
    mutate(
      lab_col = time_labels[match(time_horizon, time_labels_display)],
      crossing_threshold_flat = paste(crossing_threshold[[1]], collapse = ", "),
      n_classified_high_risk = sum(test_scores[[lab_col]] >
                                     crossing_threshold[[1]][1], na.rm = TRUE),
      n_events_among_high_risk = sum(test_scores[[lab_col]] >
                                       crossing_threshold[[1]][1] &
                                     test_scores$event == 1, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(time_horizon, crossing_threshold_flat,
           n_classified_high_risk, n_events_among_high_risk)

  log_msg("Crossing detected:")
  print(as.data.frame(crossing_with_counts))
  write.csv(crossing_with_counts,
            file.path(out_dir, "dca_crossing_diagnostics.csv"),
            row.names = FALSE)
} else {
  crossing_with_counts <- tibble()
  log_msg("no crossing observed at any horizon.")
}

# ======================================================================
# 7) DCA plot
# ======================================================================
LOW_DATA_N_CUTOFF <- 10

n_by_threshold <- map_dfr(seq_along(time_labels), function(i) {
  lab_col     <- time_labels[i]
  lab_display <- time_labels_display[i]
  tibble(
    threshold    = thresholds,
    time_horizon = lab_display,
    n_high_risk  = sapply(thresholds, function(th)
                          sum(test_scores[[lab_col]] > th, na.rm = TRUE))
  )
})

low_data_zones <- n_by_threshold %>% filter(n_high_risk < LOW_DATA_N_CUTOFF)

low_data_rects <- if (nrow(low_data_zones) > 0) {
  low_data_zones %>%
    group_by(time_horizon) %>%
    summarise(xmin = min(threshold), xmax = max(threshold), .groups = "drop")
} else NULL

p_dca <- ggplot(dca_all, aes(x = threshold, y = net_benefit, color = label)) +
  {if (!is.null(low_data_rects))
    geom_rect(
      data = low_data_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "grey85", alpha = 0.5
    )} +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  facet_wrap(~ time_horizon) +
  coord_cartesian(ylim = c(-0.05, NA)) +
  labs(
    title    = "Decision Curve Analysis — Signature v2 (Test Set)",
    subtitle = sprintf("n = %d | %d events | grey zone = threshold with n(high-risk) < %d",
                       n_patients, n_events, LOW_DATA_N_CUTOFF),
    x        = "Threshold probability",
    y        = "Net benefit",
    color    = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "dca_curves.pdf"),
       p_dca, width = 11, height = 5)
ggsave(file.path(out_dir, "dca_curves_600dpi.png"),
       p_dca, width = 11, height = 5, dpi = 600)

write.csv(dca_all,
          file.path(out_dir, "dca_net_benefit_table.csv"),
          row.names = FALSE)

# ======================================================================
# 8) Final summary
# ======================================================================
log_msg("")
log_msg("============================================================")
log_msg("Summary of STEP 17 (Signature v2)")
log_msg("============================================================")
log_msg(sprintf("n = %d | events = %d", n_patients, n_events))
log_msg(sprintf("horizons: %s", paste(time_labels_display, collapse = ", ")))
log_msg(sprintf("Crossing: %s",
                ifelse(any(crossing_diagnostics$crossing_found),
                       "yes", "no")))
log_msg(sprintf("low-data zones (n<10): %s",
                ifelse(nrow(low_data_zones) > 0, "yes", "no")))
log_msg("")
log_msg("outputs in: ", out_dir)
log_msg("STEP 17 (Signature v2) COMPLETED")
