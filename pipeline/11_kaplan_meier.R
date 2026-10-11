#!/usr/bin/env Rscript
# =============================================================================
# STEP 11 (Signature v2) — Kaplan-Meier Stratification
# -----------------------------------------------------------------------------
# Survival analysis of the v2 signature on the test set:
#   - Median split
#   - Quartile split
#   - Spline (checks linearity of the risk score - hazard relationship)
#   - Kaplan-Meier plots (600 DPI + 300 DPI + PDF)
#
# Outputs are saved in signature_v2/.
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(survminer)
  library(splines)
})

output_dir <- paste0(WORK_DIR, "/")
sig_dir    <- file.path(output_dir, "signature_v2")

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

log_msg("============================================================")
log_msg("STEP 11 (Signature v2) — Kaplan-Meier Stratification")
log_msg("============================================================")

# ======================================================================
# 1) Load the v2 predictions
# ======================================================================
log_msg("loading the v2 test predictions ...")

test_path <- file.path(sig_dir, "signature_v2_test_predictions.csv")
if (!file.exists(test_path))
  stop("not found: ", test_path)

test_scores <- read.csv(test_path, stringsAsFactors = FALSE)

# required columns
req <- c("patient_id", "risk_score_full", "survival_time", "event")
miss <- setdiff(req, names(test_scores))
if (length(miss) > 0)
  stop("required columns are missing: ", paste(miss, collapse = ", "))

# naming consistent with the main script
eval_set <- data.frame(
  patient_id    = test_scores$patient_id,
  risk_score    = test_scores$risk_score_full,
  survival_time = test_scores$survival_time,
  status        = as.integer(test_scores$event),
  stringsAsFactors = FALSE
)

# filter invalid values
keep <- is.finite(eval_set$survival_time) &
        !is.na(eval_set$status) &
        is.finite(eval_set$risk_score)
eval_set <- eval_set[keep, , drop = FALSE]
eval_set$valid_survival <- TRUE

n_eval   <- nrow(eval_set)
n_events <- sum(eval_set$status == 1)

log_msg(sprintf("   %d patients | %d events", n_eval, n_events))

surv_obj <- Surv(eval_set$survival_time, eval_set$status)

# ======================================================================
# 2) Valid time range
# ======================================================================
times_all   <- sort(unique(eval_set$survival_time))
n_risk_all  <- sapply(times_all, function(t) sum(eval_set$survival_time >= t))
valid_times <- times_all[n_risk_all >= 5]
max_valid_time <- if (length(valid_times) > 0) max(valid_times) else max(times_all)
log_msg(sprintf("   Max time with >=5 at risk: %.0f days", max_valid_time))

plot_time_limit <- min(max_valid_time, 730)   # 2 years
log_msg(sprintf("   Plot time limit: %.0f days (2 years)", plot_time_limit))

# at-risk table
risk_times <- c(180, 365, 548, 730)
risk_times <- risk_times[risk_times <= plot_time_limit]
risk_table <- data.frame(
  time_days = risk_times,
  n_risk    = sapply(risk_times, function(t) sum(eval_set$survival_time >= t)),
  n_events  = sapply(risk_times, function(t)
                     sum(eval_set$survival_time <= t & eval_set$status == 1))
)
write.csv(risk_table,
          file.path(sig_dir, "signature_v2_km_risk_table.csv"),
          row.names = FALSE)
log_msg("   signature_v2_km_risk_table.csv saved.")

# ======================================================================
# 3) Median split
# ======================================================================
log_msg("")
log_msg("Median split ...")

median_cut <- median(eval_set$risk_score)
group_med <- factor(
  ifelse(eval_set$risk_score > median_cut, "High risk", "Low risk"),
  levels = c("Low risk", "High risk")
)

fit_med   <- survfit(surv_obj ~ group_med)
sd_med    <- survdiff(surv_obj ~ group_med)
p_med     <- 1 - pchisq(sd_med$chisq, df = 1)

log_msg(sprintf("   Log-rank p = %.4g", p_med))

# ======================================================================
# 4) Quartile split
# ======================================================================
log_msg("")
log_msg("Quartile split ...")

q_cuts <- quantile(eval_set$risk_score, probs = c(0.25, 0.5, 0.75))
group_quart <- cut(
  eval_set$risk_score,
  breaks = c(-Inf, q_cuts, Inf),
  labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)")
)

fit_quart <- survfit(surv_obj ~ group_quart)
sd_quart  <- survdiff(surv_obj ~ group_quart)
df_quart  <- length(sd_quart$n) - 1
p_quart   <- 1 - pchisq(sd_quart$chisq, df = df_quart)

log_msg(sprintf("   Quartile log-rank p = %.4g", p_quart))

# ======================================================================
# 5) Spline analysis
# ======================================================================
log_msg("")
log_msg("Spline analysis ...")

cox_linear <- coxph(surv_obj ~ risk_score, data = eval_set)
cox_spline <- coxph(surv_obj ~ ns(risk_score, df = 3), data = eval_set)
lrt <- anova(cox_linear, cox_spline)
p_nonlin      <- lrt$`Pr(>|Chi|)`[2]
chisq_nonlin  <- lrt$Chisq[2]
df_nonlin     <- lrt$Df[2]
if (is.na(p_nonlin)) p_nonlin <- NA

if (!is.na(p_nonlin))
  log_msg(sprintf("   Spline LRT p = %.4g (linear vs 3-df)", p_nonlin))

# ======================================================================
# 6) Log-rank table
# ======================================================================
logrank_results <- data.frame(
  stratification = c("median_split", "quartile_based", "continuous_spline"),
  test_type = c("log-rank (2 groups)",
                sprintf("log-rank (%d groups)", length(sd_quart$n)),
                "LRT (linear vs 3-df spline)"),
  statistic = c(sd_med$chisq, sd_quart$chisq,
                ifelse(is.na(p_nonlin), NA, chisq_nonlin)),
  df = c(1, df_quart, ifelse(is.na(p_nonlin), NA, df_nonlin)),
  p_value = c(p_med, p_quart, p_nonlin),
  n = c(n_eval, n_eval, n_eval),
  n_events = c(n_events, n_events, n_events),
  notes = c("Cutpoint = median", "Cutpoints = Q1/Q2/Q3",
            "Test of non-linearity (no cutpoint)")
)
logrank_results$p_value_formatted <- sapply(logrank_results$p_value, function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.0001) return("< 0.0001")
  sprintf("%.4f", p)
})
write.csv(logrank_results,
          file.path(sig_dir, "signature_v2_logrank_results.csv"),
          row.names = FALSE)
log_msg("   signature_v2_logrank_results.csv saved.")

# ======================================================================
# 7) Helper for the p-value label
# ======================================================================
p_label <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.0001) return("p < 0.0001")
  if (p < 0.001)  return(sprintf("p = %.4f", p))
  sprintf("p = %.3f", p)
}

# ======================================================================
# 8) Plots
# ======================================================================
log_msg("")
log_msg(" building the plots ...")

# ---- Median split plot ----
p1 <- ggsurvplot(
  fit_med, data = eval_set,
  pval = p_label(p_med),
  pval.method = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  xlim = c(0, plot_time_limit),
  break.time.by = 182.5,
  palette = c("#2C7FB8", "#D7301F"),
  title = "Kaplan-Meier (Signature v2): Median split",
  legend.title = "Risk group",
  legend.labs = c("Low risk", "High risk"),
  xlab = "Time (days)",
  ylab = "Survival probability",
  ggtheme = theme_bw(base_size = 12),
  tables.theme = theme_cleantable()
)
p1$table <- p1$table + xlim(0, plot_time_limit)
p1$plot  <- p1$plot  + xlim(0, plot_time_limit)

ggsave(file.path(sig_dir, "signature_v2_km_median_600dpi.png"),
       plot = p1$plot, width = 8, height = 7, dpi = 600)
ggsave(file.path(sig_dir, "signature_v2_km_median_300dpi.png"),
       plot = p1$plot, width = 8, height = 7, dpi = 300)
log_msg("   signature_v2_km_median_600dpi.png")

# ---- Quartile split plot ----
p2 <- ggsurvplot(
  fit_quart, data = eval_set,
  pval = p_label(p_quart),
  pval.method = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  xlim = c(0, plot_time_limit),
  break.time.by = 182.5,
  palette = c("#2C7FB8", "#7FCDBB", "#FEB24C", "#D7301F"),
  title = "Kaplan-Meier (Signature v2): Quartile-based",
  legend.title = "Risk quartile",
  legend.labs = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
  xlab = "Time (days)",
  ylab = "Survival probability",
  ggtheme = theme_bw(base_size = 12),
  tables.theme = theme_cleantable()
)
p2$table <- p2$table + xlim(0, plot_time_limit)
p2$plot  <- p2$plot  + xlim(0, plot_time_limit)

ggsave(file.path(sig_dir, "signature_v2_km_quartile_600dpi.png"),
       plot = p2$plot, width = 8, height = 7, dpi = 600)
ggsave(file.path(sig_dir, "signature_v2_km_quartile_300dpi.png"),
       plot = p2$plot, width = 8, height = 7, dpi = 300)
log_msg("   signature_v2_km_quartile_600dpi.png")

# ---- Spline plot ----
p_spline <- NULL
if (!is.na(p_nonlin)) {
  risk_grid <- seq(min(eval_set$risk_score), max(eval_set$risk_score),
                   length.out = 200)
  pred <- predict(cox_spline,
                  newdata = data.frame(risk_score = risk_grid),
                  type = "lp", se.fit = TRUE, reference = "sample")
  ref_idx  <- which.min(abs(risk_grid - median(eval_set$risk_score)))
  log_hr   <- pred$fit - pred$fit[ref_idx]
  log_hr_se <- pred$se.fit
  ci_lo <- log_hr - 1.96 * log_hr_se
  ci_hi <- log_hr + 1.96 * log_hr_se

  df_spline <- data.frame(
    risk = risk_grid,
    hr   = exp(log_hr),
    ci_lo = exp(ci_lo),
    ci_hi = exp(ci_hi)
  )

  p_spline <- ggplot(df_spline, aes(x = risk, y = hr)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                fill = "#D7301F", alpha = 0.2) +
    geom_line(color = "#D7301F", linewidth = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    labs(
      x = "Risk score",
      y = "Relative HR (vs median)",
      title = sprintf("Spline (Signature v2): LRT p = %s",
                      ifelse(p_nonlin < 0.0001, "< 0.0001",
                             sprintf("%.4f", p_nonlin)))
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(file.path(sig_dir, "signature_v2_risk_spline_600dpi.png"),
         plot = p_spline, width = 8, height = 6, dpi = 600)
  log_msg("   signature_v2_risk_spline_600dpi.png")
}

# ---- PDF with all the plots ----
pdf(file.path(sig_dir, "signature_v2_km_curves_highres.pdf"),
    width = 8, height = 7)
print(p1$plot)
print(p2$plot)
if (!is.null(p_spline)) print(p_spline)
dev.off()
log_msg("   signature_v2_km_curves_highres.pdf")

# ======================================================================
# 9) Summary
# ======================================================================
log_msg("")
log_msg("============================================================")
log_msg("Summary of STEP 11 (Signature v2)")
log_msg("============================================================")
log_msg(sprintf("n = %d | events = %d", n_eval, n_events))
log_msg(sprintf("Median split    p = %s",
                ifelse(p_med < 0.0001, "< 0.0001", sprintf("%.4f", p_med))))
log_msg(sprintf("Quartile split  p = %s",
                ifelse(p_quart < 0.0001, "< 0.0001", sprintf("%.4f", p_quart))))
if (!is.na(p_nonlin))
  log_msg(sprintf("Spline LRT      p = %s",
                  ifelse(p_nonlin < 0.0001, "< 0.0001",
                         sprintf("%.4f", p_nonlin))))
log_msg(sprintf("Plot time limit = %.0f days", plot_time_limit))
log_msg("")
log_msg("outputs in: ", sig_dir)
log_msg("STEP 11 (Signature v2) COMPLETED")
