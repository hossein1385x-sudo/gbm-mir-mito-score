#!/usr/bin/env Rscript
# =============================================================================
# Supplementary Batch 2 — Table S14 + S18, Figure S6 + S7
# -----------------------------------------------------------------------------
# Table S14  — Locked model specification (reproducibility)
# Table S18  — Baseline characteristics (train / test / CGGA)
# Figure S6  — Schoenfeld residual plot
# Figure S7  — Linear vs spline comparison
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(splines)
  library(jsonlite)
})

input_dir <- paste0(WORK_DIR, "/")
sig_dir   <- file.path(input_dir, "signature_v2")
out_dir   <- file.path(input_dir, "supplementary_final")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

norm_id <- function(x) substr(toupper(gsub("\\.", "-", as.character(x))), 1, 12)

# =============================================================================
# TABLE S14 — Locked model specification
# =============================================================================
log_msg("")
log_msg("Table S14 — Locked model specification")

sig_path  <- file.path(sig_dir, "signature_v2_locked_latest.rds")
feat_path <- file.path(sig_dir, "signature_features.csv")
sum_path  <- file.path(sig_dir, "signature_v2_step12_summary.txt")

if (!file.exists(sig_path))  stop("signature_v2_locked_latest.rds is missing")

sig  <- readRDS(sig_path)
feat <- fread(feat_path)

# --- Part A: coefficients ---
coef_table <- data.frame(
  Feature     = feat$Feature,
  Coefficient = round(feat$Coefficient, 5),
  Type        = ifelse(grepl("^hsa-", feat$Feature), "miRNA", "clinical"),
  stringsAsFactors = FALSE)

# --- Part B: mean/SD on log2 ---
mirna_feats <- feat$Feature[grepl("^hsa-", feat$Feature)]
mean_sd_table <- data.frame(
  Feature = mirna_feats,
  train_mean_log2 = round(sig$scaling$mirna$center[mirna_feats], 4),
  train_sd_log2   = round(sig$scaling$mirna$scale[mirna_feats], 4),
  stringsAsFactors = FALSE)

# --- Part C: baseline survival (if stored in sig) ---
baseline_info <- NULL
if (!is.null(sig$baseline_survival)) {
  baseline_info <- data.frame(
    Horizon = c("6 months", "1 year", "2 years"),
    Days    = c(180, 365, 730),
    S0      = round(as.numeric(sig$baseline_survival), 4))
}

# --- Part D: SHA-256 ---
sha_info <- data.frame(
  Item    = c("SHA-256 fingerprint", "Locked on", "R version"),
  Value   = c(sig$meta$locked_on_sha %||% "stored_in_filename",
              as.character(sig$meta$locked_on),
              as.character(sig$meta$r_version)),
  stringsAsFactors = FALSE)

# --- save ---
write.csv(coef_table, file.path(out_dir, "TableS14a_coefficients.csv"),
          row.names = FALSE)
write.csv(mean_sd_table, file.path(out_dir, "TableS14b_mean_sd_log2.csv"),
          row.names = FALSE)

if (!is.null(baseline_info)) {
  write.csv(baseline_info, file.path(out_dir, "TableS14c_baseline_survival.csv"),
            row.names = FALSE)
}
write.csv(sha_info, file.path(out_dir, "TableS14d_manifest.csv"),
          row.names = FALSE)

log_msg("   TableS14a-d saved")
print(coef_table, row.names = FALSE)


# =============================================================================
# TABLE S18 — Baseline characteristics
# =============================================================================
log_msg("")
log_msg("Table S18 — Baseline characteristics")

clin_raw <- readRDS(file.path(input_dir, "raw_clinical_cdr_integrated.rds"))
if (!is.data.table(clin_raw)) clin_raw <- as.data.table(clin_raw)
clin_raw[, key := norm_id(patient_barcode)]

train_df <- fread(file.path(input_dir, "train_set.csv"))
test_df  <- fread(file.path(input_dir, "test_set_LOCKED.csv"))
train_keys <- norm_id(train_df$patient_barcode)
test_keys  <- norm_id(test_df$patient_barcode)

# --- CGGA ---
cgga_clin_file <- file.path(input_dir,
                            "cgga_raw/clin cgga/CGGA.microRNA_array_198_clinical.20200506.txt")
if (file.exists(cgga_clin_file)) {
  cgga_clin <- fread(cgga_clin_file)
  cgga_clin[, key := as.character(CGGA_ID)]
} else {
  cgga_clin <- data.table()
}

# --- summarising function ---
summarize_cohort <- function(clin_dt, keys, label) {
  sub <- clin_dt[key %in% keys]
  if (nrow(sub) == 0) return(NULL)

  # age
  age <- suppressWarnings(as.numeric(sub$age_at_diagnosis))
  if (median(age, na.rm = TRUE) > 200) age <- age / 365.25

  # survival time
  status <- as.integer(tolower(sub$vital_status) %in%
                         c("dead", "deceased", "1"))
  time <- ifelse(
    status == 1 & !is.na(sub$days_to_death) & sub$days_to_death > 0,
    as.numeric(sub$days_to_death),
    as.numeric(sub$days_to_last_follow_up))

  data.frame(
    Cohort = label,
    N = nrow(sub),
    Age_median = round(median(age, na.rm = TRUE), 1),
    Age_IQR    = sprintf("%.1f\u2013%.1f",
                         quantile(age, 0.25, na.rm = TRUE),
                         quantile(age, 0.75, na.rm = TRUE)),
    Male_pct   = round(100 * mean(tolower(sub$gender) == "male", na.rm = TRUE), 1),
    Dead_pct   = round(100 * mean(status == 1, na.rm = TRUE), 1),
    Followup_median = round(median(time, na.rm = TRUE), 0),
    stringsAsFactors = FALSE)
}

baseline_list <- list()
baseline_list$train <- summarize_cohort(clin_raw, train_keys, "TCGA Training")
baseline_list$test  <- summarize_cohort(clin_raw, test_keys,  "TCGA Test")

# CGGA separately (different columns)
if (nrow(cgga_clin) > 0) {
  cgga_age    <- suppressWarnings(as.numeric(cgga_clin$Age))
  cgga_status <- suppressWarnings(as.integer(cgga_clin[["Censor (alive=0; dead=1)"]]))
  cgga_os     <- suppressWarnings(as.numeric(cgga_clin$OS))

  baseline_list$cgga <- data.frame(
    Cohort = "CGGA",
    N = nrow(cgga_clin),
    Age_median = round(median(cgga_age, na.rm = TRUE), 1),
    Age_IQR    = sprintf("%.1f\u2013%.1f",
                         quantile(cgga_age, 0.25, na.rm = TRUE),
                         quantile(cgga_age, 0.75, na.rm = TRUE)),
    Male_pct   = round(100 * mean(cgga_clin$Gender == "male", na.rm = TRUE), 1),
    Dead_pct   = round(100 * mean(cgga_status == 1, na.rm = TRUE), 1),
    Followup_median = round(median(cgga_os, na.rm = TRUE), 0),
    stringsAsFactors = FALSE)
}

baseline <- rbindlist(baseline_list[!vapply(baseline_list, is.null, logical(1))])

write.csv(baseline,
          file.path(out_dir, "TableS18_baseline_characteristics.csv"),
          row.names = FALSE)

log_msg("   TableS18 saved")
print(baseline, row.names = FALSE)


# =============================================================================
# FIGURE S6 + S7 — Test-set diagnostics
# =============================================================================
log_msg("")
log_msg("Figure S6 + S7 — Test diagnostics")

pred_path <- file.path(sig_dir, "signature_v2_test_predictions.csv")
if (!file.exists(pred_path)) stop("test predictions file is missing.")

pred <- fread(pred_path)
y_test <- Surv(pred$survival_time, pred$event)
fit <- coxph(y_test ~ scale(pred$risk_score_full))

# --- Figure S6: Schoenfeld residuals ---
zph <- cox.zph(fit)
log_msg(sprintf("   cox.zph: chisq=%.3f, df=%d, p=%.4f",
                zph$table["GLOBAL", "chisq"],
                zph$table["GLOBAL", "df"],
                zph$table["GLOBAL", "p"]))

png_s6 <- file.path(out_dir, "FigureS6_Schoenfeld.png")
pdf_s6 <- file.path(out_dir, "FigureS6_Schoenfeld.pdf")

png(png_s6, width = 6, height = 5, units = "in", res = 600, bg = "white")
par(mar = c(4, 4, 3, 1))
plot(zph, var = 1, main = "Schoenfeld residuals (locked signature)",
     xlab = "Time (days)", ylab = "Beta(t) for risk score")
abline(h = coef(fit), lty = 2, col = "grey40")
dev.off()

pdf(pdf_s6, width = 6, height = 5, bg = "white")
par(mar = c(4, 4, 3, 1))
plot(zph, var = 1, main = "Schoenfeld residuals (locked signature)",
     xlab = "Time (days)", ylab = "Beta(t) for risk score")
abline(h = coef(fit), lty = 2, col = "grey40")
dev.off()

writeLines(c(
  "Figure S6. Schoenfeld residuals for the locked risk score on the test set, ",
  "with a loess smoother and 95% confidence band. The dashed horizontal line ",
  "indicates the estimated coefficient. The global proportional-hazards test ",
  sprintf("was not significant (\u03c7\u00b2 = %.2f, df = %d, P = %.2f), ",
          zph$table["GLOBAL", "chisq"], zph$table["GLOBAL", "df"],
          zph$table["GLOBAL", "p"]),
  "indicating no evidence of violation of the proportional-hazards assumption."
), file.path(out_dir, "FigureS6_Schoenfeld_caption.txt"))

log_msg("   FigureS6_Schoenfeld.png / .pdf")


# =============================================================================
# FIGURE S7 — Linear vs spline (corrected version)
# =============================================================================
log_msg("")
log_msg("Figure S7 — Linear vs spline")

# --- compute risk_z before Cox ---
pred$risk_z <- as.numeric(scale(pred$risk_score_full))
y_test <- Surv(pred$survival_time, pred$event)

# --- Cox on a single simple column ---
cox_lin <- coxph(y_test ~ risk_z, data = pred)
cox_spl <- coxph(y_test ~ ns(risk_z, df = 3), data = pred)

lrt <- anova(cox_lin, cox_spl, test = "LRT")
lrt_p     <- lrt[2, "Pr(>|Chi|)"]
lrt_chisq <- lrt$Chisq[2]
lrt_df    <- lrt$Df[2]

log_msg(sprintf("   LRT: chisq=%.3f, df=%d, p=%.4f", lrt_chisq, lrt_df, lrt_p))

# --- grid ---
risk_grid_z <- seq(min(pred$risk_z), max(pred$risk_z), length.out = 200)
new_data    <- data.frame(risk_z = risk_grid_z)

# --- predict with type="lp" and se ---
pred_spl <- predict(cox_spl, newdata = new_data, type = "lp", se.fit = TRUE)
log_msg(sprintf("   length of pred_spl$fit = %d (should be 200)", length(pred_spl$fit)))

if (length(pred_spl$fit) != length(risk_grid_z)) {
  stop("predict ran on the original data, not on new_data.")
}

# --- HR relative to the median ---
ref_idx   <- which.min(abs(risk_grid_z))
log_hr    <- pred_spl$fit - pred_spl$fit[ref_idx]
log_hr_se <- pred_spl$se.fit

df_spl <- data.frame(
  risk_z = risk_grid_z,
  hr     = exp(log_hr),
  ci_lo  = exp(log_hr - 1.96 * log_hr_se),
  ci_hi  = exp(log_hr + 1.96 * log_hr_se))

p_s7 <- ggplot(df_spl, aes(x = risk_z, y = hr)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              fill = "#D7301F", alpha = 0.20) +
  geom_line(color = "#D7301F", linewidth = 1.0) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  labs(
    title    = "Spline-based hazard-ratio curve (locked signature, test set)",
    subtitle = sprintf("LRT linear vs 3-df spline: \u03c7\u00b2 = %.2f, df = %d, P = %.3f",
                       lrt_chisq, lrt_df, lrt_p),
    x = "Risk score (standardized)",
    y = "Relative HR (vs median)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 10, color = "grey30"))

ggsave(file.path(out_dir, "FigureS7_Spline_HR.png"),
       p_s7, width = 7, height = 5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "FigureS7_Spline_HR.pdf"),
       p_s7, width = 7, height = 5)

writeLines(sprintf(paste0(
  "Figure S7. Estimated relative hazard ratio as a function of the locked ",
  "risk score in the test set, modeled using a natural cubic spline with 3 ",
  "degrees of freedom. The solid line is the point estimate; the shaded band ",
  "is the 95%% confidence interval. The dashed horizontal line marks HR = 1. ",
  "The likelihood-ratio test comparing the linear and spline models was not ",
  "significant (\u03c7\u00b2 = %.2f, df = %d, P = %.3f), consistent with a ",
  "linear effect on the log-hazard scale."),
  lrt_chisq, lrt_df, lrt_p),
  file.path(out_dir, "FigureS7_Spline_HR_caption.txt"))

log_msg("   FigureS7_Spline_HR.png / .pdf")

# =============================================================================
# summary
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("all outputs in: ", out_dir)
log_msg("   TableS14a_coefficients.csv")
log_msg("   TableS14b_mean_sd_log2.csv")
log_msg("   TableS14d_manifest.csv")
log_msg("   TableS18_baseline_characteristics.csv")
log_msg("   FigureS6_Schoenfeld.png / .pdf")
log_msg("   FigureS7_Spline_HR.png / .pdf")
log_msg("============================================================")
