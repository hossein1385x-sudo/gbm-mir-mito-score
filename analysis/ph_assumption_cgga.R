#!/usr/bin/env Rscript
# ============================================================
#  Proportional Hazards assumption test in CGGA
#  Using cgga_model_input/cgga_flat.csv
# ============================================================
rm(list = ls())

source(here::here("config.R"))  # defines WORK_DIR (see config.R)

ROOT <- WORK_DIR
OUT  <- file.path(ROOT, "PH_CGGA")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(survival); library(data.table)
})

## ---------- 1) Load ----------
d <- as.data.frame(fread(file.path(ROOT, "cgga_model_input/cgga_flat.csv")))
cat("Loaded n =", nrow(d), "\n")
cat("Columns:", paste(names(d), collapse = " | "), "\n\n")

## ---------- 2) Also load risk score ----------
# two paths: either from cgga_flat (if available) or from CGGA_sample_scores.csv
if ("risk_score" %in% names(d)) {
  cat("risk_score found in cgga_flat\n")
} else {
  rs <- read.csv(file.path(ROOT, "supplementary_final/CGGA_sample_scores.csv"),
                 stringsAsFactors = FALSE)
  d <- merge(d, rs[, c("CGGA_ID", "risk_score")],
             by.x = "id", by.y = "CGGA_ID", all.x = TRUE)
  cat("risk_score merged from CGGA_sample_scores.csv\n")
}

## ---------- 3) Prepare ----------
d$time  <- as.numeric(d$time)
d$event <- as.integer(d$status)
d$age   <- as.numeric(d$age)

d <- d[!is.na(d$time) & !is.na(d$event) & !is.na(d$risk_score) & d$time > 0, ]
cat("Analysis cohort: n =", nrow(d), " events =", sum(d$event), "\n\n")

## ---------- 4) Fit Cox models ----------
## Model A: risk score only
fit_a <- coxph(Surv(time, event) ~ risk_score, data = d)
cat("=== Model A: risk score only ===\n")
print(summary(fit_a)$coefficients)
cat("C-index:", round(summary(fit_a)$concordance[1], 4), "\n\n")

## Model B: risk score + age
d_b <- d[!is.na(d$age), ]
fit_b <- coxph(Surv(time, event) ~ risk_score + age, data = d_b)
cat("=== Model B: risk score + age ===\n")
print(summary(fit_b)$coefficients)
cat("C-index:", round(summary(fit_b)$concordance[1], 4), "\n\n")

## Model C: fully adjusted (risk + age + IDH + grade)
d_c <- d[!is.na(d$age) & !is.na(d$idh) & !is.na(d$grade), ]
fit_c <- coxph(Surv(time, event) ~ risk_score + age + idh + grade, data = d_c)
cat("=== Model C: risk + age + IDH + grade ===\n")
print(summary(fit_c)$coefficients)
cat("C-index:", round(summary(fit_c)$concordance[1], 4), "\n\n")

## ---------- 5) Schoenfeld residuals (PH test) ----------
cat("=========================================================\n")
cat("Proportional Hazards assumption (cox.zph)\n")
cat("=========================================================\n\n")

zph_a <- cox.zph(fit_a)
cat("--- Model A: risk_score only ---\n")
print(zph_a)
cat("\n")

zph_b <- cox.zph(fit_b)
cat("--- Model B: risk_score + age ---\n")
print(zph_b)
cat("\n")

zph_c <- cox.zph(fit_c)
cat("--- Model C: risk_score + age + IDH + grade ---\n")
print(zph_c)
cat("\n")

## ---------- 6) Save output ----------
sink(file.path(OUT, "PH_CGGA_summary.txt"))
cat("=== CGGA cohort: n =", nrow(d), "events =", sum(d$event), "===\n\n")
cat("--- Model A (risk only) ---\n"); print(zph_a)
cat("\n--- Model B (risk + age) ---\n"); print(zph_b)
cat("\n--- Model C (risk + age + IDH + grade) ---\n"); print(zph_c)
sink()

## Table for Supplementary
tbl_rows <- function(zph, model_name) {
  t <- as.data.frame(zph$table)
  t$covariate <- rownames(t)
  t$model     <- model_name
  rownames(t) <- NULL
  t[, c("model", "covariate", "chisq", "df", "p")]
}
ph_tab <- rbind(
  tbl_rows(zph_a, "risk only"),
  tbl_rows(zph_b, "risk + age"),
  tbl_rows(zph_c, "risk + age + IDH + grade")
)
write.csv(ph_tab, file.path(OUT, "PH_CGGA_table.csv"), row.names = FALSE)

cat("\nSaved:\n  ", file.path(OUT, "PH_CGGA_summary.txt"), "\n")
cat("  ", file.path(OUT, "PH_CGGA_table.csv"), "\n")

## ---------- 7) Interpretation ----------
cat("\n=========================================================\n")
cat("Interpretation guide\n")
cat("=========================================================\n")
cat("Global P > 0.05  ->  PH assumption holds for the model\n")
cat("Per-covariate P > 0.05  ->  PH holds for that covariate\n")
