#!/usr/bin/env Rscript
# =============================================================================
# Bivariate Cox adjustment (risk score + gender) in the TCGA test set
# =============================================================================

source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({ library(survival); library(data.table) })

input_dir <- paste0(WORK_DIR, "/")
supp_dir  <- file.path(input_dir, "supplementary_final")
norm_id   <- function(x) substr(toupper(gsub("\\.", "-", as.character(x))), 1, 12)

# --- DataS1 + risk scores ---
data_s1 <- fread(file.path(supp_dir, "DataS1_barcodes.csv"))
data_s1[, key := norm_id(patient_barcode_short)]
data_s1 <- data_s1[split == "test"]

pred <- fread(file.path(input_dir, "signature_v2/signature_v2_test_predictions.csv"))
pred[, key := norm_id(patient_id)]
pred_sub <- pred[, .(key, risk_score_full, survival_time, event)]

df <- merge(data_s1, pred_sub, by = "key", all.x = TRUE)
df <- df[!is.na(survival_time) & !is.na(event) &
           !is.na(risk_score_full) & survival_time > 0]

# --- Bivariate: gender ---
df_gender <- df[!is.na(gender)]
cat(sprintf("n for the gender adjustment: %d\n", nrow(df_gender)))

y <- Surv(df_gender$survival_time, df_gender$event)

fit <- coxph(y ~ scale(risk_score_full) + gender,
             data = df_gender)
s <- summary(fit)

cat("\n=== TCGA test: adjusted for gender ===\n")
print(s$coefficients)

# HR risk_score adjusted
b  <- s$coefficients["scale(risk_score_full)", "coef"]
se <- s$coefficients["scale(risk_score_full)", "se(coef)"]
pv <- s$coefficients["scale(risk_score_full)", "Pr(>|z|)"]

cat(sprintf("\nHR risk_score (adj gender): %.3f (95%% CI %.3f–%.3f, P = %.4g)\n",
            exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se), pv))

# save
result <- data.table(
  variable = "gender", cohort = "TCGA_test",
  n_used = nrow(df_gender),
  n_events = sum(df_gender$event),
  HR_risk = round(exp(b), 3),
  HR_low  = round(exp(b - 1.96 * se), 3),
  HR_high = round(exp(b + 1.96 * se), 3),
  p_value = signif(pv, 3),
  note = "OK")

write.csv(result,
          file.path(supp_dir, "bivariate_adjustment_test.csv"),
          row.names = FALSE)
cat("\nbivariate_adjustment_test.csv saved\n")
