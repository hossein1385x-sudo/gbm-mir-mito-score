#!/usr/bin/env Rscript
# =============================================================================
# pipeline/11_km_train_cutoff.R
# Kaplan-Meier with a cutoff computed from TRAIN only (locked)
# Train: train_full_preprocessed.rds + Age from train_set.csv
# Test and CGGA: supplementary_final/DataS2_risk_scores.csv
# =============================================================================

suppressPackageStartupMessages({
  library(survival); library(survminer)
  library(dplyr);   library(readr); library(ggplot2)
})

# ---------------------------------------------------------------- 0) Paths
source(here::here("config.R"))   # defines WORK_DIR (see README)
OUT <- WORK_DIR
FIG  <- here::here("figures")
SUPP <- here::here("supplementary")
dir.create(FIG,  showWarnings = FALSE, recursive = TRUE)
dir.create(SUPP, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- 1) Locked signature
sig  <- readRDS(file.path(OUT, "signature_locked_latest.rds"))
scal <- sig$scaling
beta <- setNames(sig$coefficients$Coefficient, sig$coefficients$Feature)
cat("Coefficients:\n"); print(sig$coefficients)

# ---------------------------------------------------------------- 2) TRAIN
tr <- readRDS(file.path(OUT, "train_full_preprocessed.rds"))
X  <- tr$X_full
y  <- tr$y_full

if (is.null(rownames(X))) stop("X_full has no rownames - cannot match Age")
bars <- rownames(X)

# --- Observed Age from train_set.csv ---
ts <- readr::read_csv(file.path(OUT, "train_set.csv"), show_col_types = FALSE)
ts$Age_years <- ts$age_at_diagnosis / 365.25
cat(sprintf("train_set Age (years): min=%.1f median=%.1f max=%.1f\n",
            min(ts$Age_years, na.rm = TRUE),
            median(ts$Age_years, na.rm = TRUE),
            max(ts$Age_years, na.rm = TRUE)))

idx <- match(bars, ts$patient_barcode)
cat(sprintf("Barcode match: matched %d / %d | unmatched: %d\n",
            sum(!is.na(idx)), length(bars), sum(is.na(idx))))
if (any(is.na(idx))) {
  cat("Unmatched (first 5):",
      paste(head(bars[is.na(idx)], 5), collapse = ", "), "\n")
}

Age_vec <- ts$Age_years[idx]
if (any(is.na(Age_vec))) {
  Age_vec[is.na(Age_vec)] <- median(ts$Age_years, na.rm = TRUE)
  cat("[warn] Imputed", sum(is.na(idx)), "missing Age with train median.\n")
}

cat(sprintf("Age used -> mean=%.2f sd=%.2f  (locked: center=%.2f scale=%.2f)\n",
            mean(Age_vec), sd(Age_vec),
            scal$clinical$center["Age"], scal$clinical$scale["Age"]))

# --- z-scores with locked parameters ---
z1 <- (X[, "hsa-miR-196a-5p"] - scal$mirna$center["hsa-miR-196a-5p"]) /
  scal$mirna$scale ["hsa-miR-196a-5p"]
z2 <- (X[, "hsa-miR-675-5p"]  - scal$mirna$center["hsa-miR-675-5p"])  /
  scal$mirna$scale ["hsa-miR-675-5p"]
zA <- (Age_vec                - scal$clinical$center["Age"])          /
  scal$clinical$scale ["Age"]

cat(sprintf("\nSanity check z-scores:\n  z1: mean=%.3f sd=%.3f\n  z2: mean=%.3f sd=%.3f\n  zA: mean=%.3f sd=%.3f\n",
            mean(z1), sd(z1), mean(z2), sd(z2), mean(zA), sd(zA)))

lp <- beta["hsa-miR-196a-5p"] * z1 +
  beta["hsa-miR-675-5p"]  * z2 +
  beta["Age"]             * zA

train_df <- data.frame(
  sample_id = bars,
  score     = as.numeric(lp),
  OS_time   = as.numeric(y[, "time"]),
  OS_event  = as.integer(y[, "status"]),
  cohort    = "Train",
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------- 3) TEST + CGGA
d2 <- readr::read_csv(file.path(OUT, "supplementary_final", "DataS2_risk_scores.csv"),
                      show_col_types = FALSE)
cat("\nCohorts in DataS2:\n"); print(table(d2$cohort))

mk <- function(df, label) df %>% transmute(
  sample_id = patient_id,
  score     = as.numeric(risk_score),
  OS_time   = as.numeric(survival_time),
  OS_event  = as.integer(event),
  cohort    = label
)

test_df <- mk(d2 %>% filter(cohort == "TCGA_Test"), "Test")
cgga_df <- mk(d2 %>% filter(cohort == "CGGA"),      "CGGA")

clean <- function(d) d %>%
  filter(!is.na(score), !is.na(OS_time), !is.na(OS_event), OS_time > 0)
train_df <- clean(train_df)
test_df  <- clean(test_df)
cgga_df  <- clean(cgga_df)

cat(sprintf("\nFinal n -> Train: %d | Test: %d | CGGA: %d\n",
            nrow(train_df), nrow(test_df), nrow(cgga_df)))

# Check that the score definition is identical between Train (computed here) and DataS2
tr_ref <- d2 %>% filter(grepl("train", cohort, ignore.case = TRUE))
if (nrow(tr_ref) > 0) {
  m <- match(train_df$sample_id, tr_ref$patient_id); ok <- !is.na(m)
  cat(sprintf("Train score here vs DataS2: matched %d | cor=%.4f | max|diff|=%.4g\n",
              sum(ok), cor(train_df$score[ok], tr_ref$risk_score[m[ok]]),
              max(abs(train_df$score[ok] - tr_ref$risk_score[m[ok]]))))
} else cat("[info] DataS2 has no Train rows; cannot cross-check score definition.\n")

# Time unit and score scale across the three cohorts
for (nm in c("train_df","test_df","cgga_df")) {
  d <- get(nm)
  cat(sprintf("%-9s OS_time median=%.1f (range %.1f-%.1f) | score mean=%.3f sd=%.3f\n",
              nm, median(d$OS_time), min(d$OS_time), max(d$OS_time),
              mean(d$score), sd(d$score)))
}
if (median(cgga_df$OS_time) < 0.2 * median(train_df$OS_time) ||
    median(cgga_df$OS_time) > 5   * median(train_df$OS_time))
  warning("CGGA OS_time is not on the same scale as Train (days/months/years?) - KM/Cox/trim would be invalid.")
if (mean(cgga_df$score >= median(train_df$score)) > 0.9 ||
    mean(cgga_df$score >= median(train_df$score)) < 0.1)
  warning("More than 90% of CGGA falls on one side of the Train cutoff: score is probably not on the Train scale (batch/platform/z-score).")

# ---------------------------------------------------------------- 4) Cutoff from TRAIN only
cut_med <- median(train_df$score)
cut_ter <- quantile(train_df$score, probs = 2/3, names = FALSE)

cat("\n=== Score ranges ===\n")
for (nm in c("train_df","test_df","cgga_df")) {
  d <- get(nm)
  cat(sprintf("%-9s n=%3d  min=%.3f  median=%.3f  max=%.3f\n",
              nm, nrow(d), min(d$score), median(d$score), max(d$score)))
}
cat(sprintf("\nTrain median  cutoff = %.6f\n", cut_med))
cat(sprintf("Train tertile cutoff = %.6f\n", cut_ter))

# ---------------------------------------------------------------- 5) KM
assign_grp <- function(df, cv) {
  df$group <- factor(ifelse(df$score >= cv, "High", "Low"), levels = c("Low","High"))
  df
}

km_plot <- function(df, cv, ttl) {
  df <- assign_grp(df, cv)
  n_groups <- length(unique(df$group))
  cat(sprintf("  [%s] High=%d Low=%d groups=%d\n",
              df$cohort[1], sum(df$group == "High"),
              sum(df$group == "Low"), n_groups))
  if (n_groups < 2) {
    fit <- survival::survfit(survival::Surv(OS_time, OS_event) ~ 1, data = df)
    return(survminer::ggsurvplot(
      fit, data = df,
      pval = FALSE, risk.table = TRUE, conf.int = TRUE,
      palette = if (unique(df$group) == "High") "#D95F02" else "#2C7FB8",
      xlab = "Time (days)", ylab = "Overall survival probability",
      title = paste0(ttl, " [single group: ", as.character(unique(df$group)), "]"),
      ggtheme = theme_bw()
    ))
  }
  fit <- survival::survfit(survival::Surv(OS_time, OS_event) ~ group, data = df)
  survminer::ggsurvplot(
    fit, data = df,
    pval = TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE,
    legend.labs = c("Low","High"), legend.title = "Score group",
    palette = c("#2C7FB8","#D95F02"),
    xlab = "Time (days)", ylab = "Overall survival probability",
    title = ttl, ggtheme = theme_bw()
  )
}

cat("\n=== Assigning groups ===\n")
km_tr <- km_plot(train_df, cut_med, sprintf("Train (median cut = %.3f)", cut_med))
km_te <- km_plot(test_df,  cut_med, sprintf("Test  (same cut  = %.3f)", cut_med))
km_cg <- km_plot(cgga_df,  cut_med, sprintf("CGGA  (same cut  = %.3f)", cut_med))

pdf(file.path(FIG, "Fig_KM_train_cutoff.pdf"), width = 7.5, height = 8.5)
print(km_tr$plot); print(km_te$plot); print(km_cg$plot)
dev.off()

png(file.path(FIG, "Fig_KM_train_cutoff.png"),
    width = 7.5, height = 8.5, units = "in", res = 300)
print(km_tr$plot); print(km_te$plot); print(km_cg$plot)
dev.off()

# ---------------------------------------------------------------- 6) Statistics table (robust)
stats_one <- function(df, cohort, ctype, cv) {
  df <- assign_grp(df, cv)
  nH <- sum(df$group == "High")
  nL <- sum(df$group == "Low")

  if (nH == 0 || nL == 0) {
    return(data.frame(
      cohort = cohort, cutoff_type = ctype, cutoff_value = round(cv, 6),
      n = nrow(df), n_high = nH, n_low = nL,
      HR_high_vs_low = NA_real_, HR_lower_95 = NA_real_, HR_upper_95 = NA_real_,
      Cox_p = NA_real_, logrank_p = NA_real_,
      note = "single group after applying train cutoff",
      stringsAsFactors = FALSE))
  }

  sdf  <- survival::survdiff(survival::Surv(OS_time, OS_event) ~ group, data = df)
  p_lr <- 1 - stats::pchisq(sdf$chisq, df = length(sdf$n) - 1)   # an object named sd would mask sd()

  cox <- survival::coxph(survival::Surv(OS_time, OS_event) ~ group,
                         data = df, ties = "efron")
  s   <- summary(cox)

  coef_rows <- rownames(s$conf.int)
  hit <- grep("High", coef_rows, value = TRUE, ignore.case = TRUE)
  if (length(hit) == 0) {
    return(data.frame(
      cohort = cohort, cutoff_type = ctype, cutoff_value = round(cv, 6),
      n = nrow(df), n_high = nH, n_low = nL,
      HR_high_vs_low = NA_real_, HR_lower_95 = NA_real_, HR_upper_95 = NA_real_,
      Cox_p = NA_real_, logrank_p = signif(p_lr, 4),
      note = paste0("Cox coef name not found; rows=", paste(coef_rows, collapse = ",")),
      stringsAsFactors = FALSE))
  }
  rn <- hit[1]

  data.frame(
    cohort         = cohort,
    cutoff_type    = ctype,
    cutoff_value   = round(cv, 6),
    n              = nrow(df),
    n_high         = nH,
    n_low          = nL,
    HR_high_vs_low = round(s$conf.int[rn, "exp(coef)"], 3),
    HR_lower_95    = round(s$conf.int[rn, "lower .95"], 3),
    HR_upper_95    = round(s$conf.int[rn, "upper .95"], 3),
    Cox_p          = signif(s$coefficients[rn, "Pr(>|z|)"], 4),
    logrank_p      = signif(p_lr, 4),
    note           = if (min(nH, nL) < 10)
      sprintf("UNDERPOWERED: smallest group n=%d; HR/CI/p unreliable", min(nH, nL)) else "",
    stringsAsFactors = FALSE
  )
}

# Continuous effect (no cutoff): HR per 1 SD of score (SD from Train)
sd_train <- sd(train_df$score)
stats_cont <- function(df, cohort) {
  df$s1 <- df$score / sd_train
  cox <- survival::coxph(survival::Surv(OS_time, OS_event) ~ s1, data = df, ties = "efron")
  s <- summary(cox)
  data.frame(cohort = cohort, cutoff_type = "continuous_per_TrainSD", cutoff_value = NA_real_,
             n = nrow(df), n_high = NA_integer_, n_low = NA_integer_,
             HR_high_vs_low = round(s$conf.int[1, "exp(coef)"], 3),
             HR_lower_95 = round(s$conf.int[1, "lower .95"], 3),
             HR_upper_95 = round(s$conf.int[1, "upper .95"], 3),
             Cox_p = signif(s$coefficients[1, "Pr(>|z|)"], 4), logrank_p = NA_real_,
             note = "HR per 1 Train-SD of score; cutoff-free", stringsAsFactors = FALSE)
}

stats <- rbind(
  stats_one(train_df, "Train", "median",  cut_med),
  stats_one(test_df,  "Test",  "median",  cut_med),
  stats_one(cgga_df,  "CGGA",  "median",  cut_med),
  stats_one(train_df, "Train", "tertile", cut_ter),
  stats_one(test_df,  "Test",  "tertile", cut_ter),
  stats_one(cgga_df,  "CGGA",  "tertile", cut_ter),
  stats_cont(train_df, "Train"), stats_cont(test_df, "Test"), stats_cont(cgga_df, "CGGA"),
  # Exploratory: cohort-specific cutoff (to see whether the issue is scale only; NOT a replacement for the locked cutoff)
  stats_one(test_df, "Test", "median_cohort_specific(exploratory)", median(test_df$score)),
  stats_one(cgga_df, "CGGA", "median_cohort_specific(exploratory)", median(cgga_df$score))
)

readr::write_csv(stats, file.path(SUPP, "Table_S_KM.csv"))

cat("\n=== KM summary ===\n")
print(stats, row.names = FALSE)
cat("\nDone.\n  Figure:", file.path(FIG, "Fig_KM_train_cutoff.pdf / .png"),
    "\n  Table :", file.path(SUPP, "Table_S_KM.csv"), "\n")
