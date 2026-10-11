#!/usr/bin/env Rscript
# =============================================================================
# STEP 14.3 (Signature v2) — Verhaak Subtype Analysis
# -----------------------------------------------------------------------------
# Association of the v2 signature with the Verhaak mRNA subtype:
#   1) risk-score distribution by subtype
#   2) multivariable Cox (risk + subtype + Age)
#   3) likelihood-ratio test for incremental value
#   4) Interaction test (risk × subtype)
#   5) C-index within each subtype (no refit)
#
# Inputs:
#   signature_v2/signature_v2_test_predictions.csv
#   test_set_LOCKED.csv
#   verhaak_subtype_reference.csv  (downloaded via TCGAbiolinks if absent)
#
# Outputs in signature_v2/14_3_verhaak_analysis/
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(readr)
  library(stringr)
  library(broom)
  library(ggplot2)
  library(data.table)
})

output_dir <- paste0(WORK_DIR, "/")
sig_dir    <- file.path(output_dir, "signature_v2")
out_dir    <- file.path(sig_dir, "14_3_verhaak_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

set.seed(20240827)

# ======================================================================
# 1) Inputs
# ======================================================================
log_msg("============================================================")
log_msg("STEP 14.3 (Signature v2) — Verhaak Subtype Analysis")
log_msg("============================================================")
log_msg("loading ...")

path_scores     <- file.path(sig_dir, "signature_v2_test_predictions.csv")
path_clinical   <- file.path(output_dir, "test_set_LOCKED.csv")
path_verhaak    <- file.path(output_dir, "verhaak_subtype_reference.csv")

for (f in c(path_scores, path_clinical)) {
  if (!file.exists(f)) stop("not found: ", f)
}

test_scores   <- read_csv(path_scores,   show_col_types = FALSE)
test_clinical <- read_csv(path_clinical, show_col_types = FALSE)

log_msg(sprintf("   scores: %d rows", nrow(test_scores)))
log_msg(sprintf("   clinical: %d rows", nrow(test_clinical)))

# ======================================================================
# 2) Verhaak reference — if missing, take it from PanCancerAtlas
# ======================================================================
if (!file.exists(path_verhaak)) {
  log_msg("")
  log_msg("extracting the Verhaak subtype from PanCancerAtlas ...")

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("TCGAbiolinks", update = FALSE, ask = FALSE)
  }
  library(TCGAbiolinks)

  subtypes <- PanCancerAtlas_subtypes()
  req_cols <- c("cancer.type", "pan.samplesID", "Subtype_mRNA")
  miss <- setdiff(req_cols, names(subtypes))
  if (length(miss) > 0)
    stop("required columns are missing: ", paste(miss, collapse = ", "))

  verhaak_ref <- subtypes %>%
    filter(cancer.type == "GBM",
           !is.na(Subtype_mRNA),
           Subtype_mRNA != "") %>%
    mutate(patient_id = sapply(strsplit(pan.samplesID, "-"), function(p)
             paste(p[1:3], collapse = "-"))) %>%
    select(patient_id, verhaak_subtype) %>%
    distinct()

  write_csv(verhaak_ref, path_verhaak)
  log_msg(sprintf("   verhaak_subtype_reference.csv (%d records)", nrow(verhaak_ref)))
} else {
  log_msg("   verhaak_subtype_reference.csv already exists.")
}

verhaak_ref <- read_csv(path_verhaak, show_col_types = FALSE)

# ======================================================================
# 3) Identify the ID columns
# ======================================================================
find_id_col <- function(df) {
  cand <- c("patient_id", "patient_barcode", "sample_id",
            "case_id", "bcr_patient_barcode")
  hit <- intersect(cand, names(df))
  if (length(hit) == 0) stop("ID column not found.")
  hit[1]
}

id_scores   <- find_id_col(test_scores)
id_clinical <- find_id_col(test_clinical)
id_verhaak  <- find_id_col(verhaak_ref)

log_msg(sprintf("   id → scores: %s | clinical: %s | Verhaak: %s",
                id_scores, id_clinical, id_verhaak))

# standardise the names
setnames(as.data.table(test_scores), id_scores, "patient_id")
test_scores   <- as.data.frame(test_scores)
setnames(as.data.table(test_clinical), id_clinical, "patient_id")
test_clinical <- as.data.frame(test_clinical)
setnames(as.data.table(verhaak_ref), id_verhaak, "patient_id")
verhaak_ref   <- as.data.frame(verhaak_ref)

# normalise to TCGA-XX-XXXX
norm_id <- function(x) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  substr(x, 1, 12)
}
test_scores$patient_id   <- norm_id(test_scores$patient_id)
test_clinical$patient_id <- norm_id(test_clinical$patient_id)
verhaak_ref$patient_id   <- norm_id(verhaak_ref$patient_id)

# check for duplicates
for (nm in c("test_scores", "test_clinical", "verhaak_ref")) {
  df <- get(nm)
  if (anyDuplicated(df$patient_id) > 0)
    stop(sprintf("%s has duplicate IDs.", nm))
}

# ======================================================================
# 4) risk score
# ======================================================================
risk_col <- intersect(c("risk_score_full", "risk_score", "risk", "score",
                        "linear_predictor"),
                      names(test_scores))[1]
if (is.na(risk_col)) stop("risk score column not found.")

if (risk_col != "risk_score") {
  names(test_scores)[names(test_scores) == risk_col] <- "risk_score"
}

# ======================================================================
# 5) survival time and event
# ======================================================================
if ("survival_time" %in% names(test_scores) &&
    "event" %in% names(test_scores)) {
  test_scores$time  <- as.numeric(test_scores$survival_time)
  test_scores$event <- as.integer(test_scores$event)
} else {
  stop("columns survival_time and event are missing from the test predictions.")
}

# ======================================================================
# 6) Age
# ======================================================================
test_clinical$Age <- if ("age_at_diagnosis" %in% names(test_clinical)) {
  suppressWarnings(as.numeric(test_clinical$age_at_diagnosis))
} else if ("age" %in% names(test_clinical)) {
  suppressWarnings(as.numeric(test_clinical$age))
} else NA_real_

# convert days to years (if needed)
if (median(test_clinical$Age, na.rm = TRUE) > 200)
  test_clinical$Age <- test_clinical$Age / 365.25

# ======================================================================
# 7) Merge
# ======================================================================
df <- test_scores %>%
  select(patient_id, risk_score, time, event) %>%
  inner_join(test_clinical %>% select(patient_id, Age), by = "patient_id") %>%
  inner_join(verhaak_ref, by = "patient_id")

log_msg("")
log_msg(sprintf("merged: %d shared patients", nrow(df)))

df <- df %>%
  filter(!is.na(time), !is.na(event), !is.na(risk_score),
         !is.na(verhaak_subtype), time > 0)

df$verhaak_subtype <- factor(df$verhaak_subtype)

log_msg(sprintf("analysable: %d patients | %d events",
                nrow(df), sum(df$event == 1)))
log_msg("")
log_msg("Verhaak subtype distribution:")
print(table(df$verhaak_subtype, useNA = "ifany"))

# ======================================================================
# 8) Risk distribution by subtype
# ======================================================================
subtype_distribution <- df %>%
  group_by(verhaak_subtype) %>%
  summarise(
    n_patients = n(),
    n_events   = sum(event == 1),
    event_rate = mean(event == 1),
    median_risk = median(risk_score),
    IQR_risk    = IQR(risk_score),
    .groups = "drop"
  )
write_csv(subtype_distribution,
          file.path(out_dir, "subtype_distribution.csv"))

log_msg("")
log_msg("subtype distribution:")
print(subtype_distribution)

# Kruskal-Wallis
kw <- kruskal.test(risk_score ~ verhaak_subtype, data = df)
log_msg("")
log_msg(sprintf("Kruskal-Wallis: H=%.4f, df=%d, p=%.4g",
                as.numeric(kw$statistic), as.numeric(kw$parameter), kw$p.value))

# Boxplot
p_box <- ggplot(df, aes(x = verhaak_subtype, y = risk_score,
                        fill = verhaak_subtype)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  labs(title = "Risk score by Verhaak subtype",
       subtitle = sprintf("Kruskal-Wallis p = %.4g | n = %d",
                          kw$p.value, nrow(df)),
       x = "Verhaak mRNA subtype",
       y = "Locked Signature v2 Risk Score") +
  theme_minimal() + theme(legend.position = "none")

ggsave(file.path(out_dir, "boxplot_riskscore_by_subtype.png"),
       p_box, width = 7, height = 5, dpi = 600)

# ======================================================================
# 9) Multivariable Cox (exploratory)
# ======================================================================
log_msg("")
log_msg("Cox multivariable (risk + subtype + Age) ...")

cox_multi <- coxph(Surv(time, event) ~ scale(risk_score) +
                     verhaak_subtype + Age,
                   data = df, ties = "efron")
cox_multi_tidy <- tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE)
write_csv(cox_multi_tidy,
          file.path(out_dir, "cox_multivariable.csv"))

print(cox_multi_tidy)

# ======================================================================
# 10) LRT (base vs full)
# ======================================================================
log_msg("")
log_msg("LRT (incremental value) ...")

cox_base <- coxph(Surv(time, event) ~ verhaak_subtype + Age,
                  data = df, ties = "efron")

lrt <- anova(cox_base, cox_multi, test = "LRT")

lrt_p <- NA_real_
if ("Pr(>|Chi|)" %in% colnames(lrt)) {
  lrt_p <- lrt[2, "Pr(>|Chi|)"]
} else if ("P(>|Chi|)" %in% colnames(lrt)) {
  lrt_p <- lrt[2, "P(>|Chi|)"]
}

lrt_summary <- tibble(
  model_base = "verhaak_subtype + Age",
  model_full = "verhaak_subtype + Age + risk_score",
  lrt_chisq  = if ("Chisq" %in% colnames(lrt)) lrt$Chisq[2] else NA_real_,
  lrt_df     = if ("Df"    %in% colnames(lrt)) lrt$Df[2]    else NA_real_,
  lrt_p      = lrt_p,
  aic_base   = AIC(cox_base),
  aic_full   = AIC(cox_multi),
  delta_aic  = AIC(cox_multi) - AIC(cox_base)
)
write_csv(lrt_summary, file.path(out_dir, "LRT_incremental_value.csv"))
print(lrt_summary)

# ======================================================================
# 11) Interaction test (risk × subtype)
# ======================================================================
log_msg("")
log_msg("Interaction test (risk × subtype) ...")

subtype_counts    <- table(df$verhaak_subtype)
eligible_subtypes <- names(subtype_counts[subtype_counts >= 5])

interaction_summary <- if (length(eligible_subtypes) >= 2) {
  df_sub <- df %>% filter(verhaak_subtype %in% eligible_subtypes)
  cox_int <- coxph(Surv(time, event) ~ scale(risk_score) * verhaak_subtype,
                   data = df_sub, ties = "efron")
  lrt_int <- anova(cox_multi, cox_int, test = "LRT")
  p_int <- if ("Pr(>|Chi|)" %in% colnames(lrt_int)) lrt_int[2, "Pr(>|Chi|)"] else
           if ("P(>|Chi|)"  %in% colnames(lrt_int)) lrt_int[2, "P(>|Chi|)"]  else NA_real_

  tibble(
    lrt_chisq = lrt_int$Chisq[2],
    lrt_df    = lrt_int$Df[2],
    lrt_p     = p_int,
    eligible_subtypes = paste(eligible_subtypes, collapse = ", "),
    note = ifelse(!is.na(p_int) && p_int < 0.05,
                  "Evidence of subtype-dependent effect",
                  "No evidence of subtype-dependent effect")
  )
} else {
  tibble(lrt_chisq = NA_real_, lrt_df = NA_real_, lrt_p = NA_real_,
         eligible_subtypes = paste(eligible_subtypes, collapse = ", "),
         note = "Not performed: fewer than 2 subtypes with >=5 patients")
}
write_csv(interaction_summary,
          file.path(out_dir, "interaction_test_subtype.csv"))
print(interaction_summary)

# ======================================================================
# 12) C-index within each subtype (no refit)
# ======================================================================
log_msg("")
log_msg("C-index within each subtype ...")

subtype_cindex <- df %>%
  group_by(verhaak_subtype) %>%
  group_modify(~ {
    n_ev <- sum(.x$event == 1)
    if (n_ev < 15) {
      return(tibble(n_patients = nrow(.x), n_events = n_ev,
                    cindex = NA_real_, se = NA_real_,
                    note = sprintf("Not estimated: %d events (<15)", n_ev)))
    }
    obj <- concordance(Surv(time, event) ~ risk_score, data = .x, reverse = TRUE)
    tibble(n_patients = nrow(.x), n_events = n_ev,
           cindex = as.numeric(obj$concordance),
           se     = as.numeric(obj$std.err),
           note   = "Locked risk score; no refit")
  }) %>% ungroup()

write_csv(subtype_cindex,
          file.path(out_dir, "subtype_stratified_cindex.csv"))
print(subtype_cindex)

# ======================================================================
# 13) Overall C-index (cross-check with Step 10)
# ======================================================================
overall <- concordance(Surv(time, event) ~ risk_score, data = df, reverse = TRUE)
overall_cindex <- as.numeric(overall$concordance)
overall_se     <- as.numeric(overall$std.err)

overall_summary <- tibble(
  n_patients = nrow(df),
  n_events   = sum(df$event == 1),
  cindex     = overall_cindex,
  se         = overall_se,
  source     = "Locked risk score (cross-check vs Step 10)"
)
write_csv(overall_summary,
          file.path(out_dir, "overall_cindex_crosscheck.csv"))

# ======================================================================
# 14) HR of risk_score (adjusted)
# ======================================================================
risk_row <- cox_multi_tidy %>% filter(grepl("risk_score", term))

risk_hr_summary <- if (nrow(risk_row) == 1) {
  tibble(HR = risk_row$estimate,
         CI_low  = risk_row$conf.low,
         CI_high = risk_row$conf.high,
         p_value = risk_row$p.value)
} else {
  tibble(HR = NA_real_, CI_low = NA_real_,
         CI_high = NA_real_, p_value = NA_real_)
}
write_csv(risk_hr_summary,
          file.path(out_dir, "risk_score_HR_adjusted.csv"))

# ======================================================================
# 15) Summary
# ======================================================================
log_msg("")
log_msg("============================================================")
log_msg("Summary of STEP 14.3 (Signature v2)")
log_msg("============================================================")
log_msg(sprintf("n = %d | events = %d", nrow(df), sum(df$event == 1)))
log_msg(sprintf("C-index overall: %.4f (SE=%.4f)", overall_cindex, overall_se))
log_msg(sprintf("Kruskal-Wallis (risk ~ subtype): p = %.4g", kw$p.value))
log_msg(sprintf("LRT (incremental value): p = %.4g | ΔAIC = %.2f",
                lrt_summary$lrt_p, lrt_summary$delta_aic))
log_msg(sprintf("Interaction (risk × subtype): p = %.4g",
                interaction_summary$lrt_p))
log_msg("")
log_msg("C-index by subtype:")
for (i in seq_len(nrow(subtype_cindex))) {
  log_msg(sprintf("   %-20s n=%d ev=%d C=%.4f",
                  subtype_cindex$verhaak_subtype[i],
                  subtype_cindex$n_patients[i],
                  subtype_cindex$n_events[i],
                  subtype_cindex$cindex[i]))
}
log_msg("")
log_msg("outputs in: ", out_dir)
log_msg("STEP 14.3 (Signature v2) COMPLETED")
