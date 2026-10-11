#!/usr/bin/env Rscript
# =============================================================================
# Supplementary — TRIPOD+AI Checklist (prefilled)
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
})

out_dir <- file.path(WORK_DIR, "supplementary_final")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

tripod <- data.table(
  Section = c(
    rep("TITLE", 1),
    rep("ABSTRACT", 1),
    rep("INTRODUCTION", 3),
    rep("INTRODUCTION", 1),
    rep("METHODS", 2),
    rep("METHODS", 1),
    rep("METHODS", 1),
    rep("METHODS", 1),
    rep("METHODS", 1),
    rep("METHODS", 3),
    rep("METHODS", 5),
    rep("METHODS", 1),
    rep("METHODS", 1),
    rep("METHODS", 4),
    rep("METHODS", 3),
    rep("METHODS", 1),
    rep("METHODS", 1),
    rep("RESULTS", 3),
    rep("RESULTS", 2),
    rep("RESULTS", 2),
    rep("RESULTS", 1),
    rep("RESULTS", 1),
    rep("DISCUSSION", 1),
    rep("DISCUSSION", 1),
    rep("OPEN SCIENCE", 1)
  ),
  Item = c(
    "1", "2",
    "3a", "3b", "3c", "4",
    "5a", "5b",
    "6a", "6b",
    "7a", "7b",
    "8", "9",
    "10a", "10b", "10c",
    "11",
    "12a", "12b", "12c", "12d", "12e",
    "12f", "12g",
    "13", "14",
    "15a", "15b", "15c", "15d",
    "16a", "16b", "16c",
    "17", "18a", "18b",
    "19",
    "20a", "20b", "20c",
    "21", "22",
    "23a", "23b",
    "24", "25",
    "26",
    "27"
  ),
  Type = c(
    "D;E", "D;E",
    "D;E", "D;E", "D;E", "D;E",
    "D;E", "D;E",
    "D", "D",
    "D", "D",
    "D", "D",
    "D", "D", "D",
    "D",
    "D", "D", "D", "D", "D",
    "E", "D;E",
    "D", "D;E",
    "D;E", "D;E", "D;E", "D;E",
    "D", "D", "D",
    "D", "D", "D",
    "D;E",
    "D;E", "D;E", "E",
    "D;E", "D",
    "D;E", "D;E",
    "D;E", "D",
    "D;E",
    "D;E"
  ),
  Checklist_item = c(
    "Identify the study as developing or evaluating the performance of a multivariable prediction model, the target population, and the outcome to be predicted.",
    "See TRIPOD+AI for Abstracts checklist.",
    "Explain the healthcare context (including whether diagnostic or prognostic) and rationale for developing or evaluating the prediction model, including references to existing models.",
    "Describe the target population and the intended purpose of the prediction model in the context of the care pathway, including its intended users.",
    "Describe any known health inequalities between sociodemographic groups.",
    "Specify the study objectives, including whether the study describes the development or validation of a prediction model (or both).",
    "Describe the sources of data separately for the development and evaluation datasets, the rationale for using these data, and representativeness of the data.",
    "Specify the dates of data collection.",
    "Describe the eligibility criteria for study participants.",
    "Specify the outcome to be predicted by the model, including time horizon of predictions in case of prognostic models.",
    "List all predictors considered, including their definition and measurement.",
    "Describe how predictors were assessed, including whether assessors were blinded to outcome.",
    "Explain how the sample size was determined.",
    "Describe how missing data were handled.",
    "Describe data cleaning, preprocessing, and transformation.",
    "Describe any data augmentation or synthetic data generation.",
    "Describe how predictors were selected for inclusion in the model.",
    "Specify the type of model and method for fitting.",
    "Describe all analytical methods used for model development and evaluation.",
    "Describe any interactions considered.",
    "Describe any model updating (e.g., recalibration).",
    "Describe any methods used to address class imbalance.",
    "Describe any methods used to assess fairness across sociodemographic groups.",
    "Describe internal validation methods.",
    "Describe any other analytical methods.",
    "Describe the degree of class imbalance and how it was addressed.",
    "Report fairness assessment results across sociodemographic groups.",
    "Report model discrimination (e.g., C-index) in development and evaluation datasets.",
    "Report model calibration (e.g., calibration plot, Brier score).",
    "Report any other performance metrics.",
    "Report any clinical utility analyses.",
    "Describe the flow of participants through the study.",
    "Report participant flow diagram.",
    "Report data on non-participants.",
    "Report baseline characteristics of participants.",
    "Report the final model specification (coefficients, intercept, baseline).",
    "Report model performance in the evaluation dataset.",
    "Report fairness results in the evaluation dataset.",
    "Describe limitations of the study.",
    "Provide an overall interpretation of the results.",
    "Describe open science practices (data and code availability)."
  ),
  Reported_in_manuscript = c(
    "Title page (explicitly identifies development + external validation of prognostic model for GBM OS).",
    "Abstract — structured Background/Methods/Results/Conclusions.",
    "Introduction, paragraphs 1–3.",
    "Introduction, paragraph 4.",
    "Limitations — TCGA (Caucasian) and CGGA (Chinese) only; generalizability to other ancestries unverified.",
    "Introduction, final paragraph.",
    "Methods 1.1 (TCGA-GBM via GDC; CGGA via microRNA_array_198).",
    "Methods 1.1 (GDC data release; CGGA: 2020-05-06).",
    "Methods 1.1 (primary solid tumor; miRNA-seq ∩ RNA-seq ∩ complete survival).",
    "Methods 1.5–1.7 (OS; 6-month, 1-year, 2-year horizons).",
    "Methods 1.7 (6 miRNAs + age; from MDS-based DE analysis).",
    "Methods 1.2 (miRNA-seq; age from GDC clinical data).",
    "Methods 1.7 (167 train / 137 events; EPV ≥ 10; 1,000 bootstrap).",
    "Methods (complete-case analysis).",
    "Methods 1.5 (DESeq2 median-of-ratios; log2 transformation; z-score).",
    "N/A — no data augmentation used.",
    "Methods 1.5–1.7 (prefiltering; univariate Cox; Elastic Net-Cox).",
    "Methods 1.7 (Elastic Net-penalized Cox; λ via 5-fold CV).",
    "Methods 1.7, 1.12 (nested CV; bootstrap optimism correction).",
    "Methods 1.15 (exploratory risk × subtype interaction).",
    "N/A — no model updating performed.",
    "N/A — event rate > 70%; no class imbalance handling.",
    "Limitations — no fairness-specific analyses performed.",
    "Methods 1.7 (nested CV); Methods 1.23 (CGGA external).",
    "N/A — no additional analytical methods beyond sensitivity analyses.",
    "Methods — event rate 82%; no resampling.",
    "Limitations — not performed.",
    "Results 2.4–2.6 (C-index: 0.613 nested CV; 0.641 bootstrap; 0.660 test; 0.699 CGGA).",
    "Results 2.6–2.7, 2.12 (HR, AUC, calibration, Brier).",
    "Results 2.15 (DCA).",
    "Results 2.15 (DCA).",
    "Table S1 (617 clinical → 273 miRNA → 293 RNA → 229 intersection → 223 final).",
    "N/A — no flow diagram.",
    "N/A — not reported.",
    "Table S18 (age, sex, survival status, follow-up).",
    "Table S14 (coefficients, mean/SD, baseline survival, SHA-256).",
    "Results 2.4–2.6 (test set); Results 2.18 (CGGA).",
    "N/A — not performed.",
    "Limitations (extensive discussion).",
    "Discussion (biological and clinical interpretation).",
    "Data/Code Availability — Data S1–S3; code repository."
  )
)

# save
write.csv(tripod,
          file.path(out_dir, "TRIPOD_AI_checklist_prefilled.csv"),
          row.names = FALSE)

# short version (5 columns)
tripod_short <- tripod[, .(Section, Item, Type, Checklist_item, Reported_in_manuscript)]
write.csv(tripod_short,
          file.path(out_dir, "TRIPOD_AI_checklist_short.csv"),
          row.names = FALSE)

log_msg(sprintf("TRIPOD+AI checklist saved (%d items)", nrow(tripod)))
log_msg("   TRIPOD_AI_checklist_prefilled.csv")
log_msg("   TRIPOD_AI_checklist_short.csv")
