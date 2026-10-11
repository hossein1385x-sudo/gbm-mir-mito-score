#!/usr/bin/env Rscript
# =============================================================================
# Covariate Audit + One-at-a-Time Bivariate Adjustment
# -----------------------------------------------------------------------------
# Why this script:
#   MGMT status is usable in only 6-7 patients. Building a multivariable model
#   with all covariates would reduce the sample intersection to zero.
#   The appropriate design:
#     1) a separate audit of missingness for each variable in our cohort (223 patients)
#     2) a per-variable decision (EXCLUDE / ONE-AT-A-TIME / JOINT)
#     3) for ONE-AT-A-TIME: a bivariate model (risk score + age + that variable)
#     4) automatic generation of the Limitations text
#
# Inputs:
#   signature_v2/signature_v2_test_predictions.csv   (risk score test)
#   supplementary_final/DataS1_barcodes.csv          (223-patient cohort)
#   raw_clinical_cdr_integrated.rds                  (base clinical data)
#   test_set_with_MGMT.csv                           (available MGMT)
#   any other file that may contain covariates
#
# Outputs:
#   supplementary_final/covariate_audit.csv
#   supplementary_final/bivariate_adjustment_test.csv
#   supplementary_final/bivariate_adjustment_cgga.csv
#   supplementary_final/limitations_covariate_text.txt
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

input_dir <- paste0(WORK_DIR, "/")
sig_dir   <- file.path(input_dir, "signature_v2")
supp_dir  <- file.path(input_dir, "supplementary_final")
out_dir   <- supp_dir

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

norm_id <- function(x) substr(toupper(gsub("\\.", "-", as.character(x))), 1, 12)

# =============================================================================
# 0) Candidate sources for covariates
# -----------------------------------------------------------------------------
# Important: the covariates below may be stored in different GDC files
#   - MGMT: TCGA clinical supplement (biotab) → sometimes a separate file
#   - KPS: follow_up biotab → usually not available via the GDC /cases endpoint
#   - Extent of resection: diagnosis / follow_up biotab
#   - Treatment: drug / radiation biotab
# =============================================================================

candidate_files <- list(
  list(path = file.path(input_dir, "raw_clinical_cdr_integrated.rds"),
       desc = "GDC clinical (base)"),
  list(path = file.path(input_dir, "test_set_with_MGMT.csv"),
       desc = "MGMT supplemental"),
  list(path = file.path(supp_dir, "DataS1_barcodes.csv"),
       desc = "DataS1 (cohort 223)"),
  list(path = file.path(sig_dir, "signature_v2_test_predictions.csv"),
       desc = "test predictions (v2)"),
  # Other possible files — skipped if they do not exist
  list(path = file.path(input_dir, "clinical_biotab_merged.csv"),
       desc = "merged clinical biotabs"),
  list(path = file.path(input_dir, "followup_merged.csv"),
       desc = "follow-up biotab"),
  list(path = file.path(input_dir, "treatment_merged.csv"),
       desc = "treatment biotab")
)

load_any <- function(path, desc) {
  if (!file.exists(path)) {
    log_msg(sprintf("    not found: %s", basename(path)))
    return(NULL)
  }
  x <- tryCatch({
    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
      as.data.table(readRDS(path))
    } else {
      fread(path)
    }
  }, error = function(e) {
    log_msg(sprintf("    error reading %s: %s", basename(path), conditionMessage(e)))
    NULL
  })
  if (!is.null(x)) log_msg(sprintf("   %s (%d × %d)", desc, nrow(x), ncol(x)))
  x
}

log_msg("============================================================")
log_msg("STEP 0 — loading covariate sources")
log_msg("============================================================")

sources <- lapply(candidate_files, function(f) {
  list(dt = load_any(f$path, f$desc), desc = f$desc, path = f$path)
})
names(sources) <- sapply(sources, `[[`, "path")

# =============================================================================
# 1) Identify the ID column in each source
# =============================================================================
id_candidates <- c("patient_barcode", "patient_id", "bcr_patient_barcode",
                   "submitter_id", "patient_barcode_short",
                   "patient_barcode_full", "CGGA_ID", "case_id")

find_id_col <- function(dt) {
  hit <- intersect(id_candidates, names(dt))
  if (length(hit) == 0) return(NULL)
  hit[1]
}

# Keep only sources that have a usable ID column
usable_sources <- list()
for (i in seq_along(sources)) {
  dt <- sources[[i]]$dt
  if (is.null(dt) || nrow(dt) == 0) next
  idc <- find_id_col(dt)
  if (is.null(idc)) {
    log_msg(sprintf("    no ID column found in %s.", sources[[i]]$desc))
    next
  }
  dt[, .key_short := norm_id(get(idc))]
  usable_sources[[length(usable_sources) + 1]] <- list(
    dt = dt, desc = sources[[i]]$desc, id_col = idc)
}

log_msg("")
log_msg(sprintf("usable sources: %d", length(usable_sources)))

# =============================================================================
# 2) Reference cohort (223 patients) from DataS1
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("STEP 1 — building the reference cohort")
log_msg("============================================================")

data_s1_path <- file.path(supp_dir, "DataS1_barcodes.csv")
if (!file.exists(data_s1_path)) {
  stop("DataS1_barcodes.csv not found. Run the DataS1 script first.")
}

data_s1 <- fread(data_s1_path)
data_s1[, .key_short := norm_id(patient_barcode_short)]

cohort_test <- data_s1[split == "test"]
n_test <- nrow(cohort_test)
log_msg(sprintf("   test cohort: %d patients", n_test))

# =============================================================================
# 3) List of candidate covariates
# =============================================================================
candidate_vars <- c(
  # Demographic / clinical
  "age_at_diagnosis", "age", "gender", "sex",
  # Performance status
  "karnofsky_performance_score", "kps", "karnofsky_score",
  "eastern_cooperative_oncology_group", "ecog",
  # Tumor characteristics
  "tumor_grade", "grade", "primary_diagnosis", "histology",
  "idh_status", "idh_mutation_status", "idh1_mutation",
  "mgmt_methylation", "mgmt_status", "mgmt_promoter_methylation",
  "1p19q_codeletion_status",
  # Treatment
  "extent_of_resection", "extent_of_tumor_resection",
  "treatment_type", "radiation_therapy", "pharmaceutical_therapy",
  # Molecular subtypes
  "verhaak_subtype", "subtype_mrna", "molecular_subtype",
  # Clinical scores
  "recursive_partitioning_analysis", "rpa_class"
)

# =============================================================================
# 4) Merge all sources into a single table (type-safe version)
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("STEP 2 — merging covariate sources")
log_msg("============================================================")

# Start from the reference cohort
master <- unique(data_s1[, .(.key_short,
                             patient_barcode_short,
                             patient_barcode_full,
                             split, vital_status)], by = ".key_short")

# Helper: first non-NA / non-empty value per key
first_non_na <- function(dt_sub, cols) {
  # dt_sub is a data.table with columns .key_short and cols
  # convert all cols to character so that the types are uniform
  for (col in cols) {
    if (!is.character(dt_sub[[col]])) {
      dt_sub[, (col) := as.character(get(col))]
    }
  }
  # drop empty values or "NA"
  dt_sub[get(cols[1]) == "" | get(cols[1]) == "NA", (cols[1]) := NA_character_]

  out <- dt_sub[, lapply(.SD, function(x) {
    x <- x[!is.na(x) & x != "" & x != "NA"]
    if (length(x) == 0) NA_character_ else x[1]
  }), by = ".key_short"]

  out
}

for (src in usable_sources) {
  dt <- src$dt

  # only candidate columns that exist in this source
  cols_here <- intersect(candidate_vars, names(dt))
  if (length(cols_here) == 0) {
    log_msg(sprintf("   [%s] has no known covariates.", src$desc))
    next
  }

  log_msg(sprintf("   [%s] adding %d columns: %s",
                  src$desc, length(cols_here),
                  paste(cols_here, collapse = ", ")))

  # only the required columns
  sub <- dt[, c(".key_short", cols_here), with = FALSE]

  # for each key, the first non-NA value in each column
  sub <- first_non_na(sub, cols_here)

  # merge with master — each column is merged separately
  for (col in cols_here) {
    if (col %in% names(master)) {
      # if already present, only fill in the gaps
      master <- merge(master,
                      sub[, c(".key_short", col), with = FALSE],
                      by = ".key_short", all.x = TRUE,
                      suffixes = c("", ".new"))
      new_col <- paste0(col, ".new")
      if (new_col %in% names(master)) {
        empty_mask <- is.na(master[[col]]) | master[[col]] == "" |
          master[[col]] == "NA"
        master[empty_mask, (col) := get(new_col)]
        master[, (new_col) := NULL]
      }
    } else {
      master <- merge(master,
                      sub[, c(".key_short", col), with = FALSE],
                      by = ".key_short", all.x = TRUE)
    }
  }
}

log_msg(sprintf("   master: %d × %d", nrow(master), ncol(master)))

# =============================================================================
# 5) Audit — check availability of each covariate in all of TCGA and in our cohort
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("STEP 3 — Audit availability")
log_msg("============================================================")

# all candidate columns that appeared in master
vars_in_master <- intersect(candidate_vars, names(master))

audit_rows <- lapply(vars_in_master, function(v) {
  x <- master[[v]]
  x <- as.character(x)
  available <- !is.na(x) & x != "" & x != "NA"
  n_avail <- sum(available)

  # check in the test cohort
  is_test <- master$split == "test"
  n_avail_test <- sum(available & is_test)

  data.frame(
    variable = v,
    n_available_in_cohort = n_avail,
    pct_in_cohort = round(100 * n_avail / nrow(master), 1),
    n_available_in_test = n_avail_test,
    pct_in_test = round(100 * n_avail_test / sum(is_test), 1),
    stringsAsFactors = FALSE)
})

audit <- rbindlist(audit_rows)
audit <- audit[order(-n_available_in_test)]
setcolorder(audit, c("variable", "n_available_in_cohort", "pct_in_cohort",
                     "n_available_in_test", "pct_in_test"))

log_msg("")
log_msg("Audit results (sorted by availability in the test set):")
print(audit, row.names = FALSE)

# =============================================================================
# 6) Decision tree
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("STEP 4 — per-variable decision")
log_msg("============================================================")

decide_strategy <- function(n_test, threshold_exclude = 15, threshold_joint = 50) {
  if (is.na(n_test) || n_test < threshold_exclude) {
    return("EXCLUDE")
  }
  pct <- 100 * n_test / 56  # n_test out of 56 patients
  if (n_test >= threshold_exclude & pct < threshold_joint) {
    return("ONE-AT-A-TIME")
  }
  return("JOINT-ELIGIBLE")
}

audit[, strategy := sapply(n_available_in_test, decide_strategy)]

log_msg("")
log_msg("strategy per variable:")
print(audit[, .(variable, n_available_in_test, pct_in_test, strategy)],
      row.names = FALSE)

# =============================================================================
# 7) Bivariate adjustment for ONE-AT-A-TIME variables
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("STEP 5 — Bivariate adjustment (one-at-a-time)")
log_msg("============================================================")

# --- merge with the risk score ---
pred_path <- file.path(sig_dir, "signature_v2_test_predictions.csv")
if (!file.exists(pred_path)) stop("signature_v2_test_predictions.csv is missing.")
pred <- fread(pred_path)
pred[, .key_short := norm_id(patient_id)]
pred_sub <- pred[, .(.key_short,
                     risk_score_full = risk_score_full,
                     risk_score_mirna = risk_score_mirna,
                     survival_time, event)]

df_test <- merge(master[split == "test"],
                 pred_sub,
                 by = ".key_short",
                 all.x = TRUE)

df_test <- df_test[!is.na(survival_time) & !is.na(event) &
                   !is.na(risk_score_full) & survival_time > 0]

log_msg(sprintf("   test cohort for bivariate: %d patients", nrow(df_test)))

# --- age from DataS1 ---
if ("age_at_diagnosis" %in% names(master)) {
  age_raw <- suppressWarnings(as.numeric(df_test$age_at_diagnosis))
  if (median(age_raw, na.rm = TRUE) > 200) age_raw <- age_raw / 365.25
  df_test[, age_num := age_raw]
} else {
  df_test[, age_num := NA_real_]
}

# --- bivariate adjustment function ---
run_bivariate <- function(var_name, data, cohort_label = "test") {

  if (!var_name %in% names(data)) {
    return(data.table(
      variable = var_name, cohort = cohort_label,
      n_used = 0L, n_events = 0L,
      HR_risk = NA_real_, HR_low = NA_real_, HR_high = NA_real_,
      p_value = NA_real_, note = "variable missing"))
  }

  sub <- data[!is.na(get(var_name)) & get(var_name) != "" &
              !is.na(survival_time) & !is.na(event) &
              !is.na(risk_score_full) & survival_time > 0]

  n_used <- nrow(sub)
  n_events <- sum(sub$event == 1, na.rm = TRUE)

  if (n_used < 10) {
    return(data.table(
      variable = var_name, cohort = cohort_label,
      n_used = n_used, n_events = n_events,
      HR_risk = NA_real_, HR_low = NA_real_, HR_high = NA_real_,
      p_value = NA_real_,
      note = sprintf("n=%d < 10 — unstable model", n_used)))
  }

  # convert the variable to factor or numeric depending on its type
  x <- sub[[var_name]]
  if (is.numeric(x)) {
    sub[, covar := scale(x)]
  } else {
    sub[, covar := as.factor(x)]
    # if more than 5 levels or any level has < 3 samples → collapse to binary only
    tbl <- table(sub$covar)
    if (length(tbl) > 5 || any(tbl < 3)) {
      return(data.table(
        variable = var_name, cohort = cohort_label,
        n_used = n_used, n_events = n_events,
        HR_risk = NA_real_, HR_low = NA_real_, HR_high = NA_real_,
        p_value = NA_real_,
        note = sprintf("categorical with %d levels or sparsely populated levels", length(tbl))))
    }
  }

  # formula
  if (all(is.na(sub$age_num))) {
    form <- Surv(survival_time, event) ~ scale(risk_score_full) + covar
  } else {
    form <- Surv(survival_time, event) ~ scale(risk_score_full) + covar + scale(age_num)
  }

  fit <- tryCatch(coxph(form, data = sub),
                  error = function(e) NULL)

  if (is.null(fit)) {
    return(data.table(
      variable = var_name, cohort = cohort_label,
      n_used = n_used, n_events = n_events,
      HR_risk = NA_real_, HR_low = NA_real_, HR_high = NA_real_,
      p_value = NA_real_, note = "coxph failed"))
  }

  s <- summary(fit)
  coef_name <- grep("risk_score_full", rownames(s$coefficients), value = TRUE)[1]
  if (is.na(coef_name)) {
    return(data.table(
      variable = var_name, cohort = cohort_label,
      n_used = n_used, n_events = n_events,
      HR_risk = NA_real_, HR_low = NA_real_, HR_high = NA_real_,
      p_value = NA_real_, note = "risk_score coef not found"))
  }

  b  <- s$coefficients[coef_name, "coef"]
  se <- s$coefficients[coef_name, "se(coef)"]
  pv <- s$coefficients[coef_name, "Pr(>|z|)"]

  data.table(
    variable = var_name, cohort = cohort_label,
    n_used = n_used, n_events = n_events,
    HR_risk = round(exp(b), 3),
    HR_low  = round(exp(b - 1.96 * se), 3),
    HR_high = round(exp(b + 1.96 * se), 3),
    p_value = signif(pv, 3),
    note = "OK")
}

# --- run for all ONE-AT-A-TIME variables ---
vars_one_at_a_time <- audit[strategy == "ONE-AT-A-TIME", variable]
log_msg(sprintf("   ONE-AT-A-TIME variables: %d", length(vars_one_at_a_time)))

if (length(vars_one_at_a_time) > 0) {
  bivariate_test <- rbindlist(lapply(vars_one_at_a_time,
                                     run_bivariate, data = df_test,
                                     cohort_label = "TCGA_test"))
  log_msg("")
  log_msg("Bivariate adjustment results (test set):")
  print(bivariate_test, row.names = FALSE)
} else {
  bivariate_test <- data.table()
  log_msg("   no ONE-AT-A-TIME variables.")
}

# =============================================================================
# 8) Bivariate adjustment in CGGA (optional, for external validation)
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("STEP 6 — bivariate adjustment in CGGA")
log_msg("============================================================")

cgga_path <- file.path(supp_dir, "CGGA_sample_scores.csv")
bivariate_cgga <- data.table()

if (!file.exists(cgga_path)) {
  log_msg("   CGGA_sample_scores.csv is missing — skipping.")
} else {
  cgga <- fread(cgga_path)

  # add the covariates available in CGGA
  cgga_clin_file <- file.path(input_dir,
    "cgga_raw/clin cgga/CGGA.microRNA_array_198_clinical.20200506.txt")

  if (!file.exists(cgga_clin_file)) {
    log_msg("   CGGA clinical file is missing — skipping.")
  } else {
    cgga_clin <- fread(cgga_clin_file)
    cgga_clin[, .key_short := norm_id(CGGA_ID)]

    # CGGA covariates
    cgga_vars <- intersect(c("Grade", "IDH_mutation_status",
                             "MGMTp_methylation_status", "PRS_type",
                             "Radio_status (treated=1;un-treated=0)",
                             "Chemo_status (TMZ treated=1;un-treated=0)",
                             "Gender", "Age"),
                           names(cgga_clin))

    log_msg(sprintf("   CGGA covariates: %s",
                    paste(cgga_vars, collapse = ", ")))

    cgga_merge <- merge(cgga,
                        cgga_clin[, c(".key_short", cgga_vars), with = FALSE],
                        by.x = "CGGA_ID", by.y = ".key_short", all.x = TRUE)

    # name the columns to match run_bivariate
    cgga_merge[, risk_score_full := risk_score]
    cgga_merge[, survival_time := OS]
    cgga_merge[, event := status]

    # one bivariate model per CGGA variable
    cgga_results <- rbindlist(lapply(cgga_vars, function(v) {
      run_bivariate(v, data = cgga_merge, cohort_label = "CGGA")
    }))

    bivariate_cgga <- cgga_results
    log_msg("")
    log_msg("Bivariate adjustment results (CGGA):")
    print(bivariate_cgga, row.names = FALSE)
  }
}

# =============================================================================
# 9) Automatic generation of the Limitations text
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("STEP 7 — generating the Limitations text")
log_msg("============================================================")

generate_limitation_text <- function(audit, bivariate_test, bivariate_cgga) {

  excluded <- audit[strategy == "EXCLUDE"]
  one_at_time <- audit[strategy == "ONE-AT-A-TIME"]

  parts <- c()

  # --- EXCLUDE section ---
  if (nrow(excluded) > 0) {
    ex_sentences <- sapply(seq_len(nrow(excluded)), function(i) {
      sprintf("%s was available for only %d of 56 test-set patients (%.1f%%)",
              excluded$variable[i],
              excluded$n_available_in_test[i],
              excluded$pct_in_test[i])
    })
    parts <- c(parts, sprintf(
      "Several candidate clinical covariates were excluded because of extensive missingness in the test cohort: %s. These variables were therefore not included in any adjustment model, and their potential confounding effect on the signature could not be evaluated.",
      paste(ex_sentences, collapse = "; ")))
  }

  # --- ONE-AT-A-TIME section ---
  if (nrow(one_at_time) > 0) {
    oa_sentences <- sapply(seq_len(nrow(one_at_time)), function(i) {
      var <- one_at_time$variable[i]
      n   <- one_at_time$n_available_in_test[i]
      res <- bivariate_test[variable == var]

      if (nrow(res) == 0 || is.na(res$HR_risk[1])) {
        sprintf("%s (n=%d, adjustment not performed)", var, n)
      } else {
        sprintf("%s (n=%d; adjusted HR = %.2f, 95%% CI %.2f\u2013%.2f, P = %.3g)",
                var, n, res$HR_risk[1], res$HR_low[1],
                res$HR_high[1], res$p_value[1])
      }
    })
    parts <- c(parts, sprintf(
      "For covariates with moderate availability, one-at-a-time bivariate adjustment was performed rather than a single multivariable model, to preserve the maximum available sample size for each comparison: %s. This approach avoids the loss of power that would result from requiring complete data across all covariates simultaneously.",
      paste(oa_sentences, collapse = "; ")))
  }

  # --- CGGA section ---
  if (nrow(bivariate_cgga) > 0) {
    cg_sentences <- sapply(seq_len(nrow(bivariate_cgga)), function(i) {
      sprintf("%s (n=%d; HR = %.2f, 95%% CI %.2f\u2013%.2f, P = %.3g)",
              bivariate_cgga$variable[i],
              bivariate_cgga$n_used[i],
              bivariate_cgga$HR_risk[i],
              bivariate_cgga$HR_low[i],
              bivariate_cgga$HR_high[i],
              bivariate_cgga$p_value[i])
    })
    parts <- c(parts, sprintf(
      "In the CGGA external cohort, bivariate adjustment for the most clinically relevant covariates was feasible: %s.",
      paste(cg_sentences, collapse = "; ")))
  }

  paste(parts, collapse = "\n\n")
}

limitation_text <- generate_limitation_text(audit, bivariate_test, bivariate_cgga)

log_msg("")
log_msg("Limitations text generated:")
log_msg("")
cat(limitation_text)
cat("\n")

# =============================================================================
# 10) Saving the outputs
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("STEP 8 — saving the outputs")
log_msg("============================================================")

fwrite(audit,
       file.path(out_dir, "covariate_audit.csv"))
fwrite(bivariate_test,
       file.path(out_dir, "bivariate_adjustment_test.csv"))
if (nrow(bivariate_cgga) > 0) {
  fwrite(bivariate_cgga,
         file.path(out_dir, "bivariate_adjustment_cgga.csv"))
}
writeLines(limitation_text,
           file.path(out_dir, "limitations_covariate_text.txt"))

# --- summary ---
summary_lines <- c(
  "============================================================",
  "Covariate Audit Summary",
  "============================================================",
  sprintf("Total cohort: %d patients", nrow(master)),
  sprintf("Test cohort:  %d patients", n_test),
  "",
  "Strategy distribution:",
  sprintf("  EXCLUDE:         %d", sum(audit$strategy == "EXCLUDE")),
  sprintf("  ONE-AT-A-TIME:   %d", sum(audit$strategy == "ONE-AT-A-TIME")),
  sprintf("  JOINT-ELIGIBLE:  %d", sum(audit$strategy == "JOINT-ELIGIBLE")),
  "",
  "Excluded variables:",
  paste(sprintf("  - %s (n=%d)",
                audit[strategy == "EXCLUDE", variable],
                audit[strategy == "EXCLUDE", n_available_in_test]),
        collapse = "\n"),
  "",
  "One-at-a-time variables:",
  paste(sprintf("  - %s (n=%d)",
                audit[strategy == "ONE-AT-A-TIME", variable],
                audit[strategy == "ONE-AT-A-TIME", n_available_in_test]),
        collapse = "\n")
)
writeLines(summary_lines,
           file.path(out_dir, "covariate_audit_summary.txt"))

log_msg("")
log_msg("outputs:")
log_msg("   covariate_audit.csv")
log_msg("   bivariate_adjustment_test.csv")
if (nrow(bivariate_cgga) > 0) log_msg("   bivariate_adjustment_cgga.csv")
log_msg("   limitations_covariate_text.txt")
log_msg("   covariate_audit_summary.txt")
log_msg("")
log_msg("DONE")
