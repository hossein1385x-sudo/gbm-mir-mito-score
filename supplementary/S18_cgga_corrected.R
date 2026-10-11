#!/usr/bin/env Rscript
# =========================================================
# fix_Table_S18_CGGA.R - corrected version
# =========================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(dplyr)
})

output_dir <- WORK_DIR

# ---------------------------------------------------------
# 1. Load the clinical data
# ---------------------------------------------------------
cgga_clin <- read.delim(
  file.path(output_dir, "cgga_raw",
            "CGGA.microRNA_array_198_clinical.20200506.txt"),
  stringsAsFactors = FALSE, check.names = FALSE)

cat("=== CGGA clinical columns ===\n")
print(colnames(cgga_clin))
cat("\n=== first 3 rows (sample) ===\n")
print(head(cgga_clin, 3))

# ---------------------------------------------------------
# 2. Load the list of 140 patients
# ---------------------------------------------------------
data_s2c <- read.csv(
  file.path(output_dir, "supplementary_final",
            "DataS2c_CGGA_analysis_data.csv"),
  stringsAsFactors = FALSE)

cat("\n=== DataS2c columns ===\n")
print(colnames(data_s2c))
cat("nrow:", nrow(data_s2c), "\n")

used_ids <- data_s2c[[1]]
cat("sample of used_ids:\n"); print(head(used_ids))

# ---------------------------------------------------------
# 3. Normalise the IDs
# ---------------------------------------------------------
normalize_id <- function(x) {
  x <- gsub("_", "-", as.character(x))
  x <- gsub(" ", "", x)
  toupper(x)
}

# find the ID column in clinical
possible_id_cols <- c("CGGA_ID", "PatientID", "patient_id",
                      "sample_id", "ID", "SampleID")
found_id_col <- possible_id_cols[possible_id_cols %in% colnames(cgga_clin)][1]

if (is.na(found_id_col)) {
  cat("available columns in clinical:\n")
  print(colnames(cgga_clin))
  stop("ID column not found")
}
cat("\nID column in clinical:", found_id_col, "\n")

# normalisation
used_ids_norm  <- normalize_id(used_ids)
cgga_clin$id_norm <- normalize_id(cgga_clin[[found_id_col]])

# ---------------------------------------------------------
# 4. Filter
# ---------------------------------------------------------
cgga_f <- cgga_clin[cgga_clin$id_norm %in% used_ids_norm, ]
cat("number matched:", nrow(cgga_f), "\n")

missing_in_clin <- setdiff(used_ids_norm, cgga_f$id_norm)
cat("IDs not in clinical:", length(missing_in_clin), "\n\n")

# ---------------------------------------------------------
# 5. Find the columns
# ---------------------------------------------------------
find_col <- function(df, candidates) {
  for (c in candidates) if (c %in% colnames(df)) return(c)
  return(NA)
}

col_age     <- find_col(cgga_f, c("Age", "age", "age_at_diagnosis"))
col_gender  <- find_col(cgga_f, c("Gender", "gender", "Sex", "sex"))
col_status  <- find_col(cgga_f, c("Censor", "censor", "VitalStatus",
                                  "vital_status", "Status", "status",
                                  "OS_status"))
col_os_time <- find_col(cgga_f, c("OS", "os", "survival_time",
                                  "OS_time", "follow_up", "FollowUp"))

cat("=== identified columns ===\n")
cat("Age     :", col_age,     "\n")
cat("Gender  :", col_gender,  "\n")
cat("Status  :", col_status,  "\n")
cat("OS_time :", col_os_time, "\n\n")

# ---------------------------------------------------------
# 6. Show the raw values
# ---------------------------------------------------------
cat("=== sample of raw values ===\n")
if (!is.na(col_age)) {
  cat("Age:\n"); print(head(cgga_f[[col_age]], 10))
}
if (!is.na(col_gender)) {
  cat("\nGender:\n"); print(table(cgga_f[[col_gender]], useNA = "ifany"))
}
if (!is.na(col_status)) {
  cat("\nStatus:\n"); print(table(cgga_f[[col_status]], useNA = "ifany"))
}
if (!is.na(col_os_time)) {
  cat("\nOS_time:\n"); print(head(cgga_f[[col_os_time]], 10))
}
cat("\n")

# ---------------------------------------------------------
# 7. Compute the statistics (with safe checks)
# ---------------------------------------------------------
# Age
if (!is.na(col_age)) {
  age <- suppressWarnings(as.numeric(cgga_f[[col_age]]))
} else {
  age <- rep(NA_real_, nrow(cgga_f))
}

# Gender
if (!is.na(col_gender)) {
  gender_clean <- tolower(trimws(as.character(cgga_f[[col_gender]])))
  is_male <- ifelse(grepl("^m", gender_clean), 1,
                    ifelse(grepl("^f", gender_clean), 0, NA))
} else {
  is_male <- rep(NA_real_, nrow(cgga_f))
}

# Status
if (!is.na(col_status)) {
  status_raw <- cgga_f[[col_status]]
  cat("=== status conversion ===\n")
  cat("type:", class(status_raw), "\n")
  cat("unique values:\n"); print(unique(status_raw))

  # convert to numeric
  if (is.character(status_raw) || is.factor(status_raw)) {
    status_num <- ifelse(
      tolower(as.character(status_raw)) %in%
        c("dead", "1", "deceased", "death", "died"),
      1, 0
    )
  } else {
    status_num <- as.numeric(status_raw)
  }

  cat("Dead count:", sum(status_num == 1, na.rm = TRUE),
      "/", length(status_num), "\n")
  cat("Dead pct:", round(100 * mean(status_num == 1, na.rm = TRUE), 1), "%\n\n")

  # if the dead percentage is below 30%, it is probably inverted
  dead_pct_raw <- mean(status_num == 1, na.rm = TRUE)
  if (!is.na(dead_pct_raw) && dead_pct_raw < 0.3) {
    cat(" dead_pct is below 30% — inverting\n")
    status_num <- 1 - status_num
  }
} else {
  status_num <- rep(NA_real_, nrow(cgga_f))
}

# OS time
if (!is.na(col_os_time)) {
  os <- suppressWarnings(as.numeric(cgga_f[[col_os_time]]))
} else {
  os <- rep(NA_real_, nrow(cgga_f))
}

# ---------------------------------------------------------
# 8. Final statistics
# ---------------------------------------------------------
N            <- nrow(cgga_f)
age_median   <- median(age, na.rm = TRUE)
age_iqr      <- quantile(age, c(0.25, 0.75), na.rm = TRUE)
male_pct     <- round(100 * mean(is_male == 1, na.rm = TRUE), 1)
dead_pct     <- round(100 * mean(status_num == 1, na.rm = TRUE), 1)
followup_med <- median(os, na.rm = TRUE)

cat("============================================================\n")
cat("Table S18 (CGGA row) — new values\n")
cat("============================================================\n")
cat(sprintf("N                 : %d\n", N))
cat(sprintf("Age_median        : %.1f\n", age_median))
cat(sprintf("Age_IQR           : %.1f–%.1f\n", age_iqr[1], age_iqr[2]))
cat(sprintf("Male_pct          : %.1f%%\n", male_pct))
cat(sprintf("Dead_pct          : %.1f%%\n", dead_pct))
cat(sprintf("Followup_median   : %.1f\n", followup_med))
cat("\n")

# ---------------------------------------------------------
# 9. Save
# ---------------------------------------------------------
table_s18_cgga <- data.frame(
  Cohort           = "CGGA",
  N                = N,
  Age_median       = round(age_median, 1),
  Age_IQR          = sprintf("%.1f–%.1f", age_iqr[1], age_iqr[2]),
  Male_pct         = male_pct,
  Dead_pct         = dead_pct,
  Followup_median  = round(followup_med, 1)
)

out_path <- file.path(output_dir, "supplementary_final",
                      "Table_S18_CGGA_corrected.csv")
write.csv(table_s18_cgga, out_path, row.names = FALSE)
cat("saved:", out_path, "\n")
