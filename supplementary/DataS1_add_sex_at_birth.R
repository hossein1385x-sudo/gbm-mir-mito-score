#!/usr/bin/env Rscript
# =============================================================================
# Merge sex_at_birth into DataS1
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({ library(data.table) })

supp_dir <- file.path(WORK_DIR, "supplementary_final")
norm_id  <- function(x) substr(toupper(gsub("\\.", "-", as.character(x))), 1, 12)

# --- 1. DataS1 ---
data_s1 <- fread(file.path(supp_dir, "DataS1_barcodes.csv"))
data_s1[, key := norm_id(patient_barcode_short)]

# --- 2. gender map from clin (obtained from TCGAbiolinks) ---
clin_dt <- as.data.table(clin)
clin_dt[, key := norm_id(bcr_patient_barcode)]
clin_dt[, gender_val := tolower(as.character(sex_at_birth))]

# normalise the values
clin_dt[gender_val == "female", gender_val := "Female"]
clin_dt[gender_val == "male",   gender_val := "Male"]
clin_dt[gender_val %in% c("unknown", "", "na", "NA"), gender_val := NA_character_]

gender_map <- unique(
  clin_dt[!is.na(gender_val), .(key, gender_val)],
  by = "key"
)

cat(sprintf("gender map: %d patients\n", nrow(gender_map)))
print(table(gender_map$gender_val))

# --- 3. drop the empty gender column and merge ---
if ("gender" %in% names(data_s1)) data_s1[, gender := NULL]
data_s1 <- merge(data_s1, gender_map, by = "key", all.x = TRUE)
setnames(data_s1, "gender_val", "gender")

# --- 4. check ---
n_ok   <- sum(!is.na(data_s1$gender))
n_test <- sum(!is.na(data_s1$gender) & data_s1$split == "test")

cat("\n============================================================\n")
cat("result\n")
cat("============================================================\n")
cat(sprintf("non-empty gender: %d of %d\n", n_ok, nrow(data_s1)))
cat(sprintf("  in test:       %d of %d\n", n_test,
            sum(data_s1$split == "test")))

cat("\ngender × split distribution:\n")
print(table(data_s1$gender, data_s1$split, useNA = "ifany"))

# --- 5. save ---
data_s1[, key := NULL]
fwrite(data_s1, file.path(supp_dir, "DataS1_barcodes.csv"))
cat("\nDataS1_barcodes.csv saved\n")
