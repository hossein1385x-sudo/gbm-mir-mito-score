#!/usr/bin/env Rscript
# Rebuild Table S18 with gender


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({ library(data.table) })

input_dir <- paste0(WORK_DIR, "/")
supp_dir  <- file.path(input_dir, "supplementary_final")
norm_id   <- function(x) substr(toupper(gsub("\\.", "-", as.character(x))), 1, 12)

data_s1 <- fread(file.path(supp_dir, "DataS1_barcodes.csv"))

# merge with age from raw_clinical
clin <- as.data.table(readRDS(file.path(input_dir, "raw_clinical_cdr_integrated.rds")))
clin[, key := norm_id(patient_barcode)]

age_map <- unique(clin[, .(key, age_days = age_at_diagnosis)], by = "key")
age_map[, age_years := suppressWarnings(as.numeric(age_days)) / 365.25]

data_s1[, key := norm_id(patient_barcode_short)]
data_s1 <- merge(data_s1, age_map[, .(key, age_years)], by = "key", all.x = TRUE)

# summary per split
summarize <- function(dt) {
  data.frame(
    N = nrow(dt),
    Age_median = round(median(dt$age_years, na.rm = TRUE), 1),
    Age_IQR = sprintf("%.1f–%.1f",
                      quantile(dt$age_years, 0.25, na.rm = TRUE),
                      quantile(dt$age_years, 0.75, na.rm = TRUE)),
    Female_pct = round(100 * mean(dt$gender == "Female", na.rm = TRUE), 1),
    Male_pct   = round(100 * mean(dt$gender == "Male", na.rm = TRUE), 1),
    Dead_pct   = round(100 * mean(tolower(dt$vital_status) %in%
                                    c("dead", "deceased", "1"), na.rm = TRUE), 1),
    stringsAsFactors = FALSE)
}

s18 <- rbind(
  cbind(Cohort = "TCGA Training", summarize(data_s1[split == "train"])),
  cbind(Cohort = "TCGA Test",     summarize(data_s1[split == "test"]))
)

# add CGGA manually (if you want)
cat("\nTable S18:\n")
print(s18, row.names = FALSE)

write.csv(s18,
          file.path(supp_dir, "TableS18_baseline_characteristics.csv"),
          row.names = FALSE)
cat("\nTableS18 saved\n")
