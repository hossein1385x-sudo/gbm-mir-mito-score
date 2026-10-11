#!/usr/bin/env Rscript
# =============================================================================
# Data S1 — Complete patient-level metadata and split assignment
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
})

input_dir <- paste0(WORK_DIR, "/")
sig_dir   <- file.path(input_dir, "signature_v2")
out_dir   <- file.path(input_dir, "supplementary_final")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

norm_id <- function(x) substr(toupper(gsub("\\.", "-", as.character(x))), 1, 12)

# =============================================================================
# 1) Load the sources
# =============================================================================
log_msg("loading the sources ...")

# --- clinical ---
clin_raw <- readRDS(file.path(input_dir, "raw_clinical.rds"))
if (!is.data.table(clin_raw)) clin_raw <- as.data.table(clin_raw)
clin_raw[, key_short := norm_id(patient_barcode)]

# --- miRNA sample IDs ---
mirna_raw <- readRDS(file.path(input_dir, "raw_mirna.rds"))
if (!is.data.table(mirna_raw)) mirna_raw <- as.data.table(mirna_raw)

mirna_meta <- unique(mirna_raw[, .(
  sample_barcode,
  patient_barcode,
  is_primary_tumor
)])
mirna_meta[, key_short := norm_id(patient_barcode)]
mirna_keys <- unique(mirna_meta[is_primary_tumor == TRUE, key_short])

# --- RNA sample IDs ---
rna_raw <- readRDS(file.path(input_dir, "raw_rnaseq_counts_matrix.rds"))
rna_keys <- unique(norm_id(colnames(rna_raw)))

# --- train / test splits ---
train_df <- fread(file.path(input_dir, "train_set.csv"))
test_df  <- fread(file.path(input_dir, "test_set_LOCKED.csv"))

train_keys <- norm_id(train_df$patient_barcode)
test_keys  <- norm_id(test_df$patient_barcode)

log_msg(sprintf("   clinical: %d", nrow(clin_raw)))
log_msg(sprintf("   miRNA:    %d", length(mirna_keys)))
log_msg(sprintf("   RNA:      %d", length(rna_keys)))
log_msg(sprintf("   train:    %d", length(train_keys)))
log_msg(sprintf("   test:     %d",  length(test_keys)))

# =============================================================================
# 2) Merge
# =============================================================================
# all patients from the three sources
all_keys <- Reduce(union, list(train_keys, test_keys))

# full matrix
master <- data.table(
  patient_barcode_short = all_keys,
  key_short             = all_keys
)

# --- information from clinical ---
clin_sub <- clin_raw[, .(
  key_short,
  patient_barcode_full = patient_barcode,
  vital_status,
  days_to_death,
  days_to_last_follow_up,
  age_at_diagnosis,
  gender
)]
master <- merge(master, clin_sub, by = "key_short", all.x = TRUE)

# --- presence flags ---
master[, in_clinical := key_short %in% clin_raw$key_short]
master[, in_miRNA    := key_short %in% mirna_keys]
master[, in_RNA      := key_short %in% rna_keys]

# --- split assignment ---
master[, split := fcase(
  key_short %in% train_keys, "train",
  key_short %in% test_keys,  "test",
  default = NA_character_
)]

master[, seed := ifelse(!is.na(split), 42, NA_integer_)]

# --- sample_type ---
master[, sample_type := ifelse(in_miRNA | in_RNA, "Primary Tumor", NA_character_)]

# --- has_survival ---
status_num <- as.integer(tolower(master$vital_status) %in%
                          c("dead", "deceased", "1"))
time_ok <- !is.na(master$days_to_death) | !is.na(master$days_to_last_follow_up)
master[, has_survival := !is.na(vital_status) & time_ok]

# =============================================================================
# 3) Sort and select the final columns
# =============================================================================
setorder(master, split, patient_barcode_short)

master_final <- master[, .(
  patient_barcode_short,
  patient_barcode_full,
  sample_type,
  split,
  seed,
  in_miRNA,
  in_RNA,
  has_survival,
  vital_status,
  age_at_diagnosis,
  gender
)]

# =============================================================================
# 4) Inspection
# =============================================================================
log_msg("")
log_msg("summary:")
log_msg(sprintf("   total rows: %d", nrow(master_final)))
log_msg(sprintf("   train:      %d", sum(master_final$split == "train", na.rm = TRUE)))
log_msg(sprintf("   test:       %d", sum(master_final$split == "test",  na.rm = TRUE)))

if (any(is.na(master_final$split))) {
  log_msg(sprintf("   patients without a split: %d",
                  sum(is.na(master_final$split))))
}

# vital_status distribution
log_msg("   vital_status distribution:")
print(table(master_final$vital_status, useNA = "ifany"))

# check for invalid barcodes
bad_barcodes <- master_final[!grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$",
                                     patient_barcode_short)]
if (nrow(bad_barcodes) > 0) {
  log_msg(sprintf("   %d invalid barcodes", nrow(bad_barcodes)))
  print(bad_barcodes[, 1:3])
}

# =============================================================================
# 5) Save
# =============================================================================
out_csv <- file.path(out_dir, "DataS1_barcodes.csv")
fwrite(master_final, out_csv)

log_msg("")
log_msg("   saved: DataS1_barcodes.csv")
log_msg(sprintf("      size: %.1f KB", file.info(out_csv)$size / 1024))
log_msg(sprintf("      number of rows: %d", nrow(master_final)))
log_msg(sprintf("      number of columns: %d", ncol(master_final)))

# --- save a summary for the DOCX ---
summary_s1 <- data.frame(
  Category = c("Training set",
               "Locked test set",
               "Total (with split)",
               "Patients with miRNA-seq",
               "Patients with RNA-seq",
               "Patients with complete survival"),
  N = c(sum(master_final$split == "train", na.rm = TRUE),
        sum(master_final$split == "test",  na.rm = TRUE),
        sum(!is.na(master_final$split)),
        sum(master_final$in_miRNA, na.rm = TRUE),
        sum(master_final$in_RNA,   na.rm = TRUE),
        sum(master_final$has_survival, na.rm = TRUE)),
  stringsAsFactors = FALSE)

write.csv(summary_s1,
          file.path(out_dir, "DataS1_summary.csv"),
          row.names = FALSE)

log_msg("   saved: DataS1_summary.csv")
log_msg("")
print(summary_s1, row.names = FALSE)
