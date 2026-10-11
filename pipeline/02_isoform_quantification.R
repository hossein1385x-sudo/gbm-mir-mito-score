#!/usr/bin/env Rscript
# =============================================================================
# 02_isoform_level_quant.R
# -----------------------------------------------------------------------------
# Step 2 - isoform-level quantification using the miRBase MIMAT -> name mapping
# (taken from mature.fa.txt, which is in FASTA format)
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# Function to read mature.fa.txt and extract the MIMAT -> NAME mapping
load_mirbase_mapping_from_fasta <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("file ", file_path, " not found!")
  }

  # read all lines of the file
  lines <- readLines(file_path, warn = FALSE)

  # find the header lines (starting with ">")
  header_lines <- grep("^>", lines, value = TRUE)

  if (length(header_lines) == 0) {
    stop("no FASTA header (starting with >) found in the file.")
  }

  # extract the name and MIMAT from each header
  mapping_list <- lapply(header_lines, function(h) {
    # strip the leading ">"
    h_clean <- sub("^>", "", h)
    # split on whitespace (or TAB)
    parts <- strsplit(h_clean, "[ \t]+")[[1]]
    if (length(parts) < 2) return(NULL)
    name <- parts[1]
    # the second part may be "MIMAT0000062" or "MIMAT:0000062"
    mimat_part <- parts[2]
    mimat <- sub("^MIMAT:?", "", mimat_part)  # strip "MIMAT" or "MIMAT:"
    mimat <- paste0("MIMAT", mimat)           # restore with the MIMAT prefix
    if (grepl("^MIMAT[0-9]+$", mimat)) {
      return(data.table(mimat = mimat, mature_name = name))
    } else {
      return(NULL)
    }
  })

  mapping_dt <- rbindlist(mapping_list, use.names = TRUE, fill = TRUE)
  if (nrow(mapping_dt) == 0) {
    stop("no valid mapping could be extracted from the file.")
  }

  # normalisation: remove duplicates by MIMAT
  mapping_dt <- unique(mapping_dt, by = "mimat")
  return(mapping_dt)
}

run_isoform_quant <- function(input_root,
                              exclude_cross_mapped = TRUE,
                              min_samples_expressed = 1) {

  in_path <- file.path(input_root, "raw_mirna.rds")
  if (!file.exists(in_path)) stop("raw_mirna.rds not found in: ", in_path)

  log_msg("==================== Step 02: Isoform-level miRNA Quantification ====================")
  log_msg(sprintf("reading %s ...", in_path))
  mirna_raw <- readRDS(in_path)
  if (!is.data.table(mirna_raw)) mirna_raw <- as.data.table(mirna_raw)
  log_msg(sprintf("  %d raw rows read.", nrow(mirna_raw)))

  # ---- 1) Filter isoforms.quantification ----
  if (!("file_type" %in% names(mirna_raw))) {
    stop("file_type column not found in raw_mirna.rds.")
  }
  iso_dt <- mirna_raw[file_type == "isoforms_quantification"]
  if (nrow(iso_dt) == 0) stop("no isoforms_quantification rows found.")
  log_msg(sprintf("  %d isoform-level rows found.", nrow(iso_dt)))

  # ---- 2) Extract annotation_type and MIMAT from miRNA_region ----
  region_col <- intersect(c("mirna_region", "miRNA_region"), names(iso_dt))
  if (length(region_col) == 0) stop("miRNA_region column not found.")
  region_col <- region_col[1]

  parts <- tstrsplit(iso_dt[[region_col]], ",", fixed = TRUE)
  iso_dt[, annotation_type := trimws(parts[[1]])]
  iso_dt[, mimat := trimws(parts[[2]])]   # the second column is the MIMAT

  # ---- 3) Load the MIMAT->name mapping from the local file ----
  # look for mature.fa.txt in input_root and, if absent, in getwd()
  mapping_file <- file.path(input_root, "mature.fa.txt")
  if (!file.exists(mapping_file)) {
    mapping_file <- file.path(getwd(), "mature.fa.txt")
  }
  log_msg(sprintf("  looking for the mapping file in: %s", mapping_file))

  if (!file.exists(mapping_file)) {
    stop(
      "mature.fa.txt was not found in any of the following paths:\n",
      "  - ", file.path(input_root, "mature.fa.txt"), "\n",
      "  - ", file.path(getwd(), "mature.fa.txt"), "\n",
      "Please download the file from https://www.mirbase.org/ftp/CURRENT/mature.fa.gz,\n",
      "unzip it, rename it to mature.fa.txt, and place it in one of the two folders."
    )
  }

  log_msg("  reading mature.fa.txt ...")
  mapping_dt <- load_mirbase_mapping_from_fasta(mapping_file)
  log_msg(sprintf("  %d MIMAT->name mappings loaded.", nrow(mapping_dt)))

  # ---- 4) Add the mapping to iso_dt ----
  # now both tables have a "mimat" column
  iso_dt <- merge(iso_dt, mapping_dt, by = "mimat", all.x = TRUE)
  iso_dt[, arm := str_extract(mature_name, "(5p|3p)$")]

  n_mapped <- sum(!is.na(iso_dt$mature_name) & iso_dt$annotation_type == "mature")
  log_msg(sprintf("  number of mature MIMATs mapped to a mature name: %d of %d",
                  n_mapped, sum(iso_dt$annotation_type == "mature" & !is.na(iso_dt$mimat))))

  # ---- 5) Filter mature rows with a defined arm ----
  mature_dt <- iso_dt[annotation_type == "mature" & !is.na(arm)]
  if (nrow(mature_dt) == 0) {
    stop("no mature row with a -5p/-3p suffix found. Mapping failed.")
  }

  n_no_arm <- sum(iso_dt$annotation_type == "mature" & is.na(iso_dt$arm), na.rm = TRUE)
  if (n_no_arm > 0) {
    log_msg(sprintf("  warning: %d mature rows without an arm suffix (failed mapping) were dropped.", n_no_arm))
  }

  # ---- 6) Cross-mapped ----
  cm_col <- intersect(c("cross_mapped", "cross-mapped"), names(mature_dt))
  if (length(cm_col) > 0 && exclude_cross_mapped) {
    cm_col <- cm_col[1]
    before <- nrow(mature_dt)
    n_cross <- sum(toupper(as.character(mature_dt[[cm_col]])) == "Y", na.rm = TRUE)
    mature_dt <- mature_dt[toupper(as.character(get(cm_col))) != "Y"]
    log_msg(sprintf("  %d/%d cross-mapped rows removed.", n_cross, before))
  }

  # ---- 7) Aggregation ----
  sample_col <- if ("sample_barcode" %in% names(mature_dt)) "sample_barcode" else "patient_barcode"
  if (!(sample_col %in% names(mature_dt))) {
    stop("sample column (sample_barcode/patient_barcode) not found.")
  }
  mature_dt[, feature_id := mature_name]   # mature name with arm suffix

  read_count_col <- intersect(c("read_count"), names(mature_dt))
  rpm_col <- intersect(c("reads_per_million_mirna_mapped", "reads_per_million_miRNA_mapped"), names(mature_dt))
  if (length(read_count_col) == 0) stop("read_count column not found.")

  agg_cols <- c(read_count_col, if (length(rpm_col) > 0) rpm_col else NULL)
  agg_dt <- mature_dt[, lapply(.SD, function(x) sum(as.numeric(x), na.rm = TRUE)),
                      by = c(sample_col, "feature_id"), .SDcols = agg_cols]

  n_samples <- length(unique(agg_dt[[sample_col]]))
  n_features <- length(unique(agg_dt$feature_id))
  n_5p <- length(unique(agg_dt$feature_id[grepl("-5p$", agg_dt$feature_id)]))
  n_3p <- length(unique(agg_dt$feature_id[grepl("-3p$", agg_dt$feature_id)]))
  log_msg(sprintf("  after aggregation: %d samples × %d features (5p=%d, 3p=%d).",
                  n_samples, n_features, n_5p, n_3p))

  # ---- 8) Minimum-expression filter ----
  expressed_counts <- agg_dt[, .(n_expressed = sum(get(read_count_col[1]) > 0)), by = feature_id]
  keep_features <- expressed_counts[n_expressed >= min_samples_expressed, feature_id]
  n_dropped <- n_features - length(keep_features)
  if (n_dropped > 0) {
    log_msg(sprintf("  %d features expressed in fewer than %d samples were removed.", n_dropped, min_samples_expressed))
  }
  agg_dt <- agg_dt[feature_id %in% keep_features]

  # ---- 9) Matrices ----
  counts_wide <- dcast(agg_dt, as.formula(paste(sample_col, "~ feature_id")),
                       value.var = read_count_col[1], fill = 0)
  rpm_wide <- NULL
  if (length(rpm_col) > 0) {
    rpm_wide <- dcast(agg_dt, as.formula(paste(sample_col, "~ feature_id")),
                      value.var = rpm_col[1], fill = 0)
  }

  # ---- 10) Reference mapping ----
  precursor_col <- intersect(c("mirna_id", "miRNA_ID"), names(mature_dt))
  mapping_ref <- unique(mature_dt[, .(
    feature_id = mature_name,
    arm,
    precursor_id = if (length(precursor_col) > 0) get(precursor_col[1]) else NA_character_
  )])
  mapping_ref <- mapping_ref[feature_id %in% keep_features]
  setorder(mapping_ref, feature_id)

  # ---- 11) Save ----
  out_counts  <- file.path(input_root, "mirna_isoform_counts_matrix.csv")
  out_rpm     <- file.path(input_root, "mirna_isoform_rpm_matrix.csv")
  out_mapping <- file.path(input_root, "isoform_mapping_reference.csv")
  out_log     <- file.path(input_root, "isoform_quant_log.txt")

  fwrite(counts_wide, out_counts)
  log_msg(sprintf("saved: %s", out_counts))
  if (!is.null(rpm_wide)) fwrite(rpm_wide, out_rpm)
  fwrite(mapping_ref, out_mapping)

  log_lines <- c(
    "==================== Isoform Quantification Log ====================",
    sprintf("input: %s", in_path),
    sprintf("number of final features: %d", length(keep_features)),
    sprintf("5p: %d, 3p: %d", n_5p, n_3p),
    sprintf("outputs: %s, %s, %s", out_counts, out_rpm, out_mapping)
  )
  writeLines(log_lines, out_log)
  log_msg("Step 02 finished.")

  invisible(list(counts = counts_wide, rpm = rpm_wide, mapping = mapping_ref))
}

## ---- Entry point ----
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  parse_args <- function(argv) {
    opts <- list(); i <- 1
    while (i <= length(argv)) {
      if (substr(argv[i], 1, 2) == "--") {
        key <- sub("^--", "", argv[i])
        if (i < length(argv) && substr(argv[i + 1], 1, 2) != "--") {
          opts[[key]] <- argv[i + 1]; i <- i + 2
        } else {
          opts[[key]] <- TRUE; i <- i + 1
        }
      } else i <- i + 1
    }
    opts
  }
  opts <- parse_args(args)
  default_root <- WORK_DIR
  run_isoform_quant(
    input_root            = opts[["input-root"]] %||% default_root,
    exclude_cross_mapped  = !isTRUE(opts[["keep-cross-mapped"]]),
    min_samples_expressed = as.integer(opts[["min-samples-expressed"]] %||% 1)
  )
}
