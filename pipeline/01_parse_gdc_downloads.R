#!/usr/bin/env Rscript
# =============================================================================
# 01b_parse_downloaded_data.R
# -----------------------------------------------------------------------------
# Part 2 of Step 1 of the miR-Mito Score pipeline for TCGA-GBM
#
# Prerequisite: the raw RNA-seq and miRNA files were already downloaded with gdc-client, using
#           manifests exported from the GDC Portal, into this layout:
#             <output_root>/01_output/rnaseq/<uuid>/<filename>
#             <output_root>/01_output/mirna/<uuid>/<filename>
#
# This script performs four tasks:
#   1) Scan the rnaseq/ and mirna/ folders and collect all downloaded UUIDs.
#   2) Query the GDC /files endpoint (batched) for these UUIDs to get each file's real sample/patient barcode and
#      project_id, because the manifest exported by hand from the GDC Portal does not contain
#      this information. It also asserts that every file belongs to TCGA-GBM:
#      if any UUID belongs to another project,
#      the script stops with an explicit error instead of continuing silently.
#   3) Parse each raw file (no filtering or normalisation) and attach the barcode and
#      tumour/normal label from parse_tcga_barcode.
#   4) Query the GDC /cases endpoint for the clinical data of the whole project (which had not been
#      downloaded up to this point).
#
# Outputs (in WORK_DIR/):
#   raw_rnaseq.rds, raw_mirna.rds, raw_clinical.rds,
#   file_uuid_metadata.csv (uuid -> barcode/project mapping, for auditing),
#   parse_log.txt
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(httr)
  library(stringr)
})

## ---- Global settings -------------------------------------------------------
EXPECTED_PROJECT_ID <- "TCGA-GBM"
GDC_API_BASE         <- "https://api.gdc.cancer.gov"

TCGA_BARCODE_REGEX <- "^(TCGA)-([A-Z0-9]{2})-([A-Z0-9]{4})-([0-9]{2})([A-Z])-([0-9]{2})([A-Z])-([0-9]{4})-([0-9]{2})$"

SAMPLE_TYPE_MAP <- c(
  "01" = "Primary Solid Tumor", "02" = "Recurrent Solid Tumor",
  "03" = "Primary Blood Derived Cancer - Peripheral Blood",
  "04" = "Recurrent Blood Derived Cancer - Bone Marrow",
  "05" = "Additional - New Primary", "06" = "Metastatic",
  "07" = "Additional Metastatic", "08" = "Human Tumor Original Cells",
  "09" = "Primary Blood Derived Cancer - Bone Marrow",
  "10" = "Blood Derived Normal", "11" = "Solid Tissue Normal",
  "12" = "Buccal Cell Normal", "13" = "EBV Immortalized Normal",
  "14" = "Bone Marrow Normal", "15" = "control analyte",
  "20" = "Control Analyte", "40" = "Recurrent Blood Derived Cancer - Peripheral Blood",
  "50" = "Cell Line Derived Xenograft Tissue", "60" = "Primary Xenograft Tissue",
  "61" = "Cell Line Derived Xenograft Tissue", "99" = "Granulocytes"
)

## ---- Basic helper functions ------------------------------------------------
log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

safe_get <- function(x, ...) {
  keys <- list(...)
  for (k in keys) { if (is.null(x)) return(NULL); x <- x[[k]] }
  x
}

with_retry <- function(fn, max_retries = 5, base_delay = 2, what = "operation") {
  attempt <- 1
  repeat {
    result <- tryCatch(list(ok = TRUE, value = fn()),
                       error = function(e) list(ok = FALSE, err = conditionMessage(e)))
    if (result$ok) return(result$value)
    if (attempt >= max_retries) stop(sprintf("'%s' failed after %d attempts: %s", what, max_retries, result$err))
    delay <- base_delay * (2 ^ (attempt - 1)) + runif(1, 0, 1)
    log_msg(sprintf("warning: '%s' failed (attempt %d/%d): %s — retrying in %.1f seconds.",
                    what, attempt, max_retries, result$err, delay))
    Sys.sleep(delay)
    attempt <- attempt + 1
  }
}

## ---- Function 1: scan the folder and collect UUIDs / file paths ------------
# expected structure: <root>/<uuid-folder>/<filename>
scan_uuid_folders <- function(root_dir) {
  if (!dir.exists(root_dir)) {
    warning("folder not found: ", root_dir)
    return(data.table(file_uuid = character(0), file_path = character(0), file_name = character(0)))
  }
  uuid_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  # only folders whose names really look like a UUID (fail-safe against extra folders such as temp)
  uuid_pattern <- "^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
  uuid_dirs <- uuid_dirs[grepl(uuid_pattern, basename(uuid_dirs))]

  rows <- lapply(uuid_dirs, function(d) {
    uuid <- basename(d)
    files <- list.files(d, full.names = TRUE, recursive = FALSE)
    # ignore incomplete (.partial) files and gdc-client logs
    files <- files[!grepl("\\.partial$|\\.logfile$|annotations\\.txt$", files)]
    if (length(files) == 0) {
      log_msg(sprintf("warning: folder %s is empty or contains only an incomplete file — ignored.", uuid))
      return(NULL)
    }
    if (length(files) > 1) {
      log_msg(sprintf("warning: folder %s contains more than one file; all are listed: %s",
                      uuid, paste(basename(files), collapse = ", ")))
    }
    data.table(file_uuid = uuid, file_path = files, file_name = basename(files))
  })
  rbindlist(rows[!vapply(rows, is.null, logical(1))], fill = TRUE)
}

## ---- Function 2: query GDC /files for a given list of UUIDs (batched) -----
gdc_query_files_by_uuid <- function(file_uuids, api_base = GDC_API_BASE, batch_size = 300) {
  file_uuids <- unique(file_uuids)
  fields <- paste(c(
    "file_id", "file_name", "data_type",
    "cases.case_id", "cases.submitter_id",
    "cases.samples.sample_id", "cases.samples.submitter_id",
    "cases.samples.sample_type_id", "cases.project.project_id"
  ), collapse = ",")

  all_rows <- list()
  for (chunk_start in seq(1, length(file_uuids), by = batch_size)) {
    chunk_end <- min(chunk_start + batch_size - 1, length(file_uuids))
    chunk <- file_uuids[chunk_start:chunk_end]

    filters <- list(op = "in", content = list(field = "file_id", value = as.list(chunk)))
    body <- list(filters = filters, fields = fields, format = "JSON", size = length(chunk))

    resp <- with_retry(function() {
      r <- POST(paste0(api_base, "/files"), body = body, encode = "json",
               timeout(60), add_headers(.headers = c("Content-Type" = "application/json")))
      if (http_error(r)) stop(sprintf("HTTP %d from /files", status_code(r)))
      r
    }, what = sprintf("GDC /files metadata query (chunk %d-%d)", chunk_start, chunk_end))

    parsed <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyDataFrame = FALSE)
    hits <- parsed$data$hits
    if (length(hits) == 0) next
    all_rows <- c(all_rows, hits)
    log_msg(sprintf("  ... metadata received for %d/%d files", chunk_end, length(file_uuids)))
  }

  if (length(all_rows) == 0) {
    stop("no metadata was returned from GDC /files for the given UUIDs — check that the UUIDs are correct.")
  }

  rbindlist(lapply(all_rows, function(h) {
    cases <- safe_get(h, "cases")
    case1 <- if (!is.null(cases) && length(cases) > 0) cases[[1]] else NULL
    proj <- safe_get(case1, "project")
    proj_id <- safe_get(proj, "project_id") %||% NA_character_

    samples <- safe_get(case1, "samples")
    sample1 <- if (!is.null(samples) && length(samples) > 0) samples[[1]] else NULL

    data.table(
      file_id            = safe_get(h, "file_id") %||% NA_character_,
      file_name_api       = safe_get(h, "file_name") %||% NA_character_,
      data_type           = safe_get(h, "data_type") %||% NA_character_,
      patient_barcode      = safe_get(case1, "submitter_id") %||% NA_character_,
      sample_id           = safe_get(sample1, "sample_id") %||% NA_character_,
      sample_barcode       = safe_get(sample1, "submitter_id") %||% NA_character_,
      sample_type_id       = safe_get(sample1, "sample_type_id") %||% NA_character_,
      project_id_seen      = as.character(proj_id)[1]
    )
  }), fill = TRUE)
}

## ---- Function 3: defensive assertion on the returned metadata --------------
assert_all_from_expected_project <- function(meta_dt, expected = EXPECTED_PROJECT_ID) {
  seen <- unique(na.omit(meta_dt$project_id_seen))
  unexpected <- setdiff(seen, expected)
  if (length(unexpected) > 0) {
    bad_uuids <- meta_dt[project_id_seen %in% unexpected, file_id]
    stop("project assertion violated: ", length(bad_uuids), " files belong to a project other than ", expected,
         " — projects found: ", paste(unexpected, collapse = ", "),
         "\nexample problematic UUIDs: ", paste(head(bad_uuids, 5), collapse = ", "),
         "\nthis is exactly the LUAD leakage pattern that must be prevented.")
  }
  n_missing <- sum(is.na(meta_dt$project_id_seen))
  if (n_missing > 0) log_msg(sprintf("warning: project_id could not be extracted for %d/%d files.", n_missing, nrow(meta_dt)))
  invisible(TRUE)
}

## ---- Function 4: parse a TCGA barcode --------------------------------------
parse_tcga_barcode <- function(barcode) {
  if (length(barcode) != 1L || is.na(barcode)) return(NA)
  m <- str_match(barcode, TCGA_BARCODE_REGEX)
  if (is.na(m[1, 1])) {
    m2 <- str_match(barcode, "^(TCGA)-([A-Z0-9]{2})-([A-Z0-9]{4})-([0-9]{2})([A-Z])")
    if (is.na(m2[1, 1])) return(NA)
    sample_code <- m2[1, 5]
    return(data.frame(
      sample_code = sample_code, sample_type = SAMPLE_TYPE_MAP[sample_code] %||% NA_character_,
      is_tumor = sample_code %in% as.character(c(1:9, 40, 50, 60, 61)),
      is_normal = sample_code %in% as.character(c(10:15, 20)),
      is_primary_tumor = sample_code == "01", is_recurrent_tumor = sample_code == "02",
      is_metastatic = sample_code == "06", stringsAsFactors = FALSE
    ))
  }
  sample_code <- m[1, 5]
  data.frame(
    sample_code = sample_code, sample_type = SAMPLE_TYPE_MAP[sample_code] %||% NA_character_,
    is_tumor = sample_code %in% as.character(c(1:9, 40, 50, 60, 61)),
    is_normal = sample_code %in% as.character(c(10:15, 20)),
    is_primary_tumor = sample_code == "01", is_recurrent_tumor = sample_code == "02",
    is_metastatic = sample_code == "06", stringsAsFactors = FALSE
  )
}

## ---- Function 5: read an RNA-seq file (STAR Counts / HTSeq legacy) ---------
read_single_rnaseq <- function(file_path) {
  stopifnot(file.exists(file_path))
  bn <- basename(file_path)
  if (grepl("augmented_star_gene_counts|star_gene_counts", bn, ignore.case = TRUE)) {
    first_lines <- readLines(file_path, n = 6, warn = FALSE)
    comment_lines <- grep("^#", first_lines)
    skip_n <- if (length(comment_lines) > 0) max(comment_lines) else 0
    dt <- fread(file_path, skip = skip_n, header = TRUE, data.table = TRUE, check.names = FALSE)
    dt[, format := "STAR_Counts"]
    return(dt)
  }
  head_line <- readLines(file_path, n = 1, warn = FALSE)
  if (grepl("^gene_id\\s", head_line) || grepl("^#", head_line)) {
    comment_lines <- grep("^#", readLines(file_path, n = 8, warn = FALSE))
    skip_n <- if (length(comment_lines) > 0) max(comment_lines) else 0
    dt <- fread(file_path, skip = skip_n, header = TRUE, data.table = TRUE)
    dt[, format := "STAR_Counts"]
    return(dt)
  }
  dt <- fread(file_path, header = FALSE, col.names = c("gene_id", "unstranded"), data.table = TRUE)
  dt[, format := "HTSeq_Counts_legacy_unknown"]
  dt
}

## ---- Function 6: read a miRNA-seq quantification file ----------------------
read_single_mirna <- function(file_path) {
  stopifnot(file.exists(file_path))
  bn <- basename(file_path)
  dt <- fread(file_path, header = TRUE, data.table = TRUE)
  setnames(dt, tolower(names(dt)))
  if ("cross-mapped" %in% names(dt)) setnames(dt, "cross-mapped", "cross_mapped")
  dt[, file_type := if (grepl("isoforms\\.quantification", bn)) "isoforms_quantification"
                     else if (grepl("mirnas\\.quantification", bn)) "mirnas_quantification"
                     else "mirna_unknown"]
  dt
}

## ---- Function 7: query clinical data from GDC /cases (paginated) -----------
gdc_query_clinical <- function(project_id, api_base = GDC_API_BASE, page_size = 500) {
  fields <- paste(c(
    "case_id", "submitter_id",
    "demographic.vital_status", "demographic.days_to_death", "demographic.gender",
    "demographic.race", "demographic.ethnicity",
    "diagnoses.days_to_last_follow_up", "diagnoses.age_at_diagnosis",
    "diagnoses.primary_diagnosis", "diagnoses.tumor_grade",
    "diagnoses.morphology", "diagnoses.site_of_resection_or_biopsy"
  ), collapse = ",")
  filters <- list(op = "in", content = list(field = "project.project_id", value = list(project_id)))

  all_rows <- list(); from <- 0
  repeat {
    body <- list(filters = filters, fields = fields, format = "JSON",
                 size = page_size, from = from, expand = "demographic,diagnoses")
    resp <- with_retry(function() {
      r <- POST(paste0(api_base, "/cases"), body = body, encode = "json",
               timeout(60), add_headers(.headers = c("Content-Type" = "application/json")))
      if (http_error(r)) stop(sprintf("HTTP %d from /cases", status_code(r)))
      r
    }, what = sprintf("GDC /cases query (from=%d)", from))

    parsed <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyDataFrame = FALSE)
    hits <- parsed$data$hits
    total <- safe_get(parsed, "data", "pagination", "total") %||% 0
    if (length(hits) == 0) break
    all_rows <- c(all_rows, hits)
    from <- from + length(hits)
    if (from >= total) break
  }
  if (length(all_rows) == 0) { warning("no clinical records found."); return(data.table()) }

  rbindlist(lapply(all_rows, function(h) {
    demo <- safe_get(h, "demographic")
    diag_list <- safe_get(h, "diagnoses")
    diag1 <- if (!is.null(diag_list) && length(diag_list) > 0) diag_list[[1]] else NULL
    data.table(
      case_id = safe_get(h, "case_id") %||% NA_character_,
      patient_barcode = safe_get(h, "submitter_id") %||% NA_character_,
      vital_status = safe_get(demo, "vital_status") %||% NA_character_,
      days_to_death = safe_get(demo, "days_to_death") %||% NA_integer_,
      gender = safe_get(demo, "gender") %||% NA_character_,
      race = safe_get(demo, "race") %||% NA_character_,
      ethnicity = safe_get(demo, "ethnicity") %||% NA_character_,
      days_to_last_follow_up = safe_get(diag1, "days_to_last_follow_up") %||% NA_integer_,
      age_at_diagnosis = safe_get(diag1, "age_at_diagnosis") %||% NA_integer_,
      primary_diagnosis = safe_get(diag1, "primary_diagnosis") %||% NA_character_,
      tumor_grade = safe_get(diag1, "tumor_grade") %||% NA_character_,
      morphology = safe_get(diag1, "morphology") %||% NA_character_,
      site_of_resection_or_biopsy = safe_get(diag1, "site_of_resection_or_biopsy") %||% NA_character_
    )
  }), fill = TRUE)
}

## ---- Main function: run_parse_downloaded -----------------------------------
run_parse_downloaded <- function(output_root,
                                 project_id = EXPECTED_PROJECT_ID,
                                 rnaseq_subdir = "01_output/rnaseq",
                                 mirna_subdir  = "01_output/mirna") {

  output_dir <- file.path(output_root, "01_output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  log_msg("==================== Step 01 (part 2): Parse the downloaded files ====================")

  # ---- 1) Scan the folders -----------------------------------------------------
  rnaseq_root <- file.path(output_root, rnaseq_subdir)
  mirna_root  <- file.path(output_root, mirna_subdir)
  log_msg(sprintf("scanning %s ...", rnaseq_root))
  rnaseq_scan <- scan_uuid_folders(rnaseq_root)
  log_msg(sprintf("  %d RNA-seq files found.", nrow(rnaseq_scan)))
  log_msg(sprintf("scanning %s ...", mirna_root))
  mirna_scan <- scan_uuid_folders(mirna_root)
  log_msg(sprintf("  %d miRNA files found.", nrow(mirna_scan)))

  if (nrow(rnaseq_scan) == 0 && nrow(mirna_scan) == 0) {
    stop("no files found in ", rnaseq_root, " or ", mirna_root, ". Check the paths.")
  }

  # ---- 2) Query metadata (barcode + project assertion) for all UUIDs together
  all_uuids <- unique(c(rnaseq_scan$file_uuid, mirna_scan$file_uuid))
  log_msg(sprintf("querying metadata of %d UUIDs from GDC /files ...", length(all_uuids)))
  meta_dt <- gdc_query_files_by_uuid(all_uuids)
  assert_all_from_expected_project(meta_dt, expected = project_id)
  log_msg(sprintf("metadata received for %d files and project_id confirmed.", nrow(meta_dt)))

  # save for audit
  fwrite(meta_dt, file.path(output_dir, "file_uuid_metadata.csv"))

  # ---- 3) Join the metadata to each file and parse -----------------------------
  enrich_and_read <- function(scan_dt, reader_fn) {
    if (nrow(scan_dt) == 0) return(data.table())
    merged <- merge(scan_dt, meta_dt, by.x = "file_uuid", by.y = "file_id", all.x = TRUE)
    n_no_meta <- sum(is.na(merged$patient_barcode))
    if (n_no_meta > 0) {
      log_msg(sprintf("warning: %d/%d files have no barcode metadata (not returned by GDC).", n_no_meta, nrow(merged)))
    }
    out_list <- vector("list", nrow(merged))
    for (i in seq_len(nrow(merged))) {
      dt <- tryCatch(reader_fn(merged$file_path[i]), error = function(e) {
        log_msg(sprintf("error reading %s: %s", merged$file_path[i], conditionMessage(e)))
        NULL
      })
      if (is.null(dt)) next
      dt[, `:=`(file_uuid = merged$file_uuid[i], source_file = merged$file_name[i],
                patient_barcode = merged$patient_barcode[i], sample_barcode = merged$sample_barcode[i])]
      bc_info <- parse_tcga_barcode(merged$sample_barcode[i])
      if (is.data.frame(bc_info)) {
        dt[, `:=`(sample_type = bc_info$sample_type, is_tumor = bc_info$is_tumor,
                  is_normal = bc_info$is_normal, is_primary_tumor = bc_info$is_primary_tumor,
                  is_recurrent_tumor = bc_info$is_recurrent_tumor, is_metastatic = bc_info$is_metastatic)]
      }
      out_list[[i]] <- dt
      if (i %% 50 == 0) log_msg(sprintf("  ... %d/%d parsed", i, nrow(merged)))
    }
    rbindlist(out_list, fill = TRUE, use.names = TRUE)
  }

  log_msg("parsing RNA-seq files ...")
  rnaseq_combined <- enrich_and_read(rnaseq_scan, read_single_rnaseq)
  log_msg(sprintf("RNA-seq parsed: %d rows.", nrow(rnaseq_combined)))

  log_msg("parsing miRNA files ...")
  mirna_combined <- enrich_and_read(mirna_scan, read_single_mirna)
  log_msg(sprintf("miRNA parsed: %d rows.", nrow(mirna_combined)))

  # ---- 4) Clinical data ----------------------------------------------------------
  log_msg("querying clinical data from GDC /cases ...")
  clinical_dt <- tryCatch(gdc_query_clinical(project_id), error = function(e) {
    warning("clinical query failed: ", conditionMessage(e)); data.table()
  })
  log_msg(sprintf("clinical: %d patients.", nrow(clinical_dt)))

  # ---- 5) Save -------------------------------------------------------------------
  out_rnaseq   <- file.path(output_dir, "raw_rnaseq.rds")
  out_mirna    <- file.path(output_dir, "raw_mirna.rds")
  out_clinical <- file.path(output_dir, "raw_clinical.rds")
  out_log      <- file.path(output_dir, "parse_log.txt")

  saveRDS(rnaseq_combined, out_rnaseq, compress = "xz")
  saveRDS(mirna_combined,  out_mirna,  compress = "xz")
  saveRDS(clinical_dt,     out_clinical, compress = "xz")

  n_patients_rnaseq <- length(unique(na.omit(rnaseq_combined$patient_barcode)))
  n_patients_mirna  <- length(unique(na.omit(mirna_combined$patient_barcode)))

  log_lines <- c(
    "==================== Parse Downloaded Data Log (Step 01, part 2) ====================",
    sprintf("run_started_at   : %s", as.character(Sys.time())),
    sprintf("project_id       : %s", project_id),
    sprintf("output_root      : %s", output_root),
    "---- RNA-seq ----",
    sprintf("files_parsed     : %d", nrow(rnaseq_scan)),
    sprintf("rows_total       : %d", nrow(rnaseq_combined)),
    sprintf("unique_patients  : %d", n_patients_rnaseq),
    "---- miRNA ----",
    sprintf("files_parsed     : %d", nrow(mirna_scan)),
    sprintf("rows_total       : %d", nrow(mirna_combined)),
    sprintf("unique_patients  : %d", n_patients_mirna),
    "---- Clinical ----",
    sprintf("cases_found      : %d", nrow(clinical_dt)),
    "---- Outputs ----",
    sprintf("raw_rnaseq.rds   : %s", out_rnaseq),
    sprintf("raw_mirna.rds    : %s", out_mirna),
    sprintf("raw_clinical.rds : %s", out_clinical),
    sprintf("file_uuid_metadata.csv : %s", file.path(output_dir, "file_uuid_metadata.csv")),
    "========================================================================================"
  )
  writeLines(log_lines, out_log)
  log_msg(sprintf("log saved: %s", out_log))
  log_msg("continue to Step 02 (isoform-level quantification).")

  invisible(list(rnaseq = rnaseq_combined, mirna = mirna_combined,
                clinical = clinical_dt, log_path = out_log))
}

## ---- Command-line entry point ----------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  parse_args <- function(argv) {
    opts <- list(); i <- 1
    while (i <= length(argv)) {
      if (substr(argv[i], 1, 2) == "--") {
        key <- sub("^--", "", argv[i])
        if (i < length(argv) && substr(argv[i + 1], 1, 2) != "--") { opts[[key]] <- argv[i + 1]; i <- i + 2 }
        else { opts[[key]] <- TRUE; i <- i + 1 }
      } else i <- i + 1
    }
    opts
  }
  opts <- parse_args(args)
  run_parse_downloaded(
    output_root = opts[["output-root"]] %||% getwd(),
    project_id  = opts[["project-id"]] %||% EXPECTED_PROJECT_ID
  )
}
