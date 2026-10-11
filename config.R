# =============================================================================
# config.R - the ONLY place where paths are defined.
# Every script starts with:  source(here::here("config.R"))
#
# WORK_DIR is the folder that holds ALL data and results (raw GDC downloads, intermediate
# files, signature_v2/, supplementary_final/, figures_final/, cgga_raw/, reference/ ...).
# It is NOT part of the repository.
#
# Set it once, outside the code, in ~/.Renviron (or the project-level .Renviron):
#     GBM_WORK_DIR="D:/gbm_work"
# If the variable is not set, <repo>/work is used.
# =============================================================================
if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required: install.packages('here')")
}

PROJECT_DIR <- here::here()
WORK_DIR    <- Sys.getenv("GBM_WORK_DIR", unset = file.path(PROJECT_DIR, "work"))
WORK_DIR    <- normalizePath(WORK_DIR, winslash = "/", mustWork = FALSE)
dir.create(WORK_DIR, showWarnings = FALSE, recursive = TRUE)

if (!nzchar(Sys.getenv("GBM_WORK_DIR"))) {
  message("[config] GBM_WORK_DIR is not set; using ", WORK_DIR)
}
