#!/usr/bin/env Rscript
# combine_benchmarks.R
# Read the outputs of the 5 benchmarks and merge them into shared CSV files
#
# Inputs (each in a sub-folder of WORK_DIR):
#   benchmark_out              -> Srinivasan2011 (first run)
#   benchmark_out_hayes        → Hayes2015
#   benchmark_out_li2014       -> Li2014 (5 subgroups)
#   benchmark_out_bafiti2022   → Bafiti2022
#   benchmark_out_niyazi2016   → Niyazi2016
#
# Outputs:
#   COMBINED_performance.csv        -> C-index of all models in all cohorts
#   COMBINED_performance_wide.csv   -> wide table (cohort x model) for the paper
#   COMBINED_coverage.csv           -> coverage of each signature
#   COMBINED_delta.csv              -> OWN minus each competitor
#   COMBINED_OWN_summary.csv        -> OWN C-index in the three cohorts

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})

source(here::here("config.R"))   # defines WORK_DIR (see README)
OUT <- WORK_DIR
cat("Working in:", OUT, "\n\n")

# Order matters: the first occurrence wins in the de-duplication
# (newest runs first, the old Srinivasan run last)
benches <- c(
  "benchmark_out_niyazi2016",
  "benchmark_out_bafiti2022",
  "benchmark_out_li2014",
  "benchmark_out_hayes",
  "benchmark_out"
)

cat("=== Folders ===\n")
for (b in benches) {
  p <- file.path(OUT, b)
  cat(sprintf("  %-30s %s\n", b, if (dir.exists(p)) "OK" else "MISSING"))
}
cat("\n")

# ---------------------------------------------------------------- helpers
read_one <- function(dir_name, file_name) {
  path <- file.path(OUT, dir_name, file_name)
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(read_csv(path, show_col_types = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$source_run <- dir_name
  df
}

combine <- function(file_name, key_cols, out_name) {
  df <- bind_rows(lapply(benches, read_one, file_name = file_name))
  if (nrow(df) == 0) {
    cat("  (no data):", file_name, "\n")
    return(invisible(NULL))
  }
  # Drop duplicates; the first occurrence (by the order above) is kept
  df <- df[!duplicated(df[, key_cols, drop = FALSE]), ]
  write_csv(df, file.path(OUT, out_name))
  cat(sprintf("  %-34s %d rows\n", out_name, nrow(df)))
  invisible(df)
}

# ---------------------------------------------------------------- 1-3: Collect
cat("=== Combining ===\n")
perf  <- combine("benchmark_performance.csv", c("cohort","model"),
                 "COMBINED_performance.csv")
cov   <- combine("benchmark_coverage.csv",    c("cohort","signature"),
                 "COMBINED_coverage.csv")
delta <- combine("benchmark_delta.csv",       c("cohort","comparison"),
                 "COMBINED_delta.csv")

# ---------------------------------------------------------------- 4: OWN summary
own_raw <- bind_rows(lapply(benches, function(b) {
  bind_rows(read_one(b, "OWN_summary.csv"),
            read_one(b, "OWN_out_of_sample_summary.csv"))
}))
if (nrow(own_raw) > 0) {
  own <- own_raw[!duplicated(own_raw$cohort), ]
  write_csv(own, file.path(OUT, "COMBINED_OWN_summary.csv"))
  cat(sprintf("  %-34s %d rows\n", "COMBINED_OWN_summary.csv", nrow(own)))
}

# ---------------------------------------------------------------- 5: wide table
if (!is.null(perf) && nrow(perf) > 0) {
  wide <- perf %>%
    select(cohort, model, n, events, C, C_lo, C_hi) %>%
    distinct(cohort, model, .keep_all = TRUE) %>%
    pivot_wider(
      id_cols     = c(cohort, n, events),
      names_from  = model,
      values_from = c(C, C_lo, C_hi),
      names_sep   = "_"
    )
  write_csv(wide, file.path(OUT, "COMBINED_performance_wide.csv"))
  cat(sprintf("  %-34s %d rows x %d cols\n",
              "COMBINED_performance_wide.csv", nrow(wide), ncol(wide)))
}

# ---------------------------------------------------------------- 6: pretty print
if (!is.null(perf) && nrow(perf) > 0) {
  cat("\n=========================================================\n")
  cat("C-index by cohort (sorted descending)\n")
  cat("=========================================================\n")
  for (co in unique(perf$cohort)) {
    cat("\n###", co, "###\n")
    d <- perf %>%
      filter(cohort == co) %>%
      arrange(desc(C)) %>%
      select(model, n, events, C, C_lo, C_hi)
    print(as.data.frame(d), row.names = FALSE)
  }

  if (!is.null(delta) && nrow(delta) > 0) {
    cat("\n=========================================================\n")
    cat("Paired differences: OWN - competitor\n")
    cat("=========================================================\n")
    d <- delta %>%
      select(cohort, comparison, delta_C, lo, hi, p_boot_2sided) %>%
      arrange(cohort, desc(delta_C))
    print(as.data.frame(d), row.names = FALSE)
  }
}

cat("\n=========================================================\n")
cat("Done. Combined files written to:\n  ", OUT, "\n")
cat("=========================================================\n")
cat("  COMBINED_performance.csv\n")
cat("  COMBINED_performance_wide.csv\n")
cat("  COMBINED_coverage.csv\n")
cat("  COMBINED_delta.csv\n")
cat("  COMBINED_OWN_summary.csv\n")
