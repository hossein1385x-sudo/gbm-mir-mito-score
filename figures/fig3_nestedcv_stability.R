#!/usr/bin/env Rscript
# =============================================================================
# STEP 18c (Signature v2) — Figure 3 (A + B Combined)
# -----------------------------------------------------------------------------
# Panel A: Distribution of outer-fold C-index (nested CV)
# Panel B: miRNA selection frequency (top 20)
#
# Output: Fig3_NestedCV_and_Stability_600dpi.png/.pdf
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(forcats)
  library(patchwork)
  library(scales)
})

# =============================================================================
# 1) Paths
# =============================================================================
INPUT_DIR  <- paste0(WORK_DIR, "/")
SIG_DIR    <- file.path(INPUT_DIR, "signature_v2")
OUTPUT_DIR <- file.path(INPUT_DIR, "figures_final")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

out_path <- function(f) file.path(OUTPUT_DIR, f)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

# =============================================================================
# 2) Shared theme
# =============================================================================
theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = "sans") %+replace%
    theme(
      panel.grid.major  = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.ticks        = element_line(color = "black", linewidth = 0.4),
      axis.title        = element_text(face = "bold", size = 10),
      axis.text         = element_text(color = "black", size = 9),
      plot.title        = element_text(face = "bold", size = 11, hjust = 0,
                                       margin = margin(b = 4)),
      plot.subtitle     = element_text(size = 9, hjust = 0, color = "grey30",
                                       margin = margin(b = 6)),
      legend.position   = "top",
      legend.title      = element_blank(),
      legend.text       = element_text(size = 8.5),
      legend.key        = element_blank(),
      plot.margin       = margin(t = 8, r = 8, b = 8, l = 8, unit = "pt")
    )
}

log_msg("============================================================")
log_msg("Figure 3 (Signature v2) — Panel A + B")
log_msg("============================================================")

# =============================================================================
# PANEL A — Distribution of the outer-fold C-index
# =============================================================================
log_msg("")
log_msg("Panel A — Outer-fold C-index distribution")

# possible paths for the per-fold C-index file
f_a_candidates <- c(
  file.path(SIG_DIR, "signature_v2_cv_performance.csv"),
  file.path(SIG_DIR, "nested_cv_performance_repeated.csv"),
  file.path(INPUT_DIR, "signature_v2_cv_performance.csv"),
  file.path(INPUT_DIR, "nested_cv_performance_repeated.csv")
)
f_a_path <- f_a_candidates[file.exists(f_a_candidates)][1]
if (is.na(f_a_path)) stop("nested CV performance file not found.")

log_msg(sprintf("   reading: %s", basename(f_a_path)))
df_a <- read_csv(f_a_path, show_col_types = FALSE)

# C-index column
ci_col <- intersect(c("c_index", "C_index", "cindex", "C"),
                    names(df_a))[1]
if (is.na(ci_col)) stop("C-index column not found.")
ci_vals <- as.numeric(df_a[[ci_col]])
ci_vals <- ci_vals[is.finite(ci_vals)]

mean_ci <- mean(ci_vals)
sd_ci   <- sd(ci_vals)
n_ci    <- length(ci_vals)

log_msg(sprintf("   n = %d | mean = %.3f | SD = %.3f", n_ci, mean_ci, sd_ci))

df_a_plot <- data.frame(C_index = ci_vals)

pA <- ggplot(df_a_plot, aes(x = C_index)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 15,
                 fill = "grey75", color = "white", linewidth = 0.3) +
  geom_density(linewidth = 0.9, color = "#2166AC") +
  geom_vline(xintercept = mean_ci, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  annotate("text", x = mean_ci, y = Inf, label = sprintf("mean = %.3f", mean_ci),
           vjust = 1.5, hjust = -0.05, size = 3, fontface = "bold") +
  labs(
    title    = "Nested Cross-Validation Discrimination",
    subtitle = sprintf("Outer-fold C-index (n = %d valid folds)", n_ci),
    x = "Outer-fold C-index",
    y = "Density"
  ) +
  theme_pub()

# =============================================================================
# PANEL B — miRNA selection frequency
# =============================================================================
log_msg("")
log_msg("Panel B — miRNA selection frequency")

f_b_path <- file.path(SIG_DIR, "signature_v2_selection_frequency.csv")
if (!file.exists(f_b_path)) {
  f_b_path <- file.path(SIG_DIR, "nested_cv_selection_frequency.csv")
}
if (!file.exists(f_b_path)) {
  f_b_path <- file.path(INPUT_DIR, "nested_cv_selection_frequency.csv")
}
if (!file.exists(f_b_path)) stop("selection frequency file not found.")

log_msg(sprintf("   reading: %s", basename(f_b_path)))
df_b_raw <- read_csv(f_b_path, show_col_types = FALSE)

# miRNA and frequency columns
mi_col <- intersect(c("miRNA", "mirna", "feature", "isoform"),
                    names(df_b_raw))[1]
fr_col <- intersect(c("frequency", "selection_frequency", "freq",
                      "selection_pct", "pct_selected"),
                    names(df_b_raw))[1]
if (any(is.na(c(mi_col, fr_col)))) stop("required Panel B columns not found.")

TOP_N <- 20
SIG_THRESHOLD <- 60   # v2: 60% (not 70%)

df_b <- df_b_raw %>%
  rename(miRNA = all_of(mi_col), freq = all_of(fr_col)) %>%
  mutate(
    freq_pct     = if (max(freq, na.rm = TRUE) <= 1) freq * 100 else freq,
    in_signature = freq_pct >= SIG_THRESHOLD
  ) %>%
  arrange(desc(freq_pct)) %>%
  slice_head(n = TOP_N) %>%
  mutate(miRNA = fct_reorder(miRNA, freq_pct))

n_total_candidates <- nrow(df_b_raw)
n_in_signature <- sum(
  (if (max(df_b_raw[[fr_col]], na.rm = TRUE) <= 1)
     df_b_raw[[fr_col]] * 100
   else df_b_raw[[fr_col]]) >= SIG_THRESHOLD,
  na.rm = TRUE
)

log_msg(sprintf("   Top %d of %d candidates | %d in the signature (\u2265%d%%)",
                TOP_N, n_total_candidates, n_in_signature, SIG_THRESHOLD))

pB <- ggplot(df_b, aes(x = miRNA, y = freq_pct, fill = in_signature)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = SIG_THRESHOLD, linetype = "dashed",
             color = "grey20", linewidth = 0.5) +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "#2166AC", "FALSE" = "grey70"),
    labels = c("TRUE"  = sprintf("In final signature (\u2265%d%%)", SIG_THRESHOLD),
               "FALSE" = "Not selected"),
    name = NULL
  ) +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(
    title    = "miRNA Selection Stability",
    subtitle = sprintf("Top %d of %d candidate isoforms", TOP_N, n_total_candidates),
    x = NULL,
    y = "Selection frequency (%)"
  ) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())

# =============================================================================
# combine A + B with patchwork
# =============================================================================
log_msg("")
log_msg("combining Panel A + B ...")

combined <- pA + pB + plot_layout(widths = c(1, 1.1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

# =============================================================================
# save
# =============================================================================
ggsave(out_path("Fig3_NestedCV_and_Stability_600dpi.png"),
       combined, width = 12, height = 5.5, dpi = 600, bg = "white")
ggsave(out_path("Fig3_NestedCV_and_Stability_600dpi.pdf"),
       combined, width = 12, height = 5.5, device = cairo_pdf, bg = "white")

log_msg("")
log_msg("============================================================")
log_msg("Figure 3 saved:")
log_msg("   Fig3_NestedCV_and_Stability_600dpi.png")
log_msg("   Fig3_NestedCV_and_Stability_600dpi.pdf")
log_msg("============================================================")

# =============================================================================
# Caption

# image paths
img1_path <- file.path(WORK_DIR, "signature_v2/signature_v2_km_median_600dpi.png")
img2_path <- file.path(WORK_DIR, "signature_v2/signature_v2_km_quartile_600dpi.png")
out_path  <- file.path(WORK_DIR, "signature_v2/signature_v2_combined_600dpi.png")

# read the images
img1 <- image_read(img1_path)
img2 <- image_read(img2_path)

# paste side by side (horizontally)
combined <- image_append(c(img1, img2), stack = FALSE)

# save at 600 DPI
image_write(
  combined,
  path = out_path,
  format = "png",
  density = "600x600"   # ← set the DPI when saving
)

# check the final DPI
info <- image_info(combined)
print(info$density)
# =============================================================================
caption <- paste0(
  "Figure 3. Nested cross-validation performance and miRNA selection stability. ",
  "(A) Distribution of outer-fold C-index across valid nested-CV folds (n = ", n_ci,
  "); the dashed vertical line indicates the mean (", sprintf("%.3f", mean_ci),
  " \u00b1 ", sprintf("%.3f", sd_ci), "). ",
  "(B) Selection frequency of the top ", TOP_N, " miRNA isoforms ",
  "(of ", n_total_candidates, " total candidate isoforms). Blue bars denote ",
  "miRNAs included in the final consensus signature (\u2265", SIG_THRESHOLD,
  "% selection frequency); grey bars denote candidates not retained."
)
writeLines(caption, out_path("Fig3_NestedCV_and_Stability_caption.txt"))
log_msg("   caption saved")
