#!/usr/bin/env Rscript
# =============================================================================
# STEP 18 (Signature v2) — Four Priority Figures
# -----------------------------------------------------------------------------
# Four final figures for the manuscript:
#   Fig 1: Selection Frequency (Top-N)
#   Fig 2: Bootstrap Optimism Distribution
#   Fig 3: Forest Plot (multivariable Cox)
#   Fig 4: Decision Curve Analysis
#
# All PNGs are saved at 600 dpi, plus PDFs.
# Each figure has its own caption .txt file.
# =============================================================================

# 0) Packages -----------------------------------------------------------

source(here::here("config.R"))  # defines WORK_DIR (see config.R)

required_pkgs <- c("ggplot2", "dplyr", "readr", "tidyr", "forcats",
                   "scales", "patchwork", "grid")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(sprintf(
    "The following package(s) are not installed: %s\n  install.packages(c(%s))",
    paste(missing_pkgs, collapse = ", "),
    paste(sprintf('"%s"', missing_pkgs), collapse = ", ")
  ))
}
suppressMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(forcats); library(scales); library(patchwork); library(grid)
})

# 1) Paths (unified, with spaces) ----------------------------------------
OUT_DIR <- paste0(WORK_DIR, "/")
FIG_DIR <- file.path(OUT_DIR, "figures_final")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

path_in  <- function(fname) file.path(OUT_DIR, fname)
path_out <- function(fname) file.path(FIG_DIR, fname)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

# 2) Shared theme -----------------------------------------------------------
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      axis.line   = element_line(color = "black", linewidth = 0.4),
      axis.ticks  = element_line(color = "black", linewidth = 0.3),
      plot.title  = element_text(face = "bold", size = rel(1.05), hjust = 0),
      plot.subtitle = element_text(size = rel(0.85), color = "grey30"),
      legend.position = "bottom",
      legend.title = element_text(size = rel(0.85))
    )
}

# 3) Helper functions --------------------------------------------------------
pick_col <- function(df, candidates, source_file, purpose) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    stop(sprintf(
      "In file '%s', no suitable column was found for '%s'.\n  candidates: %s\n  available columns: %s",
      source_file, purpose,
      paste(candidates, collapse = " / "),
      paste(names(df), collapse = ", ")
    ))
  }
  hit[1]
}

write_caption <- function(fname, text) {
  writeLines(text, path_out(fname))
}

log_msg("============================================================")
log_msg("STEP 18 (Signature v2) — Four Priority Figures")
log_msg("============================================================")


# FIGURE 1 — Selection Frequency (Top-N)

log_msg("")
log_msg("Figure 1 — Selection Frequency")

TOP_N <- 20

f1_raw <- read_csv(path_in("nested_cv_selection_frequency.csv"),
                   show_col_types = FALSE)

mirna_col_f1 <- pick_col(f1_raw,
  c("miRNA", "mirna", "feature", "isoform"),
  "nested_cv_selection_frequency.csv", "miRNA identifier")

freq_col_f1 <- pick_col(f1_raw,
  c("frequency", "selection_frequency", "freq", "selection_pct",
    "pct_selected", "selection_frequency_pct"),
  "nested_cv_selection_frequency.csv", "selection frequency")

f1_data <- f1_raw %>%
  rename(miRNA = all_of(mirna_col_f1), freq = all_of(freq_col_f1)) %>%
  mutate(
    freq_pct = if (max(freq, na.rm = TRUE) <= 1) freq * 100 else freq,
    in_signature = freq_pct >= 60
  ) %>%
  arrange(desc(freq_pct)) %>%
  slice_head(n = TOP_N) %>%
  mutate(miRNA = fct_reorder(miRNA, freq_pct))

n_total_candidates <- nrow(f1_raw)
freq_threshold <- if (max(f1_raw[[freq_col_f1]], na.rm = TRUE) <= 1) 0.6 else 60
n_in_signature <- sum(f1_raw[[freq_col_f1]] >= freq_threshold, na.rm = TRUE)

p1 <- ggplot(f1_data, aes(x = miRNA, y = freq_pct, fill = in_signature)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 60, linetype = "dashed",
             color = "grey20", linewidth = 0.5) +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "#2166AC", "FALSE" = "grey70"),
    labels = c("TRUE" = "In final signature (\u226560%)",
               "FALSE" = "Not selected"),
    name = NULL) +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = sprintf("miRNA Selection Stability (Top %d of %d)",
                    TOP_N, n_total_candidates),
    x = NULL, y = "Selection frequency (%)"
  ) +
  theme_pub(11) + theme(legend.position = "top")

ggsave(path_out("Fig1_SelectionFrequency_600dpi.png"), p1,
       width = 6.5, height = 5.5, dpi = 600, bg = "white")
ggsave(path_out("Fig1_SelectionFrequency_600dpi.pdf"), p1,
       width = 6.5, height = 5.5, bg = "white")

write_caption("Fig1_SelectionFrequency_caption.txt", paste0(
  "Figure 1. Selection stability of the top ", TOP_N, " (of ",
  n_total_candidates, " total candidate) miRNA isoforms across repeated ",
  "nested cross-validation on the training set. Candidate screening and ",
  "Elastic Net–Cox fitting were re-executed independently within each fold. ",
  n_in_signature, " miRNA(s) selected in \u226560% of folds (dashed line) ",
  "were retained in the final consensus signature (blue bars). The full ",
  "selection-frequency table for all ", n_total_candidates, " candidates ",
  "is provided in Supplementary Table S_."
))
log_msg(sprintf("   Fig1 saved (Top %d)", TOP_N))


# FIGURE 2 — Bootstrap Optimism Distribution

log_msg("")
log_msg("Figure 2 — Bootstrap Optimism")

boot_path <- path_in("signature_v2/signature_v2_bootstrap_results.csv")
if (!file.exists(boot_path)) {
  # fallback: old file name
  boot_path <- path_in("bootstrap_optimism_corrected_cindex.csv")
}
if (!file.exists(boot_path)) {
  stop("bootstrap file not found.")
}
log_msg(sprintf("   reading: %s", basename(boot_path)))

boot_raw <- read_csv(boot_path, show_col_types = FALSE)

# --- optimism column ---
optimism_col <- pick_col(boot_raw,
  c("optimism", "bootstrap_optimism", "optimism_value"),
  basename(boot_path), "optimism")

f2_data <- boot_raw %>%
  transmute(optimism = as.numeric(.data[[optimism_col]])) %>%
  filter(is.finite(optimism))

if (nrow(f2_data) < 30) {
  stop(sprintf("only %d optimism values — at least 30 are required.", nrow(f2_data)))
}

mean_optimism    <- mean(f2_data$optimism)
median_optimism  <- median(f2_data$optimism)
optimism_sd      <- sd(f2_data$optimism)
optimism_q025    <- quantile(f2_data$optimism, 0.025)
optimism_q975    <- quantile(f2_data$optimism, 0.975)

# --- apparent C-index ---
apparent_col <- pick_col(boot_raw,
  c("apparent_cindex", "apparent_C", "apparent_c", "C_apparent"),
  basename(boot_path), "apparent C-index (optional)")

apparent_cindex <- unique(as.numeric(boot_raw[[apparent_col]]))
apparent_cindex <- apparent_cindex[is.finite(apparent_cindex)][1]

# --- corrected C-index ---
corrected_cindex <- if (is.finite(apparent_cindex)) {
  apparent_cindex - mean_optimism
} else {
  NA_real_
}

p2 <- ggplot(f2_data, aes(x = optimism)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "grey70",
                 color = "white", linewidth = 0.25) +
  geom_density(linewidth = 0.8) +
  geom_vline(xintercept = mean_optimism, linewidth = 0.7) +
  geom_vline(xintercept = c(optimism_q025, optimism_q975),
             linetype = "dashed", linewidth = 0.5) +
  labs(
    title = "Bootstrap Optimism Distribution",
    subtitle = if (is.finite(corrected_cindex)) {
      sprintf("B = %d | Mean optimism = %.3f | Optimism-corrected C = %.3f",
              nrow(f2_data), mean_optimism, corrected_cindex)
    } else {
      sprintf("B = %d | Mean optimism = %.3f",
              nrow(f2_data), mean_optimism)
    },
    x = "Bootstrap optimism in Harrell's C-index",
    y = "Density"
  ) +
  theme_classic(base_size = 11)

ggsave(path_out("Fig2_BootstrapOptimism_600dpi.png"), p2,
       width = 6.5, height = 5, dpi = 600, bg = "white")
ggsave(path_out("Fig2_BootstrapOptimism_600dpi.pdf"), p2,
       width = 6.5, height = 5, bg = "white")

write_caption("Fig2_BootstrapOptimism_caption.txt", paste0(
  "Figure 2. Distribution of bootstrap optimism in Harrell's C-index ",
  "across ", nrow(f2_data), " bootstrap resamples of the training set. ",
  "Optimism is defined as the difference between model performance on the ",
  "bootstrap sample and on the original training sample. The solid vertical ",
  "line indicates the mean optimism (", sprintf("%.3f", mean_optimism),
  "); dashed vertical lines indicate the empirical 2.5th and 97.5th ",
  "percentiles [", sprintf("%.3f", optimism_q025), ", ",
  sprintf("%.3f", optimism_q975), "]. ",
  if (is.finite(corrected_cindex)) {
    paste0("The resulting optimism-corrected C-index was ",
           sprintf("%.3f", corrected_cindex), ". ")
  } else "",
  "The percentile interval shown here describes the bootstrap optimism ",
  "distribution and should not be interpreted as a formal 95% confidence ",
  "interval for the optimism-corrected C-index."
))
log_msg(sprintf("   Fig2 saved (mean optimism = %.3f)", mean_optimism))


# FIGURE 3 — Forest Plot (multivariable Cox)

log_msg("")
log_msg("Figure 3 — Forest Plot")

f3_raw <- read_csv(path_in("signature_v2/14_3_verhaak_analysis/cox_multivariable.csv"),
                   show_col_types = FALSE)

var_col_f3 <- pick_col(f3_raw, c("term", "variable", "covariate", "predictor"),
  "cox_multivariable.csv", "variable name")
hr_col_f3 <- pick_col(f3_raw,
  c("estimate", "HR", "hr", "hazard_ratio", "exp_coef"),
  "cox_multivariable.csv", "hazard ratio")
lower_col_f3 <- pick_col(f3_raw,
  c("conf.low", "CI_lower", "lower_95", "conf_low", "lower_ci", "HR_lower"),
  "cox_multivariable.csv", "CI lower bound")
upper_col_f3 <- pick_col(f3_raw,
  c("conf.high", "CI_upper", "upper_95", "conf_high", "upper_ci", "HR_upper"),
  "cox_multivariable.csv", "CI upper bound")
pval_col_f3 <- pick_col(f3_raw, c("p.value", "p_value", "pvalue", "P"),
  "cox_multivariable.csv", "p-value")

f3_data <- f3_raw %>%
  rename(variable = all_of(var_col_f3), HR = all_of(hr_col_f3),
         CI_lower = all_of(lower_col_f3), CI_upper = all_of(upper_col_f3),
         p_value = all_of(pval_col_f3)) %>%
  mutate(
    p_label = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)),
    hr_label = sprintf("%.2f", HR),
    ci_label = sprintf("%.2f\u2013%.2f", CI_lower, CI_upper),
    significant = CI_lower > 1 | CI_upper < 1
  ) %>%
  arrange(desc(HR)) %>%
  mutate(variable = fct_reorder(variable, HR),
         y_pos = as.numeric(variable))

p3_plot <- ggplot(f3_data, aes(x = HR, y = variable, color = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.18, linewidth = 0.6) +
  geom_point(size = 2.6, shape = 15) +
  scale_x_log10() +
  scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey50"),
                     guide = "none") +
  labs(title = "Multivariable Cox Model",
       subtitle = "Adjusted for molecular subtype",
       x = "Hazard Ratio (log scale)", y = NULL) +
  theme_pub(11) +
  theme(panel.grid.major.y = element_blank(),
        plot.margin = margin(t = 10, r = 4, b = 10, l = 10))

table_long <- f3_data %>%
  select(variable, y_pos, hr_label, ci_label, p_label) %>%
  pivot_longer(cols = c(hr_label, ci_label, p_label),
               names_to = "col", values_to = "text") %>%
  mutate(x_pos = case_when(
    col == "hr_label" ~ 1,
    col == "ci_label" ~ 2,
    col == "p_label"  ~ 3
  ))

header_df <- data.frame(
  x_pos = 1:3,
  label = c("HR", "95% CI", "p-value"),
  y_pos = max(f3_data$y_pos) + 0.9)

p3_table <- ggplot(table_long, aes(x = x_pos, y = y_pos, label = text)) +
  geom_text(size = 3.2, hjust = 0.5, family = "mono") +
  geom_text(data = header_df, aes(x = x_pos, y = y_pos, label = label),
            fontface = "bold", size = 3.2, hjust = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0.4, 3.6)) +
  scale_y_continuous(limits = c(0.4, max(f3_data$y_pos) + 1.3)) +
  theme_void() +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 0))

p3_combined <- p3_plot + p3_table + plot_layout(widths = c(2.2, 1))

ggsave(path_out("Fig3_ForestPlot_600dpi.png"), p3_combined,
       width = 8.5, height = max(3.5, 0.42 * nrow(f3_data) + 1),
       dpi = 600, bg = "white")
ggsave(path_out("Fig3_ForestPlot_600dpi.pdf"), p3_combined,
       width = 8.5, height = max(3.5, 0.42 * nrow(f3_data) + 1),
       bg = "white")

write_caption("Fig3_ForestPlot_caption.txt", paste0(
  "Figure 3. Forest plot of hazard ratios (HR) with 95% confidence intervals ",
  "from the multivariable Cox proportional hazards model, adjusted for ",
  "Verhaak molecular subtype. Squares denote point estimates; horizontal ",
  "lines denote 95% CIs on a log scale. Red markers indicate covariates for ",
  "which the 95% CI excludes HR = 1. Exact HR, 95% CI, and p-values are ",
  "listed in the adjoining table."
))
log_msg("   Fig3 saved")


# FIGURE 4 — Decision Curve Analysis (DCA)

log_msg("")
log_msg("Figure 4 — DCA")

f4_raw <- read_csv(path_in("signature_v2/17_dca/dca_net_benefit_table.csv"),
                   show_col_types = FALSE)

required_f4 <- c("variable", "label", "n", "threshold",
                 "net_benefit", "time_horizon")
missing_f4 <- setdiff(required_f4, names(f4_raw))
if (length(missing_f4) > 0) {
  stop(sprintf("DCA missing columns: %s\n  Available: %s",
               paste(missing_f4, collapse = ", "),
               paste(names(f4_raw), collapse = ", ")))
}

f4_data <- f4_raw %>%
  mutate(
    threshold   = as.numeric(threshold),
    net_benefit = as.numeric(net_benefit),
    time_horizon = as.character(time_horizon),
    strategy = case_when(
      variable == "all"  ~ "Treat All",
      variable == "none" ~ "Treat None",
      grepl("mir|mito|score|signature",
            paste(variable, label), ignore.case = TRUE) ~ "miR-Mito Score",
      TRUE ~ as.character(label)
    )
  )

if (anyNA(f4_data$threshold))   stop("threshold contains NA.")
if (anyNA(f4_data$net_benefit)) stop("net_benefit contains NA.")
if (any(f4_data$threshold < 0 | f4_data$threshold > 1, na.rm = TRUE)) {
  stop("threshold must be probability in [0, 1].")
}

f4_data <- f4_data %>% mutate(threshold_pct = threshold * 100)

# --- horizon labels (English) ---
f4_data <- f4_data %>%
  mutate(horizon_label = case_when(
    time_horizon %in% c("6-month", "6 month", "6mo", "\u06F6 \u0645\u0627\u0647") ~ "6 months",
    time_horizon %in% c("1-year", "1 year", "1yr", "\u06F1 \u0633\u0627\u0644")  ~ "1 year",
    time_horizon %in% c("2-year", "2 years", "2yr", "\u06F2 \u0633\u0627\u0644")  ~ "2 years",
    TRUE ~ NA_character_
  ))

if (anyNA(f4_data$horizon_label)) {
  print(unique(f4_data$time_horizon))
  stop("unknown time_horizon values — see print above.")
}

f4_data$horizon_label <- factor(f4_data$horizon_label,
                                levels = c("6 months", "1 year", "2 years"))

expected_strategies <- c("miR-Mito Score", "Treat All", "Treat None")
missing_strategies <- setdiff(expected_strategies, unique(f4_data$strategy))
if (length(missing_strategies) > 0) {
  stop(sprintf("missing strategies: %s",
               paste(missing_strategies, collapse = ", ")))
}

dup_check <- f4_data %>%
  count(horizon_label, strategy, threshold) %>%
  filter(n > 1)
if (nrow(dup_check) > 0) {
  print(dup_check)
  stop("duplicate strategy×horizon×threshold rows.")
}

f4_data$strategy <- factor(f4_data$strategy, levels = expected_strategies)
f4_data <- f4_data %>% arrange(horizon_label, strategy, threshold)

dca_colors <- c(
  "miR-Mito Score" = "#2166AC",
  "Treat All"      = "grey40",
  "Treat None"     = "black"
)
dca_linetypes <- c(
  "miR-Mito Score" = "solid",
  "Treat All"      = "dashed",
  "Treat None"     = "dotted"
)

y_min <- min(f4_data$net_benefit, na.rm = TRUE)
y_max <- max(f4_data$net_benefit, na.rm = TRUE)
y_padding <- max(0.02, 0.08 * (y_max - y_min))
y_limits  <- c(y_min - y_padding, y_max + y_padding)

p4 <- ggplot(f4_data, aes(x = threshold_pct, y = net_benefit,
                          color = strategy, linetype = strategy,
                          group = interaction(strategy, horizon_label))) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.45) +
  geom_line(linewidth = 0.8, lineend = "round") +
  facet_wrap(~ horizon_label, nrow = 1, scales = "fixed") +
  scale_color_manual(values = dca_colors, breaks = expected_strategies,
                     drop = FALSE) +
  scale_linetype_manual(values = dca_linetypes, breaks = expected_strategies,
                        drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 50, by = 10), limits = c(0, 50),
                     labels = function(x) sprintf("%.0f", x),
                     expand = expansion(mult = c(0, 0.015))) +
  scale_y_continuous(limits = y_limits,
                     breaks = pretty_breaks(n = 5),
                     labels = function(x) sprintf("%.2f", x),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Decision Curve Analysis",
    subtitle = "Net benefit across threshold probabilities",
    x = "Threshold probability (%)",
    y = "Net benefit",
    color = NULL, linetype = NULL
  ) +
  theme_pub(11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(1.4, "cm"),
    legend.text = element_text(size = 9.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11, color = "black"),
    panel.spacing = unit(1.2, "lines"),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0),
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 9),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

ggsave(path_out("Fig4_DCA_600dpi.png"), p4,
       width = 10, height = 4.8, dpi = 600, bg = "white")
ggsave(path_out("Fig4_DCA_600dpi.pdf"), p4,
       width = 10, height = 4.8, bg = "white")

write_caption("Fig4_DCA_caption.txt", paste0(
  "Figure 4. Decision curve analysis (DCA) of the miR-Mito Score compared ",
  "with treat-all and treat-none strategies across threshold probabilities ",
  "from 0% to 50%. Panels show net benefit at 6 months, 1 year, and 2 years. ",
  "The miR-Mito Score is shown as a solid blue line, the treat-all strategy ",
  "as a dashed grey line, and the treat-none strategy as a dotted black line. ",
  "Net benefit reflects the clinical utility of the prediction model by ",
  "accounting for the relative consequences of true-positive and ",
  "false-positive classifications at each threshold."
))
log_msg("   Fig4 saved")


# summary
log_msg("")
log_msg("============================================================")
log_msg("all 4 figures saved in:")
log_msg("   ", FIG_DIR)
log_msg("   the full caption of each figure is in the adjacent .txt file — not on the image.")
log_msg("============================================================")
