#!/usr/bin/env Rscript
# ============================================================
#  Figure S2 — Calibration of miR-Mito Score (locked test set)
#  Publication-quality, 600 DPI, PNG + PDF
# ============================================================
rm(list = ls())

source(here::here("config.R"))  # defines WORK_DIR (see config.R)


ROOT <- WORK_DIR
OUT  <- file.path(ROOT, "figures_final")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

## ---------- 1) Load test predictions ----------
pred_file <- file.path(ROOT, "signature_v2/14_5_brier/test_survival_predictions.csv")
stopifnot(file.exists(pred_file))

pred <- read_csv(pred_file, show_col_types = FALSE) %>% as.data.frame()

test <- pred
test$time  <- as.numeric(test$survival_time)
test$event <- as.integer(test$event)
test <- test[!is.na(test$time) & !is.na(test$event), ]

cat(sprintf("Test set: n=%d | events=%d\n",
            nrow(test), sum(test$event == 1)))

## ---------- 2) Horizon specs ----------
H <- list(
  "6 months" = list(t = 180, S = "S_full_6mo"),
  "1 year"   = list(t = 365, S = "S_full_1yr"),
  "2 years"  = list(t = 730, S = "S_full_2yr")
)

## ---------- 3) Build calibration bins (quartiles) with 95% CI ----------
bins <- do.call(rbind, lapply(names(H), function(hz) {
  h  <- H[[hz]]
  t  <- h$t
  S  <- as.numeric(test[[h$S]])

  # quartile bins
  q <- cut(S, breaks = quantile(S, probs = 0:4/4),
           include.lowest = TRUE, labels = FALSE)

  do.call(rbind, lapply(sort(unique(q)), function(k) {
    idx <- which(q == k)
    km  <- survfit(Surv(time, event) ~ 1, data = test[idx, ])
    ss  <- summary(km, times = t, extend = TRUE)
    data.frame(
      horizon   = hz,
      bin       = k,
      n         = length(idx),
      pred_mean = mean(S[idx]),
      obs       = ss$surv[1],
      lo        = if (!is.null(ss$lower)) ss$lower[1] else NA_real_,
      hi        = if (!is.null(ss$upper)) ss$upper[1] else NA_real_
    )
  }))
}))

# lock horizon order
bins$horizon <- factor(bins$horizon,
                       levels = c("6 months", "1 year", "2 years"))

cat("\nCalibration bins:\n")
print(bins, row.names = FALSE, digits = 3)

## ---------- 4) Compose the figure ----------
p <- ggplot(bins, aes(x = pred_mean, y = obs)) +
  # perfect-calibration diagonal
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50", linewidth = 0.5) +
  # CI error bars
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.025, linewidth = 0.55,
                color = "#1f77b4") +
  # connecting line
  geom_line(color = "#1f77b4", linewidth = 0.7) +
  # points
  geom_point(size = 2.6, shape = 21,
             fill = "#1f77b4", color = "white", stroke = 0.6) +
  # panel per horizon
  facet_wrap(~ horizon, nrow = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(
    title = NULL,               # caption goes in the legend, not the panel
    x     = "Predicted survival probability",
    y     = "Observed survival (Kaplan\u2013Meier, 95% CI)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"),
    panel.border     = element_rect(colour = "black", linewidth = 0.6),
    strip.background = element_rect(fill = "grey95", colour = "black",
                                    linewidth = 0.4),
    strip.text       = element_text(face = "bold", size = 11),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10, colour = "black"),
    plot.margin      = margin(8, 10, 8, 8)
  )

## ---------- 5) Save at 600 DPI ----------
png_file <- file.path(OUT, "FigS2_Calibration_600dpi.png")
pdf_file <- file.path(OUT, "FigS2_Calibration_600dpi.pdf")

ggsave(png_file, p, width = 9.5, height = 3.6,
       dpi = 600, units = "in", bg = "white")
ggsave(pdf_file, p, width = 9.5, height = 3.6,
       units = "in", bg = "white", device = cairo_pdf)

cat("\nSaved:\n  ", png_file, "\n  ", pdf_file, "\n")

## ---------- 6) Also save the underlying bin data ----------
write.csv(bins, file.path(OUT, "FigS2_calibration_bins.csv"),
          row.names = FALSE)
cat("  ", file.path(OUT, "FigS2_calibration_bins.csv"), "\n")
