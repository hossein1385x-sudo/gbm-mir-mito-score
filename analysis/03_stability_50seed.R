#!/usr/bin/env Rscript
# =============================================================================
# analysis/03_stability_50seed.R
# Stability of elastic-net/LASSO-Cox feature selection over 50 independent seeds
# TRAIN only (the test set is never touched).
# Output: figures/Fig_stability_jaccard.{pdf,png}
#         supplementary/Table_S_stability_frequency.csv
#         supplementary/Table_S_stability_selection_matrix.csv
# =============================================================================

suppressPackageStartupMessages({
  library(glmnet); library(survival)
  library(dplyr);  library(readr); library(ggplot2)
})

# ---------------------------------------------------------------- 0) Paths
source(here::here("config.R"))   # defines WORK_DIR (see README)
OUT <- WORK_DIR
FIG  <- here::here("figures")
SUPP <- here::here("supplementary")
dir.create(FIG,  showWarnings = FALSE, recursive = TRUE)
dir.create(SUPP, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- 1) Locked signature
sig   <- readRDS(file.path(OUT, "signature_locked_latest.rds"))
scal  <- sig$scaling
ALPHA <- sig$alpha
cat("alpha (locked):", ALPHA, "\n")
cat("Locked features:", paste(sig$coefficients$Feature, collapse = ", "), "\n\n")

# ---------------------------------------------------------------- 2) Train data
tr <- readRDS(file.path(OUT, "train_full_preprocessed.rds"))
X  <- as.matrix(tr$X_full)
y  <- tr$y_full

# Age from train_set.csv (days -> years)
ts <- readr::read_csv(file.path(OUT, "train_set.csv"), show_col_types = FALSE)
ts$Age_years <- ts$age_at_diagnosis / 365.25
idx <- match(rownames(X), ts$patient_barcode)
Age <- ts$Age_years[idx]
if (any(is.na(Age))) {
  cat("[warn] imputing", sum(is.na(Age)), "missing Age with train median\n")
  Age[is.na(Age)] <- median(ts$Age_years, na.rm = TRUE)
}

# z-score Age with the locked parameters
zA <- (Age - scal$clinical$center["Age"]) / scal$clinical$scale["Age"]
X  <- cbind(X, Age = zA)

# NA/Inf handling
# NA/Inf is imputed with the column median (0 is not neutral on the raw miRNA scale)
bad <- !is.finite(X)
if (any(bad)) {
  cat(sprintf("[warn] %d non-finite entries in %d columns -> column-median imputation\n",
              sum(bad), sum(colSums(bad) > 0)))
  for (j in which(colSums(bad) > 0)) X[bad[, j], j] <- median(X[!bad[, j], j])
}

cat("X dim:", paste(dim(X), collapse = "x"), "\n")
cat("Events:", sum(y[, "status"] == 1), " / ", nrow(X), "\n\n")

# ---------------------------------------------------------------- 3) 50 seeds
seeds    <- 1:50
# With SUBSAMPLE_FRAC = 1 the seeds only change the CV folds (data fixed) => stability w.r.t. fold assignment.
# SUBSAMPLE_FRAC < 1 (e.g. 0.8) subsamples patients in every seed = true sampling stability.
SUBSAMPLE_FRAC <- 1
selected <- vector("list", length(seeds))
names(selected) <- paste0("seed_", seeds)

t0 <- Sys.time()
for (i in seq_along(seeds)) {
  s <- seeds[i]
  set.seed(s)
  rows <- if (SUBSAMPLE_FRAC < 1) sort(sample(nrow(X), floor(SUBSAMPLE_FRAC * nrow(X)))) else seq_len(nrow(X))
  fit <- tryCatch(
    glmnet::cv.glmnet(
      x = X[rows, , drop = FALSE], y = y[rows, ],
      family      = "cox",
      alpha       = ALPHA,
      nfolds      = 10,
      standardize = TRUE      # glmnet standardizes the miRNAs itself
    ),
    error = function(e) {
      message("seed ", s, " failed: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(fit)) { selected[[i]] <- character(0); next }

  coefs <- coef(fit, s = "lambda.min")
  nz    <- rownames(coefs)[which(coefs[, 1] != 0)]
  selected[[i]] <- nz
  cat(sprintf("seed %2d: %d features  ->  %s\n",
              s, length(nz), paste(nz, collapse = ", ")))
}
cat(sprintf("\nElapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---------------------------------------------------------------- 4) Selection frequency
all_feats <- sort(unique(unlist(selected)))
if (length(all_feats) == 0) stop("No feature was selected in any seed")

freq <- sapply(all_feats, function(f)
  mean(sapply(selected, function(s) f %in% s)))
freq_df <- data.frame(Feature = names(freq), Frequency = as.numeric(freq)) %>%
  arrange(desc(Frequency))
cat("\n=== Selection frequency (top 15) ===\n")
print(head(freq_df, 15), row.names = FALSE)

# Selection frequency of the locked-signature features themselves
locked_feats <- setdiff(sig$coefficients$Feature, "Age")
locked_freq  <- sapply(locked_feats, function(f) if (f %in% names(freq)) unname(freq[f]) else 0)
cat("\n=== Locked-signature features: selection frequency ===\n"); print(locked_freq)
if (any(locked_freq < 0.5))
  warning("Locked feature(s) selected in <50% of seeds: ", paste(names(locked_freq)[locked_freq < 0.5], collapse = ", "),
          " - the input matrix/candidate space differs from the signature-selection procedure, or the signature is unstable.")

# ---------------------------------------------------------------- 5) Jaccard
pair_jaccard <- function(A, B) {
  if (length(A) == 0 && length(B) == 0) return(NA_real_)
  u <- length(union(A, B)); if (u == 0) return(NA_real_)
  length(intersect(A, B)) / u
}
n <- length(seeds)
jacc <- numeric(0)
for (i in 1:(n - 1)) for (j in (i + 1):n)
  jacc <- c(jacc, pair_jaccard(selected[[i]], selected[[j]]))
jacc <- jacc[!is.na(jacc)]

cat(sprintf("\n=== Jaccard stability (n_pairs=%d) ===\n", length(jacc)))
cat(sprintf("mean   : %.3f\n", mean(jacc)))
cat(sprintf("median : %.3f\n", median(jacc)))
cat(sprintf("IQR    : [%.3f, %.3f]\n",
            quantile(jacc, 0.25), quantile(jacc, 0.75)))
cat(sprintf("min/max: %.3f / %.3f\n", min(jacc), max(jacc)))

# ---------------------------------------------------------------- 6) Outputs
sel_mat <- sapply(seeds, function(s)
  as.integer(all_feats %in% selected[[paste0("seed_", s)]]))
rownames(sel_mat) <- all_feats
colnames(sel_mat) <- paste0("seed_", seeds)
write.csv(data.frame(Feature = rownames(sel_mat), sel_mat, check.names = FALSE),
          file.path(SUPP, "Table_S_stability_selection_matrix.csv"), row.names = FALSE)  # header of the first column
write.csv(freq_df, file.path(SUPP, "Table_S_stability_frequency.csv"),
          row.names = FALSE)

# Summary in a txt file
summ <- c(
  sprintf("50-seed %s of elastic-net Cox (alpha=%.2f) on TRAIN (n=%d, events=%d); subsample frac=%.2f",
          if (SUBSAMPLE_FRAC < 1) "subsampling stability" else "CV-fold-seed stability (data fixed)",
          ALPHA, nrow(X), sum(y[, "status"] == 1), SUBSAMPLE_FRAC),   # alpha = 0.5 => elastic net, not pure LASSO
  sprintf("alpha = %.3f | lambda.min in each seed | 10-fold CV",
          ALPHA),
  "",
  sprintf("Jaccard mean   : %.3f", mean(jacc)),
  sprintf("Jaccard median : %.3f", median(jacc)),
  sprintf("Jaccard IQR    : [%.3f, %.3f]",
          quantile(jacc, 0.25), quantile(jacc, 0.75)),
  sprintf("Jaccard range  : [%.3f, %.3f]", min(jacc), max(jacc)),
  "",
  "Features selected in >=50% of seeds:",
  paste0("  ", freq_df$Feature[freq_df$Frequency >= 0.5],
         "  (", sprintf("%.0f%%", 100 * freq_df$Frequency[freq_df$Frequency >= 0.5]), ")"),
  "",
  "Features selected in 100% of seeds:",
  paste0("  ", freq_df$Feature[freq_df$Frequency == 1]),
  "",
  "Locked-signature features (selection frequency):",
  paste0("  ", names(locked_freq), "  (", sprintf("%.0f%%", 100 * locked_freq), ")")
)
writeLines(summ, file.path(SUPP, "Table_S_stability_summary.txt"))

# Jaccard histogram
jacc_df <- data.frame(Jaccard = jacc)
p <- ggplot(jacc_df, aes(x = Jaccard)) +
  geom_histogram(bins = 30, fill = "#2C7FB8", color = "white") +
  geom_vline(xintercept = mean(jacc), linetype = "dashed",
             color = "red", linewidth = 0.8) +
  annotate("text", x = mean(jacc), y = Inf,
           label = sprintf("mean = %.3f", mean(jacc)),
           hjust = -0.1, vjust = 1.5, color = "red") +
  labs(title = "50-seed stability of feature selection",
       subtitle = "Pairwise Jaccard similarity between selected feature sets",
       x = "Jaccard index", y = "Number of seed pairs") +
  theme_bw()
ggsave(file.path(FIG, "Fig_stability_jaccard.pdf"), p, width = 6.5, height = 4.2)
ggsave(file.path(FIG, "Fig_stability_jaccard.png"), p, width = 6.5, height = 4.2, dpi = 300)

cat("\nDone.\n")
cat("  Figure : figures/Fig_stability_jaccard.pdf / .png\n")
cat("  Tables : supplementary/Table_S_stability_selection_matrix.csv\n")
cat("           supplementary/Table_S_stability_frequency.csv\n")
cat("           supplementary/Table_S_stability_summary.txt\n")
