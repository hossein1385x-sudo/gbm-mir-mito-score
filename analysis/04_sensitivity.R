#!/usr/bin/env Rscript
# =============================================================================
# analysis/04_sensitivity.R
# Four sensitivity analyses of the miR-Mito model (fitted on train only; test stays locked)
#
#   1) Trim short follow-up: OS_time >= 30 days
#   2) Cox ties: breslow vs efron
#   3) Alternative splits: 50 random 70/30 splits, C-index distribution
#   4) Remove Age: 2-miRNA-only model
#
# Output: supplementary/Table_S_sensitivity.csv
#         supplementary/Table_S_sensitivity_split.csv
#         figures/Fig_sensitivity.{pdf,png}
# =============================================================================

suppressPackageStartupMessages({
  library(survival); library(glmnet); library(survminer)
  library(dplyr);   library(readr);  library(ggplot2)
})

# =============================================================================
# C-index with fixed direction (risk score: higher = worse).
# =============================================================================
cindex_fixed <- function(time, event, score) {
  ok <- is.finite(time) & is.finite(event) & is.finite(score) & time > 0
  if (sum(ok) < 10 || length(unique(event[ok])) < 2) return(NA_real_)
  df <- data.frame(t = time[ok], e = event[ok], s = score[ok])
  as.numeric(survival::concordance(survival::Surv(t, e) ~ s, data = df,
                                   reverse = TRUE)$concordance)
}

# ---------------------------------------------------------------- 0) Paths
source(here::here("config.R"))   # defines WORK_DIR (see README)
OUT <- WORK_DIR
FIG  <- here::here("figures")
SUPP <- here::here("supplementary")
dir.create(FIG,  showWarnings = FALSE, recursive = TRUE)
dir.create(SUPP, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- 1) Locked signature
sig  <- readRDS(file.path(OUT, "signature_locked_latest.rds"))
scal <- sig$scaling
ALPHA <- sig$alpha
beta <- setNames(sig$coefficients$Coefficient, sig$coefficients$Feature)

# ---------------------------------------------------------------- 2) Train
tr <- readRDS(file.path(OUT, "train_full_preprocessed.rds"))
X  <- as.matrix(tr$X_full)
y  <- tr$y_full

ts <- readr::read_csv(file.path(OUT, "train_set.csv"), show_col_types = FALSE)
ts$Age_years <- ts$age_at_diagnosis / 365.25
idx <- match(rownames(X), ts$patient_barcode)
Age <- ts$Age_years[idx]
if (any(is.na(Age))) Age[is.na(Age)] <- median(ts$Age_years, na.rm = TRUE)

zA <- (Age - scal$clinical$center["Age"]) / scal$clinical$scale["Age"]

# z-score miRNAs with the locked parameters
z1 <- (X[, "hsa-miR-196a-5p"] - scal$mirna$center["hsa-miR-196a-5p"]) /
  scal$mirna$scale ["hsa-miR-196a-5p"]
z2 <- (X[, "hsa-miR-675-5p"]  - scal$mirna$center["hsa-miR-675-5p"])  /
  scal$mirna$scale ["hsa-miR-675-5p"]

lp <- beta["hsa-miR-196a-5p"]*z1 + beta["hsa-miR-675-5p"]*z2 + beta["Age"]*zA

train_df <- data.frame(
  sample_id = rownames(X),
  score     = as.numeric(lp),
  score_miRNA_only = as.numeric(beta["hsa-miR-196a-5p"]*z1 + beta["hsa-miR-675-5p"]*z2),
  z1 = z1, z2 = z2, zA = zA,
  Age = Age,
  OS_time  = as.numeric(y[, "time"]),
  OS_event = as.integer(y[, "status"]),
  stringsAsFactors = FALSE
)
train_df <- train_df %>%
  filter(!is.na(score), !is.na(OS_time), !is.na(OS_event), OS_time > 0)

# ---------------------------------------------------------------- 3) Test + CGGA
d2 <- readr::read_csv(file.path(OUT, "supplementary_final", "DataS2_risk_scores.csv"),
                      show_col_types = FALSE)
mk <- function(df, label) df %>% transmute(
  sample_id = patient_id,
  score     = as.numeric(risk_score),
  score_miRNA_only = as.numeric(risk_score_miRNA_only),
  OS_time   = as.numeric(survival_time),
  OS_event  = as.integer(event),
  cohort    = label)
test_df <- mk(d2 %>% filter(cohort == "TCGA_Test"), "Test")
cgga_df <- mk(d2 %>% filter(cohort == "CGGA"),      "CGGA")

test_df <- test_df %>% filter(!is.na(score), OS_time > 0)
cgga_df <- cgga_df %>% filter(!is.na(score), OS_time > 0)

cat(sprintf("Cohort sizes -> Train %d | Test %d | CGGA %d\n\n",
            nrow(train_df), nrow(test_df), nrow(cgga_df)))
cat(sprintf("OS_time median -> Train %.0f | Test %.0f | CGGA %.0f\n",
            median(train_df$OS_time), median(test_df$OS_time), median(cgga_df$OS_time)))
if (median(cgga_df$OS_time) < 0.2 * median(train_df$OS_time) ||
    median(cgga_df$OS_time) > 5   * median(train_df$OS_time))
  warning("CGGA OS_time is not on the same scale as Train; the >=30d trim would be meaningless.")
cat(sprintf("non-NA miRNA-only score -> Test %d/%d | CGGA %d/%d\n",
            sum(!is.na(test_df$score_miRNA_only)), nrow(test_df),
            sum(!is.na(cgga_df$score_miRNA_only)), nrow(cgga_df)))
if (all(is.na(cgga_df$score_miRNA_only)))
  warning("risk_score_miRNA_only is all NA in CGGA. The CGGA 'main' column is also not equivalent to Test/Train - check the script that builds DataS2.")

# ---------------------------------------------------------------- Helper: C-index
cindex_cox <- function(d, score_col = "score") {
  if (!score_col %in% names(d)) return(NA_real_)
  if (all(is.na(d[[score_col]]))) return(NA_real_)
  cindex_fixed(d$OS_time, d$OS_event, d[[score_col]])
}

# =============================================================================
# C-index of the main model (locked score) in the three cohorts.
# (Required by main_tbl below.)
# cindex_fixed is rank-based and therefore insensitive to ties.
# =============================================================================
c_main <- c(
  Train = cindex_cox(train_df, "score"),
  Test  = cindex_cox(test_df,  "score"),
  CGGA  = cindex_cox(cgga_df,  "score"))
cat("Main (locked score):",
    paste(sprintf("%s=%.3f", names(c_main), c_main), collapse = "  "), "\n")

# ---------------------------------------------------------------- 1) Trim short follow-up
cat("\n=== 1) Trim follow-up (OS_time >= 30 days) ===\n")
tr30 <- train_df %>% filter(OS_time >= 30)
te30 <- test_df  %>% filter(OS_time >= 30)
cg30 <- cgga_df  %>% filter(OS_time >= 30)
cat(sprintf("n after trim -> Train %d | Test %d | CGGA %d\n",
            nrow(tr30), nrow(te30), nrow(cg30)))

c_trim <- c(
  Train = cindex_cox(tr30, "score"),
  Test  = cindex_cox(te30, "score"),
  CGGA  = cindex_cox(cg30, "score"))

# ---------------------------------------------------------------- 2) Cox ties
cat("\n=== 2) Cox ties: breslow vs efron (REFIT in TRAIN) ===\n")
f_ef <- survival::coxph(survival::Surv(OS_time, OS_event) ~ z1 + z2 + zA, data = train_df, ties = "efron")
f_br <- survival::coxph(survival::Surv(OS_time, OS_event) ~ z1 + z2 + zA, data = train_df, ties = "breslow")
ties_tbl <- data.frame(term = names(coef(f_ef)), beta_efron = coef(f_ef), beta_breslow = coef(f_br),
                       beta_locked = as.numeric(beta[c("hsa-miR-196a-5p","hsa-miR-675-5p","Age")]))
print(ties_tbl, row.names = FALSE)
readr::write_csv(ties_tbl, file.path(SUPP, "Table_S_ties_coefficients.csv"))
c_efron_refit   <- cindex_fixed(train_df$OS_time, train_df$OS_event, predict(f_ef, type = "lp"))
c_breslow_refit <- cindex_fixed(train_df$OS_time, train_df$OS_event, predict(f_br, type = "lp"))
cat(sprintf("Train C (refit): efron=%.3f  breslow=%.3f\n", c_efron_refit, c_breslow_refit))
# Test/CGGA: miRNA z-scores are not in DataS2 => not computable (intentional NA, not an error)
c_breslow <- c(Train = c_breslow_refit, Test = NA_real_, CGGA = NA_real_)

# ---------------------------------------------------------------- 3) 50 alternative splits
cat("\n=== 3) 50x random 70/30 hold-out WITHIN TRAIN, full pipeline refit ===\n")
keep <- match(train_df$sample_id, rownames(X))
Xall <- cbind(X[keep, , drop = FALSE], Age = train_df$zA)
badx <- !is.finite(Xall)
for (j in which(colSums(badx) > 0)) Xall[badx[, j], j] <- median(Xall[!badx[, j], j])
yy   <- y[keep, ]
n_rep <- 50
c_vec <- rep(NA_real_, n_rep); n_sel <- integer(n_rep)
set.seed(2026)
for (i in seq_len(n_rep)) {
  idx_tr <- sample(nrow(Xall), floor(0.7 * nrow(Xall)))
  fit <- tryCatch(glmnet::cv.glmnet(Xall[idx_tr, ], yy[idx_tr, ], family = "cox",
                                    alpha = ALPHA, nfolds = 10, standardize = TRUE),
                  error = function(e) NULL)
  if (is.null(fit)) next
  n_sel[i] <- sum(as.numeric(coef(fit, s = "lambda.min")) != 0)
  if (n_sel[i] == 0) next
  lp_i <- as.numeric(predict(fit, newx = Xall[-idx_tr, ], s = "lambda.min", type = "link"))
  c_vec[i] <- cindex_fixed(yy[-idx_tr, "time"], yy[-idx_tr, "status"], lp_i)
}
cat(sprintf("splits with empty model (dropped): %d | valid: %d\n", sum(n_sel == 0), sum(!is.na(c_vec))))
c_vec <- c_vec[!is.na(c_vec)]
cat(sprintf("C-index over %d splits -> mean=%.3f  median=%.3f  sd=%.3f  range=[%.3f, %.3f]\n",
            length(c_vec), mean(c_vec), median(c_vec), sd(c_vec), min(c_vec), max(c_vec)))

# ---------------------------------------------------------------- 4) Remove Age
cat("\n=== 4) Remove Age (2-miRNA-only model) ===\n")
c_mirna <- c(
  Train = cindex_cox(train_df, "score_miRNA_only"),
  Test  = cindex_cox(test_df,  "score_miRNA_only"),
  CGGA  = cindex_cox(cgga_df,  "score_miRNA_only"))
cat("miRNA-only:", paste(sprintf("%s=%.3f", names(c_mirna), c_mirna), collapse = "  "), "\n")

# ---------------------------------------------------------------- Main table
main_tbl <- data.frame(
  analysis = c(
    "Main model (locked score)",
    "Trim follow-up >= 30d",
    "Ties = breslow (refit)",
    "Remove Age (2-miRNA only)"),
  C_index_Train = round(c(c_main["Train"], c_trim["Train"],
                          c_breslow["Train"], c_mirna["Train"]), 3),
  C_index_Test  = round(c(c_main["Test"],  c_trim["Test"],
                          c_breslow["Test"],  c_mirna["Test"]), 3),
  C_index_CGGA  = round(c(c_main["CGGA"],  c_trim["CGGA"],
                          c_breslow["CGGA"],  c_mirna["CGGA"]), 3),
  note = c("locked score; fixed-direction C",
           "locked score; fixed-direction C",
           "Train = refit (breslow); Test/CGGA not computable (no miRNA z-scores in DataS2)",
           "miRNA-only score from DataS2; NA = column missing for that cohort"),
  stringsAsFactors = FALSE)
readr::write_csv(main_tbl, file.path(SUPP, "Table_S_sensitivity.csv"))
cat("\n=== Table_S_sensitivity.csv ===\n"); print(main_tbl, row.names = FALSE)

# ---------------------------------------------------------------- Split table
split_tbl <- data.frame(
  metric = c("n_splits", "mean_C", "median_C", "sd_C", "min_C", "max_C"),
  value  = round(c(length(c_vec), mean(c_vec), median(c_vec),
                   sd(c_vec), min(c_vec), max(c_vec)), 3))
readr::write_csv(split_tbl, file.path(SUPP, "Table_S_sensitivity_split.csv"))

# ---------------------------------------------------------------- Figure
p_split <- ggplot(data.frame(C = c_vec), aes(x = C)) +
  geom_histogram(bins = 15, fill = "#2C7FB8", color = "white") +
  geom_vline(xintercept = median(c_vec), color = "red",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = median(c_vec), y = Inf,
           label = sprintf("median C = %.3f", median(c_vec)),
           hjust = -0.05, vjust = 1.5, color = "red", size = 3.3) +
  labs(title = "Sensitivity: 50x random 70/30 splits (within TRAIN, full refit)",
       subtitle = "Distribution of test C-index across alternative splits",
       x = "C-index", y = "Number of splits") +
  theme_bw()
ggsave(file.path(FIG, "Fig_sensitivity.pdf"), p_split, width = 6.5, height = 4.2)
ggsave(file.path(FIG, "Fig_sensitivity.png"), p_split, width = 6.5, height = 4.2, dpi = 300)

cat("\nDone.\n")
cat("  Tables : supplementary/Table_S_sensitivity.csv\n")
cat("           supplementary/Table_S_sensitivity_split.csv\n")
cat("  Figure : figures/Fig_sensitivity.pdf / .png\n")
