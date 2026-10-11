#!/usr/bin/env Rscript
# ============================================================
#  BRIER — CORRECTED (Graf formula, right orientation)
# ============================================================
rm(list = ls())

source(here::here("config.R"))  # defines WORK_DIR (see config.R)

ROOT <- WORK_DIR
OUT  <- file.path(ROOT, "BRIER_FIX_v2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(survival); library(dplyr); library(readr); library(ggplot2)
})

## ---------- Load ----------
pred <- read_csv(file.path(ROOT, "signature_v2/14_5_brier/test_survival_predictions.csv"),
                 show_col_types = FALSE) %>% as.data.frame()
train_raw <- read_csv(file.path(ROOT, "train_set.csv"),
                      show_col_types = FALSE) %>% as.data.frame()

# Build train time/event
train <- train_raw
train$time  <- ifelse(train$vital_status == "Dead",
                      train$days_to_death,
                      train$days_to_last_follow_up)
train$event <- as.integer(train$vital_status == "Dead")
train <- train[!is.na(train$time) & !is.na(train$event) & train$time > 0,
               c("time","event")]

# Test
test       <- pred
test$time  <- as.numeric(test$survival_time)
test$event <- as.integer(test$event)
test       <- test[!is.na(test$time) & !is.na(test$event), ]

cat(sprintf("train: n=%d events=%d censored=%d\n",
            nrow(train), sum(train$event==1), sum(train$event==0)))
cat(sprintf("test : n=%d events=%d censored=%d\n",
            nrow(test), sum(test$event==1), sum(test$event==0)))

## ---------- Ĝ ----------
# use the train set only if it contains censoring.
# otherwise use the test set or a combination.
build_G <- function(fit_df, label) {
  if (sum(fit_df$event == 0) < 5) {
    warning(sprintf("only %d censored observations in %s — Ĝ is unreliable",
                    sum(fit_df$event == 0), label))
  }
  cen <- survfit(Surv(time, 1 - event) ~ 1, data = fit_df)
  function(u) {
    s <- summary(cen, times = u, extend = TRUE)$surv
    pmax(s, 1e-6)
  }
}

G_train <- build_G(train, "train")
G_test  <- build_G(test,  "test")
G_all   <- build_G(rbind(train[,c("time","event")], test[,c("time","event")]),
                   "train+test")

cat("\nG values (train / test / all):\n")
for (tt in c(180, 365, 730)) {
  cat(sprintf("  G(%d) = %.4f / %.4f / %.4f\n",
              tt, G_train(tt), G_test(tt), G_all(tt)))
}

## ---------- CORRECT Graf Brier ----------
brier_graf <- function(time, event, S_hat, t, G) {
  stopifnot(length(time) == length(event),
            length(time) == length(S_hat))
  survive_past_t <- as.numeric(time > t)
  event_by_t     <- as.numeric(event == 1 & time <= t)

  # term1: survivors → (1 - S_hat)^2 / G(t)
  # term2: events → S_hat^2 / G(T_i)
  term1 <- survive_past_t * (1 - S_hat)^2 / G(t)
  term2 <- event_by_t     * S_hat^2       / G(time)

  mean(term1 + term2)
}

## ---------- Horizons ----------
H <- list(
  `6mo` = list(t=180, Sf="S_full_6mo", Sm="S_mir_6mo"),
  `1yr` = list(t=365, Sf="S_full_1yr", Sm="S_mir_1yr"),
  `2yr` = list(t=730, Sf="S_full_2yr", Sm="S_mir_2yr")
)

## ---------- Compute for both G options ----------
compute_all <- function(Gfun, Gname) {
  do.call(rbind, lapply(names(H), function(hz) {
    h  <- H[[hz]]; t <- h$t
    Sf <- as.numeric(test[[h$Sf]])
    Sm <- as.numeric(test[[h$Sm]])

    km_test <- survfit(Surv(time, event) ~ 1, data = test)
    S_km    <- summary(km_test, times = t, extend = TRUE)$surv

    b_full  <- brier_graf(test$time, test$event, Sf, t, Gfun)
    b_mir   <- brier_graf(test$time, test$event, Sm, t, Gfun)
    b_km    <- brier_graf(test$time, test$event, rep(S_km, nrow(test)), t, Gfun)

    # theoretical reference: S(1-S)
    theo <- S_km * (1 - S_km)

    # calibration slope
    lp  <- qlogis(pmin(pmax(Sf,1e-3), 1-1e-3))
    I   <- as.numeric(test$event == 1 & test$time <= t)
    cal <- tryCatch(glm(I ~ lp, family=binomial), error=function(e) NULL)
    slope <- if (!is.null(cal)) unname(coef(cal)[2]) else NA

    data.frame(
      G_source = Gname,
      horizon = hz, time_days = t,
      n_at_risk = sum(test$time >= t),
      n_events_by_t = sum(test$event==1 & test$time <= t),
      S_km = round(S_km, 4),
      S_km_x_1_minus = round(theo, 4),      # theoretical (correct)
      Brier_full = round(b_full, 4),
      Brier_mirna = round(b_mir, 4),
      Brier_KM_ref = round(b_km, 4),
      IPS_full = round(1 - b_full/b_km, 4),
      IPS_mirna = round(1 - b_mir /b_km, 4),
      cal_slope = round(slope, 3)
    )
  }))
}

res_G_train <- compute_all(G_train, "train")
res_G_test  <- compute_all(G_test,  "test")
res_G_all   <- compute_all(G_all,   "train+test")

cat("\n=========================================================\n")
cat("A) Ĝ from TRAIN\n");   print(res_G_train, row.names=FALSE)
cat("\nB) Ĝ from TEST\n");    print(res_G_test,  row.names=FALSE)
cat("\nC) Ĝ from TRAIN+TEST\n"); print(res_G_all, row.names=FALSE)

## ---------- Save ----------
write.csv(res_G_train, file.path(OUT, "brier_G_train.csv"), row.names=FALSE)
write.csv(res_G_test,  file.path(OUT, "brier_G_test.csv"),  row.names=FALSE)
write.csv(res_G_all,   file.path(OUT, "brier_G_all.csv"),   row.names=FALSE)
cat("\nsaved to:", OUT, "\n")
