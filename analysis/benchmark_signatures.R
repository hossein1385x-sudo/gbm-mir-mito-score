#!/usr/bin/env Rscript
# =============================================================================
# benchmark_signatures.R
# Benchmark of your model against published signatures (fixed published coefficients, no refitting)
#
# Outputs:
#   coverage      : fraction of each signature's miRNAs available in each cohort
#   performance   : C-index (fixed direction: higher score = higher risk) + bootstrap CI
#   delta         : paired difference C(your model) - C(competitor), bootstrap CI over the same patients
#   null (optional): C distribution of random size-matched signatures (fitted on the same Train)
#
# Requires only the survival package (base R + survival).
#
# Input cohorts: named list; each element is a list with:
#   expr  : sample x miRNA matrix, rownames = patient id (log2 scale comparable to Train)
#   time  : survival time (same order as expr rows)      event : 0/1
#   age   : (optional) age in years
# Input own_scores: named list (same names as cohorts); numeric vectors named by patient id
# Input signatures: data.frame with columns signature, feature, coef, [alt_names], [note]
#   (higher = more risk; take coefficients verbatim from the paper)
# =============================================================================

suppressPackageStartupMessages(library(survival))

# ---------------------------------------------------------------- C-index with fixed direction
cidx <- function(time, event, score) {
  ok <- is.finite(time) & is.finite(event) & is.finite(score)
  if (sum(ok) < 10 || sum(event[ok] == 1) < 5) return(NA_real_)
  df <- data.frame(t = time[ok], e = event[ok], s = score[ok])
  as.numeric(survival::concordance(survival::Surv(t, e) ~ s, data = df,
                                   reverse = TRUE)$concordance)
}

# ---------------------------------------------------------------- Signatures
read_signatures <- function(path) {
  s <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("signature", "feature", "coef")
  if (!all(need %in% names(s))) stop("signatures CSV needs columns: ", paste(need, collapse = ", "))
  if (!"alt_names" %in% names(s)) s$alt_names <- ""
  s$alt_names[is.na(s$alt_names)] <- ""
  s
}

# Name matching (case-insensitive). Exact name or alt_names only; the 3p/5p arm is never guessed.
resolve_features <- function(sig_df, available) {
  al <- tolower(available)
  matched <- vapply(seq_len(nrow(sig_df)), function(i) {
    alts <- if (nzchar(sig_df$alt_names[i])) trimws(strsplit(sig_df$alt_names[i], ";")[[1]]) else character(0)
    cn <- tolower(unique(c(sig_df$feature[i], alts)))
    hit <- match(cn, al); hit <- hit[!is.na(hit)]
    if (length(hit)) available[hit[1]] else NA_character_
  }, character(1))
  data.frame(feature = sig_df$feature, coef = sig_df$coef, matched = matched,
             stringsAsFactors = FALSE)
}

score_signature <- function(expr, res, variant = c("as_published", "zscore_within_cohort")) {
  variant <- match.arg(variant)
  m <- res[!is.na(res$matched), , drop = FALSE]
  X <- expr[, m$matched, drop = FALSE]
  if (variant == "zscore_within_cohort") X <- scale(X)
  s <- drop(X %*% m$coef); names(s) <- rownames(expr); s
}

check_scale <- function(cohorts) {
  cat("--- expression scale check ---\n")
  for (nm in names(cohorts)) {
    x <- cohorts[[nm]]$expr
    if (is.null(x) || length(x) == 0 || all(is.na(x))) { cat(sprintf("%-8s EMPTY / all-NA — skipped\n", nm)); next }
    q <- quantile(x, c(0, .5, .99, 1), na.rm = TRUE)
    cat(sprintf("%-8s n=%d p=%d | min=%.2f median=%.2f p99=%.2f max=%.2f | NA=%.1f%%\n",
                nm, nrow(x), ncol(x), q[1], q[2], q[3], q[4], 100 * mean(is.na(x))))
    if (q[4] > 50) warning(nm, ": max expression > 50 — probably NOT log2-scale; apply log2(x+1) consistently in ALL cohorts.")
  }
  invisible(NULL)
}

# ---------------------------------------------------------------- Main run
run_benchmark <- function(cohorts, own_scores, signatures,
                          train_cohort = NULL,        # if given: a 'signature + Age' version with coefficients refitted on Train is also built
                          eval_cohorts = NULL,        # default: all except train_cohort
                          include_age_baseline = TRUE,
                          variants = c("as_published", "zscore_within_cohort"),
                          extra_scores = NULL,        # list(cohort = list(NAME = named numeric vector)): additional models of your own (e.g. OWN_miRNA_only)
                          min_coverage = 0.8, B = 1000, out_dir = NULL, seed = 1) {
  set.seed(seed)
  check_scale(cohorts)
  if (is.null(eval_cohorts)) eval_cohorts <- setdiff(names(cohorts), train_cohort)
  sig_names <- unique(signatures$signature)

  # ---- 1) score of each signature in each cohort + coverage
  raw <- list(); cov_rows <- list()
  for (cn in names(cohorts)) {
    ex <- cohorts[[cn]]$expr; raw[[cn]] <- list()
    for (sn in sig_names) {
      sd_ <- signatures[signatures$signature == sn, , drop = FALSE]
      avail <- colnames(ex)[colSums(!is.na(ex)) > 0]      # all-NA columns do not count as available
      res <- resolve_features(sd_, avail)
      cov_n <- mean(!is.na(res$matched))
      cov_w <- sum(abs(res$coef[!is.na(res$matched)])) / sum(abs(res$coef))
      cov_rows[[length(cov_rows) + 1]] <- data.frame(
        cohort = cn, signature = sn, n_features = nrow(res), n_matched = sum(!is.na(res$matched)),
        coverage = round(cov_n, 3), coverage_weighted = round(cov_w, 3),
        missing = paste(res$feature[is.na(res$matched)], collapse = ";"),
        used = cov_n >= min_coverage, stringsAsFactors = FALSE)
      if (cov_n >= min_coverage)
        for (v in variants) raw[[cn]][[paste0(sn, " [", v, "]")]] <- score_signature(ex, res, v)
    }
  }
  coverage <- do.call(rbind, cov_rows)
  cat("\n--- coverage ---\n"); print(coverage[, c("cohort","signature","n_matched","n_features","coverage","used")], row.names = FALSE)

  # ---- 2) signature + Age version (coefficients refitted on Train): fair comparison with a model that includes Age
  plus_age <- list()
  if (!is.null(train_cohort)) {
    tr <- cohorts[[train_cohort]]
    if (is.null(tr$age)) {
      warning("train_cohort has no age; skipping 'signature + Age' variants.")
    } else {
      mu_a <- mean(tr$age, na.rm = TRUE); sd_a <- sd(tr$age, na.rm = TRUE)
      for (col in names(raw[[train_cohort]])) {
        s_tr <- raw[[train_cohort]][[col]]; mu_s <- mean(s_tr, na.rm = TRUE); sd_s <- sd(s_tr, na.rm = TRUE)
        d <- data.frame(t = tr$time, e = tr$event, s = (s_tr - mu_s) / sd_s, a = (tr$age - mu_a) / sd_a)
        fit <- tryCatch(survival::coxph(survival::Surv(t, e) ~ s + a, data = d), error = function(e) NULL)
        if (is.null(fit)) next
        b <- coef(fit)
        for (cn in eval_cohorts) {
          if (is.null(cohorts[[cn]]$age) || is.null(raw[[cn]][[col]])) next
          s <- raw[[cn]][[col]]
          plus_age[[cn]][[sub(" \\[(.*)\\]$", " + Age(refit on Train) [\\1]", col)]] <-
            b["s"] * (s - mu_s) / sd_s + b["a"] * (cohorts[[cn]]$age - mu_a) / sd_a
        }
      }
    }
  }

  # ---- 3) evaluation in each cohort
  perf <- list(); delt <- list(); keep_scores <- list()
  for (cn in eval_cohorts) {
    co <- cohorts[[cn]]; ids <- rownames(co$expr)
    if (is.null(ids)) stop(cn, ": expr needs rownames (patient ids)")
    own <- own_scores[[cn]]
    if (is.null(own)) { warning(cn, ": no own_scores; skipping"); next }
    if (is.null(names(own))) stop(cn, ": own_scores must be a NAMED vector (names = patient ids) to avoid misalignment")
    mm <- match(ids, names(own))
    if (anyNA(mm)) warning(cn, ": ", sum(is.na(mm)), " patients have no own score; they are dropped")
    S <- list(OWN = as.numeric(own[mm]))
    if (!is.null(extra_scores[[cn]]))
      for (nm in names(extra_scores[[cn]])) {
        v <- extra_scores[[cn]][[nm]]
        if (is.null(names(v))) stop(cn, ": extra score '", nm, "' must be a NAMED vector")
        S[[nm]] <- as.numeric(v[match(ids, names(v))])
      }
    if (include_age_baseline && !is.null(co$age)) S[["AGE_ONLY"]] <- as.numeric(co$age)
    for (col in names(raw[[cn]])) S[[col]] <- as.numeric(raw[[cn]][[col]])
    for (col in names(plus_age[[cn]])) S[[col]] <- as.numeric(plus_age[[cn]][[col]])
    S <- as.data.frame(S, check.names = FALSE)
    ok <- is.finite(co$time) & co$time > 0 & co$event %in% c(0, 1) & complete.cases(S)
    if (sum(ok) < nrow(S)) cat(sprintf("[%s] common complete cases: %d / %d\n", cn, sum(ok), nrow(S)))
    S <- S[ok, , drop = FALSE]; t_ <- co$time[ok]; e_ <- co$event[ok]; n <- nrow(S)
    cat(sprintf("\n[%s] n=%d events=%d | models compared: %d | bootstrap B=%d\n",
                cn, n, sum(e_), ncol(S) - 1, B))

    C0 <- vapply(S, function(s) cidx(t_, e_, s), numeric(1))
    bs <- matrix(NA_real_, B, ncol(S), dimnames = list(NULL, names(S)))
    for (b in seq_len(B)) {
      i <- sample.int(n, n, replace = TRUE)
      if (sum(e_[i]) < 5) next
      for (j in seq_len(ncol(S))) bs[b, j] <- cidx(t_[i], e_[i], S[i, j])
    }
    q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
    perf[[cn]] <- data.frame(
      cohort = cn, model = names(S), n = n, events = sum(e_), C = round(C0, 3),
      C_lo = round(apply(bs, 2, q, .025), 3), C_hi = round(apply(bs, 2, q, .975), 3),
      stringsAsFactors = FALSE, row.names = NULL)
    perf[[cn]]$inverted_flag <- perf[[cn]]$C_hi < 0.5   # CI entirely below 0.5: score direction is probably inverted
    others <- setdiff(names(S), "OWN")
    delt[[cn]] <- do.call(rbind, lapply(others, function(m) {
      d <- bs[, "OWN"] - bs[, m]
      data.frame(cohort = cn, comparison = paste0("OWN - ", m),
                 delta_C = round(C0["OWN"] - C0[m], 3),
                 lo = round(q(d, .025), 3), hi = round(q(d, .975), 3),
                 p_boot_2sided = round(min(1, 2 * min(mean(d <= 0, na.rm = TRUE), mean(d >= 0, na.rm = TRUE))), 3),
                 stringsAsFactors = FALSE, row.names = NULL)
    }))
    keep_scores[[cn]] <- S
  }
  performance <- do.call(rbind, perf); delta <- do.call(rbind, delt)
  rownames(performance) <- rownames(delta) <- NULL
  cat("\n--- performance (fixed direction; CI = patient bootstrap) ---\n"); print(performance, row.names = FALSE)
  cat("\n--- paired difference OWN - other ---\n"); print(delta, row.names = FALSE)

  out <- list(coverage = coverage, performance = performance, delta = delta, scores = keep_scores)
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(coverage,    file.path(out_dir, "benchmark_coverage.csv"),    row.names = FALSE)
    write.csv(performance, file.path(out_dir, "benchmark_performance.csv"), row.names = FALSE)
    write.csv(delta,       file.path(out_dir, "benchmark_delta.csv"),       row.names = FALSE)
    plot_benchmark(performance, file.path(out_dir, "benchmark_forest.pdf"))
    cat("\nSaved to:", out_dir, "\n")
  }
  invisible(out)
}

# ---------------------------------------------------------------- Forest plot (base R)
plot_benchmark <- function(performance, file) {
  pdf(file, width = 8, height = 6)
  for (cn in unique(performance$cohort)) {
    d <- performance[performance$cohort == cn, ]
    d <- d[order(d$C), ]
    par(mar = c(4, 22, 3, 1))
    plot(NA, xlim = c(min(0.4, d$C_lo, na.rm = TRUE), max(0.9, d$C_hi, na.rm = TRUE)),
         ylim = c(0.5, nrow(d) + 0.5), yaxt = "n", xlab = "C-index (95% bootstrap CI)", ylab = "",
         main = paste0("Benchmark — ", cn))
    abline(v = 0.5, lty = 2, col = "grey50")
    cols <- ifelse(d$model == "OWN", "#D95F02", ifelse(d$model == "AGE_ONLY", "grey40", "#2C7FB8"))
    segments(d$C_lo, seq_len(nrow(d)), d$C_hi, seq_len(nrow(d)), col = cols, lwd = 2)
    points(d$C, seq_len(nrow(d)), pch = 19, col = cols)
    axis(2, at = seq_len(nrow(d)), labels = d$model, las = 1, cex.axis = 0.65)
  }
  dev.off()
}

# ---------------------------------------------------------------- null: random size-matched signatures
# Each draw: k random miRNAs from the pool, z-scored with Train mean/SD, coxph on Train (+ optional Age), applied to the target cohorts.
# Answers: does your selection beat k random miRNAs fitted by the same method?
random_signature_null <- function(train, targets, k, own_C = NULL, n_draw = 500,
                                  include_age = TRUE, pool = NULL, seed = 1) {
  set.seed(seed)
  Xtr <- train$expr
  ok_f <- colnames(Xtr)[colSums(is.na(Xtr)) == 0 & apply(Xtr, 2, sd, na.rm = TRUE) > 0]
  for (tg in targets) ok_f <- intersect(ok_f, colnames(tg$expr))
  if (!is.null(pool)) ok_f <- intersect(ok_f, pool)
  if (length(ok_f) < k) stop("pool smaller than k")
  mu <- colMeans(Xtr[, ok_f]); sdv <- apply(Xtr[, ok_f], 2, sd)
  use_age <- include_age && !is.null(train$age) && all(vapply(targets, function(t) !is.null(t$age), logical(1)))
  if (use_age) { mu_a <- mean(train$age, na.rm = TRUE); sd_a <- sd(train$age, na.rm = TRUE) }
  zf <- function(M, f) sweep(sweep(M[, f, drop = FALSE], 2, mu[f], "-"), 2, sdv[f], "/")
  res <- vector("list", n_draw)
  for (d in seq_len(n_draw)) {
    f <- sample(ok_f, k)
    Z <- zf(Xtr, f); if (use_age) Z <- cbind(Z, AGE = (train$age - mu_a) / sd_a)
    fit <- tryCatch(suppressWarnings(survival::coxph(survival::Surv(train$time, train$event) ~ Z)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    b <- coef(fit); b[is.na(b)] <- 0
    res[[d]] <- vapply(names(targets), function(tn) {
      tg <- targets[[tn]]; Zt <- zf(tg$expr, f)
      if (use_age) Zt <- cbind(Zt, AGE = (tg$age - mu_a) / sd_a)
      cidx(tg$time, tg$event, drop(Zt %*% b))
    }, numeric(1))
  }
  null <- do.call(rbind, res); colnames(null) <- names(targets)
  summ <- do.call(rbind, lapply(names(targets), function(tn) {
    x <- null[, tn]; x <- x[is.finite(x)]
    r <- data.frame(cohort = tn, n_draw = length(x), null_mean = round(mean(x), 3), null_sd = round(sd(x), 3),
                    q2.5 = round(quantile(x, .025), 3), q97.5 = round(quantile(x, .975), 3),
                    stringsAsFactors = FALSE, row.names = NULL)
    if (!is.null(own_C) && !is.na(own_C[tn])) {
      r$own_C <- own_C[[tn]]
      r$p_emp_one_sided <- round((1 + sum(x >= own_C[[tn]])) / (1 + length(x)), 4)   # (1+#)/(1+n)
    }
    r
  }))
  cat("\n--- random same-size signatures (k =", k, if (use_age) "+ Age" else "", ") ---\n")
  print(summ, row.names = FALSE)
  invisible(list(draws = null, summary = summ))
}
