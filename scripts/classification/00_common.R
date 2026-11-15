# ============================================================
# SCRIPT 00: shared setup (sourced by every step)
# Dataset-agnostic: every dataset-specific value comes from the project's config.yaml
# (`snakemake@config`). Scripts run through the Snakemake `script:` directive.
# ============================================================
if (!exists("snakemake")) {
  stop("Run this script through Snakemake (`snakemake --use-conda --cores 2`).")
}

# R messages and output go to the rule log.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

cfg <- snakemake@config
out <- snakemake@output
inp <- snakemake@input

for (f in unlist(out)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

options(timeout = 600)
SEED <- cfg$seed %||% 123
set.seed(SEED)
if (!interactive()) pdf(NULL)

CLASS_LEVELS <- c(cfg$classes$negative, cfg$classes$positive)
CLASS_COLORS <- unlist(cfg$classes$colors)

theme_set(theme_bw(base_size = 11))

save_plot <- function(path, plot, ...) {
  tryCatch(ggsave(path, plot = plot, dpi = 300, bg = "white", ...),
           error = function(e) message("Plot ", path, " failed: ", conditionMessage(e)))
}

write_tsv <- function(x, path) fwrite(x, path, sep = "\t")

# --- Input contract ---------------------------------------------------------------
# The prepare step (01_prepare.R + inputs/<format>.R) is the only dataset-specific step. It writes
#   data/expression.rds  numeric matrix, samples x features, on the analysis scale (e.g. log2)
#   data/samples.csv     sample_id, group, class [, set] [, any other columns]
#   data/features.csv    feature_id, symbol, entrez (NA if unknown)
# and every later step reads only these files. See README.md in this folder.

# Stops with a message listing every violation of the contract.
validate_inputs <- function(X, samples, features) {
  problems <- character()
  check <- function(ok, msg) if (!isTRUE(ok)) problems <<- c(problems, msg)
  check(length(CLASS_LEVELS) == 2 && !anyDuplicated(CLASS_LEVELS),
        "config classes must define two different classes (negative, positive)")
  check(is.matrix(X) && is.numeric(X), "expression must be a numeric matrix")
  check(all(is.finite(X)), "expression contains NA, NaN or Inf")
  check(!anyDuplicated(rownames(X)) && !anyDuplicated(colnames(X)), "duplicated sample or feature ids")
  check(all(c("sample_id", "group", "class") %in% names(samples)), "samples needs sample_id, group, class")
  check(identical(rownames(X), as.character(samples$sample_id)), "matrix rows != samples$sample_id (same order)")
  check(!anyNA(samples$group) && all(samples$group != ""), "samples$group has missing values")
  check(all(as.character(samples$class) %in% CLASS_LEVELS),
        paste("samples$class must be one of:", paste(CLASS_LEVELS, collapse = ", ")))
  check(all(CLASS_LEVELS %in% as.character(samples$class)), "both classes must be present")
  if ("set" %in% names(samples)) {
    check(all(samples$set %in% c("train", "test")), "samples$set must be train or test")
    check(samples[, uniqueN(set), by = group][, all(V1 == 1)], "a group has samples in both train and test")
  }
  check(all(c("feature_id", "symbol", "entrez") %in% names(features)), "features needs feature_id, symbol, entrez")
  check(identical(colnames(X), as.character(features$feature_id)), "matrix columns != features$feature_id (same order)")
  check(!anyNA(features$symbol), "features$symbol has missing values")
  if (length(problems)) stop("Input contract violated:\n  - ", paste(problems, collapse = "\n  - "))
  invisible(TRUE)
}

read_expression <- function(path) readRDS(path)

read_samples <- function(path) {
  s <- fread(path, colClasses = list(character = c("sample_id", "group")))
  s[, class := factor(class, levels = CLASS_LEVELS)]
  s
}

read_features <- function(path) {
  f <- fread(path, colClasses = list(character = c("feature_id", "symbol", "entrez")))
  f[entrez == "", entrez := NA_character_]
  f
}

# Adds symbol and entrez to a table with a feature_id column.
annotate_features <- function(dt, features, ...) {
  merge(dt, features[, .(feature_id, symbol, entrez)], by = "feature_id", ...)
}

# Name of the grouping unit in text and plots (e.g. "patient"); config `group_name`.
GROUP_NAME <- cfg$group_name %||% "group"

# --- Feature selection ----------------------------------------------------------

# Between-group / within-group sum-of-squares ratio per gene (Dudoit et al., 2002).
# Vectorized over genes; any number of classes. Constant genes are dropped.
bss_wss <- function(X, y) {
  y <- as.factor(y)
  overall <- colMeans(X)
  bss <- wss <- numeric(ncol(X))
  for (lv in levels(y)) {
    Xk <- X[y == lv, , drop = FALSE]
    mk <- colMeans(Xk)
    bss <- bss + nrow(Xk) * (mk - overall)^2
    wss <- wss + colSums(sweep(Xk, 2, mk)^2)
  }
  ratio <- bss / wss
  res <- data.table(feature_id = colnames(X), bss_wss = ratio)[is.finite(bss_wss)]
  setorder(res, -bss_wss)
  res[, rank := .I]
  res
}

# --- Cross-validation folds -----------------------------------------------------------

# Fold id per sample for internal cross-validation (e.g. tuning a penalty): whole groups (patients)
# go to one fold, and groups with a negative-class sample are spread evenly across folds.
group_folds <- function(groups, y, k) {
  g <- data.table(group = groups, neg = y == cfg$classes$negative)[, .(neg = any(neg)), by = group]
  g[, fold := sample(rep_len(seq_len(k), .N)), by = neg]
  g$fold[match(groups, g$group)]
}

# Observation weights giving both classes the same total weight.
balanced_weights <- function(y) {
  n <- table(y)
  as.numeric(length(y) / (length(n) * n[as.character(y)]))
}

# --- Performance metrics ------------------------------------------------------------
# Shared by 04_classification.R (per model) and 05_compare_models.R (bootstrap, pairwise tests).
# `actual` and `predicted` are factors or characters with class names; `prob_pos` is the score for
# the positive class (higher = more positive).

# AUC = probability that a random positive scores higher than a random negative (ties count 1/2).
auc_score <- function(score, is_pos) {
  n_pos <- sum(is_pos); n_neg <- sum(!is_pos)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[is_pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

# Average precision (area under the precision-recall curve, step-wise; tied scores form one step).
average_precision <- function(score, is_pos) {
  if (!any(is_pos)) return(NA_real_)
  o <- order(score, decreasing = TRUE)
  s <- score[o]
  tp <- cumsum(is_pos[o])
  n_called <- seq_along(s)
  last <- c(s[-1] != s[-length(s)], TRUE)          # last sample of each group of tied scores
  tp <- tp[last]; n_called <- n_called[last]
  sum(diff(c(0, tp / sum(is_pos))) * tp / n_called)
}

# Wilson score interval for a proportion x / n.
wilson_ci <- function(x, n, level = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - level) / 2)
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(max(0, centre - half), min(1, centre + half))
}

# One row of metrics. `eps` clips probabilities for the log loss (a hard 0/1 score that is wrong
# would otherwise give an infinite loss).
classification_metrics <- function(actual, predicted, prob_pos, eps = 1e-3) {
  pos <- cfg$classes$positive
  is_pos <- as.character(actual) == pos
  called_pos <- as.character(predicted) == pos
  tp <- sum(is_pos & called_pos); tn <- sum(!is_pos & !called_pos)
  fp <- sum(!is_pos & called_pos); fn <- sum(is_pos & !called_pos)
  sens <- tp / (tp + fn); spec <- tn / (tn + fp)
  mcc_den <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  p <- pmin(pmax(prob_pos, eps), 1 - eps)
  data.table(
    accuracy = (tp + tn) / length(actual),
    sensitivity = sens, specificity = spec, balanced_accuracy = (sens + spec) / 2,
    ppv = if (tp + fp > 0) tp / (tp + fp) else NA_real_,
    npv = if (tn + fn > 0) tn / (tn + fn) else NA_real_,
    mcc = if (mcc_den > 0) (as.numeric(tp) * tn - as.numeric(fp) * fn) / mcc_den else 0,
    auc = auc_score(prob_pos, is_pos),
    pr_auc_negative = average_precision(-prob_pos, !is_pos),
    brier = mean((prob_pos - is_pos)^2),
    log_loss = -mean(ifelse(is_pos, log(p), log(1 - p)))
  )
}
