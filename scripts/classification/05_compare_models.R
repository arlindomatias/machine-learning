# ============================================================
# SCRIPT 05: comparison of every model in config.yaml on the shared test set.
#   - metrics table (from the per-model metrics files)
#   - bootstrap confidence intervals: test-set patients are resampled with replacement and the same
#     resamples are used for every model (paired)
#   - pairwise tests: exact McNemar on the errors, DeLong on the AUCs; Holm-adjusted across pairs
#   - calibration: observed share of positives per bin of predicted probability
#   - misclassified samples: which samples each model gets wrong
# Metric functions are in 00_common.R.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))

pos <- cfg$classes$positive
neg <- cfg$classes$negative
ev <- cfg$evaluation
model_order <- names(cfg$models)

metrics <- rbindlist(lapply(inp$metrics, fread))
metrics[, model := factor(model, levels = model_order)]
setorder(metrics, model)
print(metrics)
write_tsv(metrics, out$comparison)

predictions <- rbindlist(lapply(inp$predictions, fread))
predictions[, model := factor(model, levels = model_order)]
setorder(predictions, model, barcode)
# Wide views: one row per test sample, one column per model.
samples <- unique(predictions[, .(barcode, patient, actual)])
setorder(samples, barcode)
prob <- dcast(predictions, barcode ~ model, value.var = "prob_positive")
called <- dcast(predictions, barcode ~ model, value.var = "predicted")
stopifnot(identical(prob$barcode, samples$barcode), identical(called$barcode, samples$barcode))
is_pos <- samples$actual == pos

# The tests and intervals below recompute metrics from the prediction files; they must match 04's.
for (m in model_order) {
  again <- classification_metrics(samples$actual, called[[m]], prob[[m]], ev$log_loss_eps)
  stored <- metrics[model == m, names(again), with = FALSE]
  if (!isTRUE(all.equal(as.list(again), as.list(stored), tolerance = 1e-6, check.attributes = FALSE)))
    stop("Metrics recomputed from the predictions of ", m, " differ from its metrics file")
}

# --- Bootstrap confidence intervals -------------------------------------------------------
set.seed(SEED)
rows_by_patient <- split(seq_len(nrow(samples)), samples$patient)
boot <- rbindlist(lapply(seq_len(ev$bootstrap), function(b) {
  idx <- unlist(rows_by_patient[sample(length(rows_by_patient), replace = TRUE)], use.names = FALSE)
  rbindlist(lapply(model_order, function(m) {
    cbind(data.table(model = m, b = b),
          classification_metrics(samples$actual[idx], called[[m]][idx], prob[[m]][idx], ev$log_loss_eps))
  }))
}))
metric_cols <- setdiff(names(boot), c("model", "b"))
alpha <- (1 - ev$ci_level) / 2
ci <- melt(boot, id.vars = c("model", "b"), measure.vars = metric_cols, variable.name = "metric")[
  , .(lower = quantile(value, alpha, na.rm = TRUE), upper = quantile(value, 1 - alpha, na.rm = TRUE)),
  by = .(model, metric)]
point <- melt(metrics, id.vars = "model", measure.vars = metric_cols, variable.name = "metric",
              value.name = "estimate")
ci <- merge(point, ci, by = c("model", "metric"))
ci[, `:=`(model = factor(model, levels = model_order), metric = factor(metric, levels = metric_cols))]
setorder(ci, metric, model)
write_tsv(ci, out$bootstrap)

# --- Pairwise tests -----------------------------------------------------------------------
# DeLong et al. (1988): AUC variance and covariance from the structural components
# V10 (per positive) and V01 (per negative), computed from ranks.
delong_components <- function(score) {
  x <- score[is_pos]; y <- score[!is_pos]
  n1 <- length(x); n0 <- length(y)
  r_all <- rank(c(x, y))
  list(v10 = (r_all[seq_len(n1)] - rank(x)) / n0,
       v01 = 1 - (r_all[n1 + seq_len(n0)] - rank(y)) / n1)
}
comp <- lapply(setNames(model_order, model_order), function(m) delong_components(prob[[m]]))
pairs <- combn(model_order, 2, simplify = FALSE)
pairwise <- rbindlist(lapply(pairs, function(pr) {
  a <- pr[1]; b <- pr[2]
  right_a <- called[[a]] == samples$actual
  right_b <- called[[b]] == samples$actual
  only_a <- sum(right_a & !right_b); only_b <- sum(!right_a & right_b)
  p_mcnemar <- if (only_a + only_b > 0) binom.test(only_a, only_a + only_b)$p.value else 1
  ca <- comp[[a]]; cb <- comp[[b]]
  auc_a <- mean(ca$v10); auc_b <- mean(cb$v10)
  var_diff <- var(ca$v10 - cb$v10) / length(ca$v10) + var(ca$v01 - cb$v01) / length(ca$v01)
  z <- if (var_diff > 0) (auc_a - auc_b) / sqrt(var_diff) else 0
  data.table(model_a = a, model_b = b, errors_a = sum(!right_a), errors_b = sum(!right_b),
             only_a_correct = only_a, only_b_correct = only_b, mcnemar_p = p_mcnemar,
             auc_a = auc_a, auc_b = auc_b, delong_z = z, delong_p = 2 * pnorm(-abs(z)))
}))
pairwise[, `:=`(mcnemar_p_holm = p.adjust(mcnemar_p, "holm"), delong_p_holm = p.adjust(delong_p, "holm"))]
setcolorder(pairwise, c("model_a", "model_b", "errors_a", "errors_b", "only_a_correct", "only_b_correct",
                        "mcnemar_p", "mcnemar_p_holm", "auc_a", "auc_b", "delong_z", "delong_p", "delong_p_holm"))
write_tsv(pairwise, out$pairwise)
message("Pairs with Holm-adjusted p < 0.05: McNemar ", sum(pairwise$mcnemar_p_holm < 0.05),
        ", DeLong ", sum(pairwise$delong_p_holm < 0.05), " of ", nrow(pairwise))

# --- Misclassified samples ----------------------------------------------------------------
wrong <- as.matrix(called[, ..model_order]) != samples$actual
errors <- cbind(samples, n_models_wrong = rowSums(wrong),
                models_wrong = apply(wrong, 1, function(w) paste(model_order[w], collapse = ", ")),
                mean_prob_positive = round(rowMeans(as.matrix(prob[, ..model_order])), 4))[n_models_wrong > 0]
setorder(errors, -n_models_wrong, barcode)
write_tsv(errors, out$errors)

# <figure name="model_comparison">
show <- c(sensitivity = "Sensitivity", specificity = "Specificity", balanced_accuracy = "Balanced accuracy",
          mcc = "MCC", auc = "ROC AUC", pr_auc_negative = paste0("PR AUC (", neg, ")"),
          npv = "NPV", brier = "Brier score (lower is better)")
d <- ci[metric %in% names(show)]
d[, metric := factor(show[as.character(metric)], levels = show)]
d[, model := factor(model, levels = rev(model_order))]
p <- ggplot(d, aes(estimate, model)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.3, orientation = "y", color = "grey45") +
  geom_point(color = "#08519c", size = 1.8) +
  facet_wrap(~metric, ncol = 4, scales = "free_x") +
  labs(x = sprintf("Test-set value (%.0f%% patient bootstrap interval)", 100 * ev$ci_level), y = NULL,
       title = "Model comparison") +
  theme(strip.text = element_text(size = 9))
save_plot(out$comparison_plot, p, width = 11, height = 6)
# </figure>

# <figure name="roc_curves">
roc <- predictions[, {
  thr <- sort(unique(c(Inf, prob_positive)), decreasing = TRUE)
  .(fpr = vapply(thr, function(t) mean(prob_positive[actual != pos] >= t), numeric(1)),
    tpr = vapply(thr, function(t) mean(prob_positive[actual == pos] >= t), numeric(1)))
}, by = model]
labels <- metrics[, .(model, label = sprintf("%s (AUC %.3f)", model, auc))]
roc <- merge(roc, labels, by = "model")
p <- ggplot(roc, aes(fpr, tpr, color = label)) +
  geom_abline(linetype = "dashed", color = "grey60") +
  geom_step(direction = "vh", linewidth = 0.7) +
  coord_equal() +
  labs(x = "False positive rate (1 - specificity)", y = "True positive rate (sensitivity)",
       color = NULL, title = "ROC curves, test set") +
  theme(legend.position = "bottom", legend.direction = "vertical")
save_plot(out$roc_plot, p, width = 5.5, height = 6.5)
# </figure>

# <figure name="calibration">
breaks <- seq(0, 1, length.out = ev$calibration_bins + 1)
cal <- predictions[, .(bin = cut(prob_positive, breaks, include.lowest = TRUE), prob_positive,
                       is_pos = actual == pos), by = model][
  , .(predicted = mean(prob_positive), observed = mean(is_pos), n = .N), by = .(model, bin)]
p <- ggplot(cal, aes(predicted, observed)) +
  geom_abline(linetype = "dashed", color = "grey60") +
  geom_point(aes(size = n), color = "#08519c", alpha = 0.8) +
  scale_size_area(max_size = 5, breaks = c(1, 10, 100), name = "Samples") +
  facet_wrap(~model, ncol = 4) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = sprintf("Mean predicted P(%s) per bin", pos), y = sprintf("Observed share of %s", pos),
       title = "Calibration, test set")
save_plot(out$calibration_plot, p, width = 9, height = 2.6 * ceiling(length(model_order) / 4) + 0.8)
# </figure>

# <figure name="misclassified">
if (nrow(errors) > 0) {
  tiles <- predictions[barcode %in% errors$barcode]
  tiles[, `:=`(correct = predicted == actual,
               sample = factor(sprintf("%s (%s)", barcode, actual),
                               levels = rev(errors[, sprintf("%s (%s)", barcode, actual)])),
               model = factor(model, levels = model_order))]
  p <- ggplot(tiles, aes(model, sample, fill = correct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", prob_positive)), size = 2.6) +
    scale_fill_manual(values = c(`TRUE` = "#e5f0e5", `FALSE` = "#f4a6a6"),
                      labels = c(`TRUE` = "Correct", `FALSE` = "Wrong"), name = NULL) +
    labs(x = NULL, y = NULL, title = "Test samples misclassified by at least one model",
         subtitle = sprintf("Cell label: predicted P(%s)", pos)) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")
  save_plot(out$errors_plot, p, width = 9, height = 1.6 + 0.28 * nrow(errors))
} else {
  save_plot(out$errors_plot, ggplot() + labs(title = "No misclassified test samples") + theme_void())
}
# </figure>
