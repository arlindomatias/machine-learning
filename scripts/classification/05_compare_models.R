# ============================================================
# SCRIPT 05: side-by-side comparison of every model in config.yaml (metrics table, bar chart,
# ROC curves on the shared test set).
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))

pos <- cfg$classes$positive
metrics <- rbindlist(lapply(inp$metrics, fread))
model_order <- names(cfg$models)
metrics[, model := factor(model, levels = model_order)]
setorder(metrics, model)
print(metrics)
write_tsv(metrics, out$comparison)

predictions <- rbindlist(lapply(inp$predictions, fread))
predictions[, model := factor(model, levels = model_order)]

# <figure name="model_comparison">
long <- melt(metrics, id.vars = "model",
             measure.vars = c("sensitivity", "specificity", "balanced_accuracy", "auc"),
             variable.name = "metric")
lower <- floor(min(long$value) * 20) / 20
p <- ggplot(long, aes(model, value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  coord_cartesian(ylim = c(lower, 1)) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = NULL, y = "Test-set value", title = "Model comparison") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
save_plot(out$comparison_plot, p, width = 7, height = 4.5)
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
