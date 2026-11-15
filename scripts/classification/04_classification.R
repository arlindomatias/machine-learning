# ============================================================
# SCRIPT 04: train and evaluate ONE model (Snakemake wildcard `model`).
# The model id is a key of `models:` in config.yaml and prefixes every output file, so each
# model keeps its own predictions, metrics, confusion matrix, figures and fitted object.
#
# Method-specific code lives in models/<method>.R (Snakemake input `method`). Each of those files
# defines
#   fit_model(X_train, y_train, X_test, mcfg, classes, groups)
# where X_* are samples x genes matrices, y_train is a factor with levels CLASS_LEVELS, mcfg is the
# model's config entry, classes = cfg$classes and groups holds the training samples' patient ids
# (for internal cross-validation folds that never split a patient). It returns a list with
#   fit           fitted object, saved to data/models/<id>.rds
#   pred          predicted class of each test sample
#   prob_pos      score of each test sample for the positive class (higher = more positive)
#   genes         optional data.table(feature_id, weight, ...) of per-gene weights
#   weight_label  axis label for the weight
#   signed        TRUE if the sign of the weight gives the direction (positive = towards the
#                 positive class)
# Everything that must be identical across models (data, gene selection, metrics, figures) is here.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))

model_id <- snakemake@wildcards$model
mcfg <- cfg$models[[model_id]]
if (is.null(mcfg)) stop("Model not defined in config.yaml: ", model_id)
message("Model: ", model_id, " (", mcfg$method, ")")
pos <- cfg$classes$positive
neg <- cfg$classes$negative

X <- readRDS(inp$expression)
split <- fread(inp$split)
split[, class := factor(class, levels = CLASS_LEVELS)]
stopifnot(identical(rownames(X), split$sample_id))
ranking <- fread(inp$ranking)

genes <- if (identical(mcfg$n_features, "all")) ranking$feature_id else head(ranking$feature_id, mcfg$n_features)
train <- split$set == "train"
X_train <- X[train, genes, drop = FALSE]
X_test <- X[!train, genes, drop = FALSE]
y_train <- split$class[train]
y_test <- split$class[!train]
rm(X)
message("Features: ", length(genes), " | train: ", nrow(X_train), " | test: ", nrow(X_test))

source(inp$method)
set.seed(SEED)
res <- fit_model(X_train, y_train, X_test, mcfg, cfg$classes, groups = split$group[train])
saveRDS(res$fit, out$model)
pred <- factor(as.character(res$pred), levels = CLASS_LEVELS)
prob_pos <- as.numeric(res$prob_pos)
stopifnot(length(pred) == nrow(X_test), !anyNA(pred), length(prob_pos) == nrow(X_test))

# --- Predictions and metrics -----------------------------------------------------
predictions <- data.table(model = model_id, sample_id = split$sample_id[!train],
                          group = split$group[!train], actual = y_test,
                          predicted = pred, prob_positive = prob_pos)  # full precision: rounding creates ties that change AUC
write_tsv(predictions, out$predictions)

cm <- table(predicted = pred, actual = y_test)
write_tsv(as.data.table(cm)[, model := model_id][, .(model, actual, predicted, N)], out$confusion)

# Metrics (classification_metrics() in 00_common.R) plus Wilson intervals for sensitivity and
# specificity. Bootstrap intervals for every metric are computed in 05_compare_models.R.
ev <- cfg$evaluation
ci_s <- wilson_ci(cm[pos, pos], sum(cm[, pos]), ev$ci_level)
ci_e <- wilson_ci(cm[neg, neg], sum(cm[, neg]), ev$ci_level)
metrics <- cbind(
  data.table(model = model_id, method = mcfg$method, n_features = length(genes),
             n_train = nrow(X_train), n_test = nrow(X_test)),
  classification_metrics(y_test, pred, prob_pos, eps = ev$log_loss_eps),
  data.table(sensitivity_wilson_lower = ci_s[1], sensitivity_wilson_upper = ci_s[2],
             specificity_wilson_lower = ci_e[1], specificity_wilson_upper = ci_e[2])
)
print(cm)
print(metrics)
write_tsv(metrics, out$metrics)

# <figure name="confusion_matrix">
cm_dt <- as.data.table(cm)
cm_dt[, share := N / sum(N), by = actual]
p <- ggplot(cm_dt, aes(actual, predicted, fill = share)) +
  geom_tile(color = "white") +
  geom_text(aes(label = N, color = share > 0.5), size = 6, show.legend = FALSE) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 1), name = "Share of\nactual class") +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  labs(title = model_id, x = "Actual class", y = "Predicted class",
       subtitle = sprintf("Balanced accuracy %.3f | AUC %.3f", metrics$balanced_accuracy, metrics$auc)) +
  coord_equal()
save_plot(out$confusion_plot, p, width = 5, height = 4.5)
# </figure>

# --- Per-gene weights (methods that have them) ----------------------------------------
if (!is.null(out$genes)) {
  if (is.null(res$genes)) stop("Method ", mcfg$method, " returned no per-gene weights")
  gw <- as.data.table(res$genes)[weight != 0]
  gw <- merge(gw, split_gene_ids(gw$feature_id), by = "feature_id")
  gw <- gw[order(-abs(weight))]
  gw[, `:=`(model = model_id, abs_rank = .I)]
  if (isTRUE(res$signed)) gw[, direction := paste("Favors", ifelse(weight > 0, pos, neg))]
  setcolorder(gw, intersect(c("model", "feature_id", "symbol", "entrez", "weight", "abs_rank", "direction"),
                            names(gw)))
  message("Genes with non-zero weight: ", nrow(gw))
  write_tsv(gw, out$genes)

  # <figure name="genes">
  top20 <- head(gw, 20)
  p <- ggplot(top20, aes(weight, reorder(symbol, abs(weight)))) +
    labs(title = sprintf("%s: top %d genes", model_id, nrow(top20)), x = res$weight_label, y = NULL)
  p <- if (isTRUE(res$signed)) {
    p + geom_col(aes(fill = direction)) + geom_vline(xintercept = 0, color = "grey40") +
      scale_fill_manual(values = setNames(CLASS_COLORS[c(pos, neg)], paste("Favors", c(pos, neg))),
                        name = NULL) +
      theme(legend.position = "bottom")
  } else {
    p + geom_col(fill = "#4c72b0")
  }
  save_plot(out$genes_plot, p, width = 5.5, height = 5.5)
  # </figure>
}
