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
#   genes         optional data.table(gene_id, weight, ...) of per-gene weights
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
stopifnot(identical(rownames(X), split$barcode))
ranking <- fread(inp$ranking)

genes <- if (identical(mcfg$n_features, "all")) ranking$gene_id else head(ranking$gene_id, mcfg$n_features)
train <- split$set == "train"
X_train <- X[train, genes, drop = FALSE]
X_test <- X[!train, genes, drop = FALSE]
y_train <- split$class[train]
y_test <- split$class[!train]
rm(X)
message("Features: ", length(genes), " | train: ", nrow(X_train), " | test: ", nrow(X_test))

source(inp$method)
set.seed(SEED)
res <- fit_model(X_train, y_train, X_test, mcfg, cfg$classes, groups = split$patient[train])
saveRDS(res$fit, out$model)
pred <- factor(as.character(res$pred), levels = CLASS_LEVELS)
prob_pos <- as.numeric(res$prob_pos)
stopifnot(length(pred) == nrow(X_test), !anyNA(pred), length(prob_pos) == nrow(X_test))

# --- Predictions and metrics -----------------------------------------------------
predictions <- data.table(model = model_id, barcode = split$barcode[!train],
                          patient = split$patient[!train], actual = y_test,
                          predicted = pred, prob_positive = round(prob_pos, 4))
write_tsv(predictions, out$predictions)

cm <- table(predicted = pred, actual = y_test)
write_tsv(as.data.table(cm)[, model := model_id][, .(model, actual, predicted, N)], out$confusion)

# AUC = probability that a random positive scores higher than a random negative (ties count 1/2).
auc <- function(score, truth) {
  r <- rank(score)
  n_pos <- sum(truth == pos); n_neg <- sum(truth == neg)
  (sum(r[truth == pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
sens <- cm[pos, pos] / sum(cm[, pos])
spec <- cm[neg, neg] / sum(cm[, neg])
metrics <- data.table(
  model = model_id, method = mcfg$method, n_features = length(genes),
  n_train = nrow(X_train), n_test = nrow(X_test),
  accuracy = sum(diag(cm)) / sum(cm), sensitivity = sens, specificity = spec,
  balanced_accuracy = (sens + spec) / 2, auc = auc(prob_pos, y_test)
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
  gw <- merge(gw, split_gene_ids(gw$gene_id), by = "gene_id")
  gw <- gw[order(-abs(weight))]
  gw[, `:=`(model = model_id, abs_rank = .I)]
  if (isTRUE(res$signed)) gw[, direction := paste("Favors", ifelse(weight > 0, pos, neg))]
  setcolorder(gw, intersect(c("model", "gene_id", "symbol", "entrez", "weight", "abs_rank", "direction"),
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
