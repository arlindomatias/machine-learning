# ============================================================
# SCRIPT 04: train and evaluate ONE model (Snakemake wildcard `model`).
# The model id is a key of `models:` in config.yaml and prefixes every output file, so each
# model keeps its own predictions, metrics, confusion matrix, figures and fitted object.
#   method knn:           class::knn on genes standardized with training-set statistics
#   method random_forest: randomForest, optionally class-balanced (strata + sampsize)
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

set.seed(SEED)
if (mcfg$method == "knn") {
  suppressPackageStartupMessages(library(class))
  mu <- colMeans(X_train)
  sigma <- apply(X_train, 2, sd)
  sigma[sigma == 0] <- 1
  pred <- knn(scale(X_train, mu, sigma), scale(X_test, mu, sigma), y_train, k = mcfg$k, prob = TRUE)
  vote <- attr(pred, "prob")                      # share of votes for the predicted class
  prob_pos <- ifelse(pred == pos, vote, 1 - vote)
  fit <- list(method = "knn", k = mcfg$k, genes = genes, center = mu, scale = sigma)
} else if (mcfg$method == "random_forest") {
  suppressPackageStartupMessages(library(randomForest))
  args <- list(x = X_train, y = y_train, ntree = mcfg$ntree, importance = TRUE)
  if (isTRUE(mcfg$balanced)) {
    n_min <- min(table(y_train))
    args <- c(args, list(strata = y_train, sampsize = c(n_min, n_min)))
  }
  fit <- do.call(randomForest, args)
  print(fit)
  pred <- predict(fit, X_test)
  prob_pos <- predict(fit, X_test, type = "prob")[, pos]
} else {
  stop("Unknown method: ", mcfg$method)
}
saveRDS(fit, out$model)

# --- Predictions and metrics -----------------------------------------------------
pred <- factor(as.character(pred), levels = CLASS_LEVELS)
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

# --- Variable importance (random forest only) -----------------------------------------
if (!is.null(out$importance)) {
  imp <- importance(fit)
  imp <- data.table(model = model_id, gene_id = rownames(imp),
                    mean_decrease_accuracy = imp[, "MeanDecreaseAccuracy"],
                    mean_decrease_gini = imp[, "MeanDecreaseGini"])
  imp <- merge(imp, split_gene_ids(imp$gene_id), by = "gene_id")
  setorder(imp, -mean_decrease_accuracy)
  setcolorder(imp, c("model", "gene_id", "symbol", "entrez"))
  write_tsv(imp, out$importance)

  # <figure name="importance">
  top20 <- head(imp, 20)
  p <- ggplot(top20, aes(mean_decrease_accuracy, reorder(symbol, mean_decrease_accuracy))) +
    geom_col(fill = "#4c72b0") +
    labs(title = sprintf("%s: top 20 genes", model_id), x = "Mean decrease in accuracy", y = NULL)
  save_plot(out$importance_plot, p, width = 5.5, height = 5)
  # </figure>
}
