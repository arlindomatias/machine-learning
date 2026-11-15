# ============================================================
# Method random_forest: randomForest, optionally class-balanced.
# Config: ntree, balanced (each tree is grown on as many positives as negatives, via strata + sampsize).
# Gene weight = mean decrease in accuracy (permutation importance).
# Sourced by 04_classification.R; see that script for the fit_model() contract.
# ============================================================
suppressPackageStartupMessages(library(randomForest))

# <figure name="method_random_forest">
fit_model <- function(X_train, y_train, X_test, mcfg, classes, groups) {
  args <- list(x = X_train, y = y_train, ntree = mcfg$ntree, importance = TRUE)
  if (isTRUE(mcfg$balanced)) {
    n_min <- min(table(y_train))
    args <- c(args, list(strata = y_train, sampsize = c(n_min, n_min)))
  }
  fit <- do.call(randomForest, args)
  print(fit)
  imp <- importance(fit)
  list(
    fit = fit,
    pred = predict(fit, X_test),
    prob_pos = predict(fit, X_test, type = "prob")[, classes$positive],
    genes = data.table(gene_id = rownames(imp),
                       weight = imp[, "MeanDecreaseAccuracy"],
                       mean_decrease_gini = imp[, "MeanDecreaseGini"]),
    weight_label = "Mean decrease in accuracy",
    signed = FALSE
  )
}
# </figure>
