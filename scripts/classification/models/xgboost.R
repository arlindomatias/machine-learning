# ============================================================
# Method xgboost: gradient-boosted trees (binary:logistic).
# Config: nrounds (default 200), eta (default 0.1), max_depth (default 3), subsample (default 0.8),
#         colsample_bytree (default 0.8), balanced (scale_pos_weight = negatives / positives;
#         default true).
# Score = predicted probability of the positive class. Gene weight = gain (share of the loss
# reduction from splits on the gene); genes never used in a split are dropped.
# Sourced by 04_classification.R; see that script for the fit_model() contract.
# ============================================================
suppressPackageStartupMessages(library(xgboost))

# <figure name="method_xgboost">
fit_model <- function(X_train, y_train, X_test, mcfg, classes, groups) {
  label <- as.integer(y_train == classes$positive)
  params <- list(objective = "binary:logistic", eval_metric = "logloss",
                 eta = mcfg$eta %||% 0.1, max_depth = mcfg$max_depth %||% 3,
                 subsample = mcfg$subsample %||% 0.8, colsample_bytree = mcfg$colsample_bytree %||% 0.8,
                 nthread = 1, seed = SEED)
  if (isTRUE(mcfg$balanced %||% TRUE)) params$scale_pos_weight <- sum(label == 0) / sum(label == 1)
  dtrain <- xgb.DMatrix(X_train, label = label)
  fit <- xgb.train(params = params, data = dtrain, nrounds = mcfg$nrounds %||% 200)
  prob_pos <- predict(fit, xgb.DMatrix(X_test))
  imp <- as.data.table(xgb.importance(model = fit))
  list(
    fit = fit,
    pred = ifelse(prob_pos >= 0.5, classes$positive, classes$negative),
    prob_pos = prob_pos,
    genes = imp[, .(feature_id = Feature, weight = Gain, cover = Cover, frequency = Frequency)],
    weight_label = "Gain",
    signed = FALSE
  )
}
# </figure>
