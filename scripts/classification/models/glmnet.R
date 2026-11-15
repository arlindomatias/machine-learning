# ============================================================
# Method glmnet: penalized logistic regression (LASSO for alpha = 1, elastic net for 0 < alpha < 1).
# Config: alpha, s (lambda.1se | lambda.min; default lambda.1se), nfolds (default 10),
#         balanced (weight each class equally; default false).
# The penalty is tuned by cross-validation on the training set, with folds that never split a
# patient. Gene weight = coefficient on standardized genes (glmnet's default standardization, so
# coefficients are comparable across genes); positive = towards the positive class.
# Sourced by 04_classification.R; see that script for the fit_model() contract.
# ============================================================
suppressPackageStartupMessages(library(glmnet))

# <figure name="method_glmnet">
fit_model <- function(X_train, y_train, X_test, mcfg, classes, groups) {
  s <- mcfg$s %||% "lambda.1se"
  w <- if (isTRUE(mcfg$balanced)) balanced_weights(y_train) else rep(1, length(y_train))
  folds <- group_folds(groups, y_train, mcfg$nfolds %||% 10)
  cv <- cv.glmnet(X_train, y_train, family = "binomial", alpha = mcfg$alpha, weights = w,
                  foldid = folds, type.measure = "deviance")
  print(cv)
  prob_pos <- as.numeric(predict(cv, X_test, s = s, type = "response"))  # P(second level)
  stopifnot(levels(y_train)[2] == classes$positive)

  # Coefficients on the standardized scale: raw coefficient x training-set sd.
  beta <- coef(cv, s = s)[-1, 1]
  sd_train <- apply(X_train, 2, sd)
  list(
    fit = cv,
    pred = ifelse(prob_pos >= 0.5, classes$positive, classes$negative),
    prob_pos = prob_pos,
    genes = data.table(feature_id = names(beta), weight = beta * sd_train[names(beta)], coefficient = beta),
    weight_label = "Standardized coefficient",
    signed = TRUE
  )
}
# </figure>
