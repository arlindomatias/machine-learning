# ============================================================
# Method knn: class::knn on genes standardized with training-set mean and sd.
# Config: k. Score = share of neighbour votes for the positive class.
# Sourced by 04_classification.R; see that script for the fit_model() contract.
# ============================================================
suppressPackageStartupMessages(library(class))

# <figure name="method_knn">
fit_model <- function(X_train, y_train, X_test, mcfg, classes, groups) {
  mu <- colMeans(X_train)
  sigma <- apply(X_train, 2, sd)
  sigma[sigma == 0] <- 1
  pred <- knn(scale(X_train, mu, sigma), scale(X_test, mu, sigma), y_train, k = mcfg$k, prob = TRUE)
  vote <- attr(pred, "prob")                      # share of votes for the predicted class
  list(
    fit = list(method = "knn", k = mcfg$k, genes = colnames(X_train), center = mu, scale = sigma),
    pred = pred,
    prob_pos = ifelse(pred == classes$positive, vote, 1 - vote)
  )
}
# </figure>
