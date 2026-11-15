# ============================================================
# Method dlda: diagonal linear discriminant analysis (Dudoit et al., 2002): Gaussian classes with a
# shared, diagonal covariance (genes treated as independent, pooled within-class variances).
# Config: prior (balanced | empirical; default balanced).
# Score = posterior probability of the positive class. Gene weight = standardized mean difference
# (positive minus negative class mean, over the pooled sd); every gene contributes.
# Sourced by 04_classification.R; see that script for the fit_model() contract.
# ============================================================

# <figure name="method_dlda">
fit_model <- function(X_train, y_train, X_test, mcfg, classes, groups) {
  lv <- levels(y_train)
  means <- sapply(lv, function(k) colMeans(X_train[y_train == k, , drop = FALSE]))
  resid <- X_train - t(means[, as.character(y_train)])
  pooled_var <- colSums(resid^2) / (nrow(X_train) - length(lv))
  pooled_var[pooled_var == 0] <- min(pooled_var[pooled_var > 0])
  prior <- if ((mcfg$prior %||% "balanced") == "balanced") setNames(rep(1 / length(lv), length(lv)), lv)
           else c(table(y_train) / length(y_train))

  # Log posterior up to a constant: -sum_j (x_j - mu_kj)^2 / (2 var_j) + log prior_k.
  log_post <- sapply(lv, function(k) {
    -colSums((t(X_test) - means[, k])^2 / (2 * pooled_var)) + log(prior[[k]])
  })
  post <- exp(log_post - apply(log_post, 1, max))
  post <- post / rowSums(post)
  list(
    fit = list(method = "dlda", means = means, pooled_var = pooled_var, prior = prior),
    pred = lv[max.col(post, ties.method = "first")],
    prob_pos = post[, classes$positive],
    genes = data.table(gene_id = rownames(means),
                       weight = (means[, classes$positive] - means[, classes$negative]) / sqrt(pooled_var)),
    weight_label = "Standardized mean difference",
    signed = TRUE
  )
}
# </figure>
