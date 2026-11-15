# ============================================================
# Method pamr: nearest shrunken centroids (Tibshirani et al., 2002).
# Config: balanced (default true), nfolds (default 10).
#
# balanced: every class is replicated up to the size of the largest one before training.
#   pamr shrinks each class centroid towards the overall centroid, which is weighted by class size.
#   With ~10 tumors per normal the overall centroid sits on the tumor centroid, so shrinkage moves
#   only the normal centroid and the boundary drifts into the tumors (on this cohort: ~10% balanced
#   CV error at 100 genes, every sample one class below ~40 genes). Balancing puts the overall
#   centroid halfway between the classes. Priors are equal in both modes.
# The shrinkage threshold is tuned by cross-validation on the training set (folds never split a
# patient; replicates of a sample share its fold): the largest threshold with the minimum
# class-balanced CV error, i.e. the fewest genes at the best error.
# Score = posterior probability of the positive class. Gene weight = shrunken standardized
# centroid difference (positive class minus overall centroid); genes shrunk to zero are dropped.
# Sourced by 04_classification.R; see that script for the fit_model() contract.
# ============================================================
suppressPackageStartupMessages(library(pamr))

# <figure name="method_pamr">
fit_model <- function(X_train, y_train, X_test, mcfg, classes, groups) {
  idx <- seq_along(y_train)
  if (isTRUE(mcfg$balanced %||% TRUE)) {
    n <- table(y_train)
    idx <- unlist(lapply(names(n), function(k) rep(which(y_train == k), round(max(n) / n[[k]]))))
  }
  y <- y_train[idx]
  data <- list(x = t(X_train[idx, , drop = FALSE]), y = y,
               geneid = colnames(X_train), genenames = colnames(X_train))
  fit <- pamr.train(data, prior = rep(1 / nlevels(y), nlevels(y)))
  folds <- patient_folds(groups[idx], y, mcfg$nfolds %||% 10)
  cv <- pamr.cv(fit, data, folds = split(seq_along(folds), folds))
  bal_err <- apply(cv$yhat, 2, function(p) mean(tapply(as.character(p) != as.character(y), y, mean)))
  threshold <- max(cv$threshold[bal_err == min(bal_err)])
  message(sprintf("Threshold %.3f, balanced CV error %.4f, %d genes",
                  threshold, min(bal_err), cv$size[match(threshold, cv$threshold)]))

  post <- pamr.predict(fit, t(X_test), threshold, type = "posterior")
  # Shrunken standardized centroid difference, as in pamr:::nsc: (class centroid - overall) / sd
  # (sd already includes the offset s0), scaled per class, then soft-thresholded.
  d <- (fit$centroids - fit$centroid.overall) / fit$sd
  d <- sweep(d, 2, fit$threshold.scale * fit$se.scale, "/")
  shrunk <- sign(d) * pmax(abs(d) - threshold, 0)
  list(
    fit = list(fit = fit, cv = data.table(threshold = cv$threshold, balanced_error = bal_err, size = cv$size),
               threshold = threshold),
    pred = pamr.predict(fit, t(X_test), threshold, type = "class"),
    prob_pos = post[, classes$positive],
    genes = data.table(gene_id = rownames(shrunk), weight = shrunk[, classes$positive]),
    weight_label = "Shrunken centroid difference",
    signed = TRUE
  )
}
# </figure>
