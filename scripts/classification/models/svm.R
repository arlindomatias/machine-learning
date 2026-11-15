# ============================================================
# Method svm: support vector machine (e1071::svm, C-classification).
# Config: kernel (linear | radial; default linear), cost (default 1), gamma (radial only;
#         default 1 / n_features), balanced (class weights inversely proportional to class size;
#         default true).
# Genes are scaled with training-set statistics (e1071 default). Score = Platt-scaled probability
# of the positive class. Linear kernel only: gene weight = hyperplane coefficient on scaled genes;
# positive = towards the positive class.
# Sourced by 04_classification.R; see that script for the fit_model() contract.
# ============================================================
suppressPackageStartupMessages(library(e1071))

# <figure name="method_svm">
fit_model <- function(X_train, y_train, X_test, mcfg, classes, groups) {
  kernel <- mcfg$kernel %||% "linear"
  args <- list(x = X_train, y = y_train, type = "C-classification", kernel = kernel,
               cost = mcfg$cost %||% 1, probability = TRUE)
  if (kernel == "radial" && !is.null(mcfg$gamma)) args$gamma <- mcfg$gamma
  if (isTRUE(mcfg$balanced %||% TRUE)) {
    n <- table(y_train)
    args$class.weights <- setNames(as.numeric(length(y_train) / (length(n) * n)), names(n))
  }
  fit <- do.call(svm, args)
  print(fit)
  pred <- predict(fit, X_test, probability = TRUE)
  res <- list(fit = fit, pred = pred, prob_pos = attr(pred, "probabilities")[, classes$positive])

  if (kernel == "linear") {
    # Decision value = w . x_scaled - rho. Its sign convention follows the order in which classes
    # appear in the data, so orient w so that the positive class has the higher training-set score.
    w <- drop(t(fit$coefs) %*% fit$SV)
    dv <- scale(X_train, fit$x.scale[[1]], fit$x.scale[[2]])[, names(w)] %*% w
    if (mean(dv[y_train == classes$positive]) < mean(dv[y_train == classes$negative])) w <- -w
    res <- c(res, list(genes = data.table(feature_id = names(w), weight = w),
                       weight_label = "Hyperplane weight (scaled genes)", signed = TRUE))
  }
  res
}
# </figure>
