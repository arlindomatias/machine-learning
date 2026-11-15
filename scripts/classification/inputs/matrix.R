# ============================================================
# Input format matrix: any delimited expression matrix plus a sample sheet and, optionally, a
# feature annotation table. Config (input:), file names under $DATA/<input.dir>/:
#   expression:   matrix file. The first column holds the row ids.
#   orientation:  features_x_samples (default; rows = features, as in GEO) | samples_x_features
#   samples:      {file, columns: {sample_id, class, group}}. group is optional (default: each sample
#                 is its own group). Every other column of the sheet is kept. A `set` column
#                 (train/test), if present, fixes the split.
#   features:     optional {file, columns: {feature_id, symbol, entrez}}. symbol and entrez are
#                 optional. Without it, symbol = feature_id and entrez = NA. Entrez ids
#                 that map to several genes ("1 /// 2", GEO style) are set to NA.
# Samples in the matrix but not in the sheet are dropped (and logged); sheet rows not in the matrix
# are an error. Class labels are mapped by classes.map in 01_prepare.R.
# Sourced by 01_prepare.R; see that script for the load_input() contract.
# ============================================================

load_input <- function(icfg, files, classes) {
  dt <- fread(files$expression, check.names = FALSE)
  ids <- as.character(dt[[1]])
  X <- as.matrix(dt[, -1])
  rownames(X) <- ids
  if ((icfg$orientation %||% "features_x_samples") == "features_x_samples") X <- t(X)
  storage.mode(X) <- "double"

  scfg <- icfg$samples
  sheet <- fread(files$samples, colClasses = "character")
  cols <- scfg$columns
  for (k in c("sample_id", "class")) if (is.null(cols[[k]])) stop("input.samples.columns needs ", k)
  samples <- copy(sheet)
  setnames(samples, cols$sample_id, "sample_id")
  setnames(samples, cols$class, "class")
  if (!is.null(cols$group)) setnames(samples, cols$group, "group") else samples[, group := sample_id]
  missing <- setdiff(samples$sample_id, rownames(X))
  if (length(missing)) stop("Samples in the sheet but not in the matrix: ", paste(head(missing), collapse = ", "))
  extra <- setdiff(rownames(X), samples$sample_id)
  if (length(extra)) message("Dropping ", length(extra), " matrix samples that are not in the sample sheet")
  X <- X[samples$sample_id, , drop = FALSE]

  features <- data.table(feature_id = colnames(X))
  if (!is.null(icfg$features)) {
    fcols <- icfg$features$columns
    ann <- fread(files$features, colClasses = "character")
    ann <- ann[, c(fcols$feature_id, fcols$symbol, fcols$entrez), with = FALSE]
    setnames(ann, c(fcols$feature_id, fcols$symbol, fcols$entrez),
             c("feature_id", if (!is.null(fcols$symbol)) "symbol", if (!is.null(fcols$entrez)) "entrez"))
    features <- merge(features, unique(ann, by = "feature_id"), by = "feature_id", all.x = TRUE, sort = FALSE)
    features <- features[match(colnames(X), feature_id)]
  }
  if (!"symbol" %in% names(features)) features[, symbol := feature_id]
  features[is.na(symbol) | symbol == "", symbol := feature_id]
  if (!"entrez" %in% names(features)) features[, entrez := NA_character_]
  features[entrez == "", entrez := NA_character_]
  features[grepl("///", entrez, fixed = TRUE), entrez := NA_character_]   # ambiguous (GEO multi-gene probes)
  list(X = X, samples = samples, features = features)
}
