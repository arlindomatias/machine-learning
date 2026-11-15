# ============================================================
# SCRIPT 01: the only dataset-specific step. Reads the raw input with the adapter
# inputs/<input.format>.R (Snakemake input `adapter`) and writes the three files of the input
# contract (see README.md in this folder and validate_inputs() in 00_common.R):
#   data/expression.rds, data/samples.csv, data/features.csv
#
# Each adapter defines
#   load_input(icfg, files, classes) -> list(X, samples, features)
# where icfg is the config's `input:` block, files the named raw input paths and classes
# cfg$classes. X is samples x features on the raw scale; samples has sample_id, group, class
# (raw labels) and any extra columns; features has feature_id, symbol, entrez.
# This script then maps class labels (classes.map), drops samples of other classes, applies
# input.transform, validates and writes.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))

icfg <- cfg$input
files <- snakemake@input[setdiff(names(snakemake@input), c("adapter", ""))]
source(inp$adapter)
d <- load_input(icfg, files, cfg$classes)
X <- d$X
samples <- as.data.table(d$samples)
features <- as.data.table(d$features)
message("Raw input, samples x features: ", paste(dim(X), collapse = " x "))

# --- Classes ---------------------------------------------------------------------
cmap <- unlist(cfg$classes$map)
raw <- as.character(samples$class)
samples[, class := ifelse(raw %in% names(cmap), cmap[raw], raw)]
print(samples[, .N, by = .(class)])
keep <- samples$class %in% CLASS_LEVELS
if (any(!keep)) {
  message("Dropping ", sum(!keep), " samples whose class is not ", paste(CLASS_LEVELS, collapse = " / "), ":")
  print(samples[!keep, .N, by = class])
}
samples <- samples[keep]
X <- X[keep, , drop = FALSE]
samples[, class := factor(class, levels = CLASS_LEVELS)]
setcolorder(samples, intersect(c("sample_id", "group", "class", "set"), names(samples)))

# --- Transform -------------------------------------------------------------------
transform <- icfg$transform %||% "none"
X <- switch(transform,
  log2p1 = log2(X + 1),
  none = X,
  stop("Unknown input.transform: ", transform, " (log2p1 | none)")
)
message("Transform: ", transform, " | samples x features: ", paste(dim(X), collapse = " x "))

validate_inputs(X, samples, features)
saveRDS(X, out$expression)
fwrite(samples, out$samples)
fwrite(features, out$features)
