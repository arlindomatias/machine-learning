# ============================================================
# SCRIPT 02: group-level train/test split + BSS/WSS feature ranking on the training set.
# All samples of a group (e.g. a patient's tumor and normal) fall on the same side. The split is
# stratified by whether the group has a negative-class sample, so the rarer paired groups are
# spread over both sets. If samples.csv has a `set` column, that split is used instead.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))

samples <- read_samples(inp$samples)
features <- read_features(inp$features)
X <- read_expression(inp$expression)
stopifnot(identical(rownames(X), samples$sample_id))

if ("set" %in% names(samples)) {
  message("Using the predefined split in samples.csv")
} else {
  strata <- samples[, .(has_negative = any(class == cfg$classes$negative)), by = group]
  set.seed(SEED)
  test_groups <- unlist(lapply(split(strata$group, strata$has_negative), function(g) {
    sample(g, round(cfg$split$test_fraction * length(g)))
  }))
  samples[, set := ifelse(group %in% test_groups, "test", "train")]
}
print(samples[, .N, by = .(set, class)][order(set, class)])
fwrite(samples[, .(sample_id, group, class, set)], out$split)

train <- samples$set == "train"
ranking <- annotate_features(bss_wss(X[train, ], samples$class[train]), features, sort = FALSE)
setorder(ranking, rank)
setcolorder(ranking, c("rank", "feature_id", "symbol", "entrez", "bss_wss"))
write_tsv(ranking, out$ranking)
message("Features ranked on the training set: ", nrow(ranking))
