# ============================================================
# SCRIPT 02: patient-level train/test split + BSS/WSS gene ranking on the training set
# A patient's tumor and normal samples always fall on the same side; the split is stratified
# by whether the patient contributed a normal sample.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))

samples <- read_samples(inp$samples)
X <- readRDS(inp$expression)
stopifnot(identical(rownames(X), samples$barcode))

has_normal <- samples[, .(has_normal = any(class == cfg$classes$negative)), by = patient]
set.seed(SEED)
test_patients <- unlist(lapply(split(has_normal$patient, has_normal$has_normal), function(p) {
  sample(p, round(cfg$split$test_fraction * length(p)))
}))
samples[, set := ifelse(patient %in% test_patients, "test", "train")]
print(samples[, .N, by = .(set, class)][order(set, class)])
fwrite(samples[, .(barcode, patient, class, set)], out$split)

train <- samples$set == "train"
ranking <- bss_wss(X[train, ], samples$class[train])
ranking <- merge(ranking, split_gene_ids(ranking$gene_id), by = "gene_id", sort = FALSE)
setorder(ranking, rank)
setcolorder(ranking, c("rank", "gene_id", "symbol", "entrez", "bss_wss"))
write_tsv(ranking, out$ranking)
message("Genes ranked on the training set: ", nrow(ranking))
