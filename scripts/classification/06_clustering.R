# ============================================================
# SCRIPT 06: unsupervised clustering of the top full-cohort genes with ONE method
# (Snakemake wildcard `method`: kmeans | pam). Each method writes its own files.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages(library(cluster))

method <- snakemake@wildcards$method
ccfg <- cfg$clustering
samples <- read_samples(inp$samples)
top <- fread(inp$top_genes)
X <- scale(readRDS(inp$expression)[, top$feature_id])

set.seed(SEED)
cl <- switch(method,
  kmeans = kmeans(X, centers = ccfg$k, nstart = ccfg$nstart)$cluster,
  pam = pam(X, ccfg$k, cluster.only = TRUE),
  stop("Unknown clustering method: ", method)
)

write_tsv(data.table(method = method, sample_id = samples$sample_id, class = samples$class, cluster = cl),
          out$assignments)
tab <- table(cluster = cl, class = samples$class)
print(tab)
ct <- as.data.table(tab)[, method := method]
setcolorder(ct, "method")
write_tsv(ct, out$contingency)

# <figure name="clusplot">
png(out$clusplot, width = 2400, height = 2400, res = 300)
clusplot(X, cl, color = TRUE, shade = TRUE, labels = 4, lines = 0,
         col.clus = grey.colors(ccfg$k, 0.4, 0.75), col.p = CLASS_COLORS[as.character(samples$class)],
         cex = 0.5, main = sprintf("%s (k = %d), points colored by true class", method, ccfg$k))
legend("bottomright", legend = names(CLASS_COLORS), col = CLASS_COLORS, pch = 1, bty = "n", cex = 0.8)
invisible(dev.off())
# </figure>
