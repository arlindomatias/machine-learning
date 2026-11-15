# ============================================================
# SCRIPT 03: descriptive BSS/WSS ranking on the full cohort + heatmap of the top genes.
# Not used by the classifiers (they rank genes on the training set only, script 02).
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))

samples <- read_samples(inp$samples)
X <- readRDS(inp$expression)
n_top <- cfg$feature_selection$n_top
neg <- cfg$classes$negative
pos <- cfg$classes$positive

ranking <- bss_wss(X, samples$class)
top <- merge(ranking[seq_len(n_top)], split_gene_ids(ranking$gene_id[seq_len(n_top)]),
             by = "gene_id", sort = FALSE)
top[, mean_log2_negative := colMeans(X[samples$class == neg, gene_id])]
top[, mean_log2_positive := colMeans(X[samples$class == pos, gene_id])]
top[, log2_fold_change := mean_log2_positive - mean_log2_negative]
setnames(top, c("mean_log2_negative", "mean_log2_positive"),
         paste0("mean_log2_", tolower(c(neg, pos))))
setorder(top, rank)
setcolorder(top, c("rank", "gene_id", "symbol", "entrez"))
write_tsv(top, out$top_genes)
print(head(top[, .(rank, symbol, bss_wss, log2_fold_change)], 20))

# <figure name="heatmap_top_genes">
# All normals plus an equal-sized random subset of tumors; row z-scores clipped to [-3, 3].
set.seed(SEED)
normals <- which(samples$class == neg)
tumors <- sample(which(samples$class == pos), length(normals))
idx <- c(normals, tumors)
mat <- t(X[idx, top$gene_id])
rownames(mat) <- top$symbol
mat <- pmin(pmax(t(scale(t(mat))), -3), 3)
png(out$heatmap, width = 2400, height = 3200, res = 300)
heatmap(mat, scale = "none", col = hcl.colors(100, "Blue-Red 3"), zlim = c(-3, 3),
        ColSideColors = CLASS_COLORS[as.character(samples$class[idx])],
        labCol = FALSE, cexRow = 0.45, main = sprintf("Top %d BSS/WSS genes", n_top))
legend("topright", legend = names(CLASS_COLORS), fill = CLASS_COLORS, border = NA, bty = "n", cex = 0.8)
invisible(dev.off())
# </figure>
