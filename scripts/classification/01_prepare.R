# ============================================================
# SCRIPT 01: RNA-seq expression -> log2(RSEM + 1) matrix (samples x genes) + sample table
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))

path <- inp$expression
# Line 1: "Hybridization REF" + barcodes; line 2: "gene_id" + "normalized_count" labels.
barcodes <- scan(path, what = "", sep = "\t", nlines = 1, quiet = TRUE)[-1]
dt <- fread(path, skip = 2, header = FALSE, sep = "\t")
genes <- dt[[1]]
X <- t(as.matrix(dt[, -1]))
dimnames(X) <- list(barcodes, genes)
rm(dt)
X <- log2(X + 1)
message("Samples x genes: ", paste(dim(X), collapse = " x "))

samples <- sample_metadata(barcodes)
print(samples[, .N, by = .(type_code, class)][order(type_code)])

saveRDS(X, out$expression)
fwrite(samples, out$samples)
