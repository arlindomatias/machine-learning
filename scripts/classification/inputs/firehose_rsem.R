# ============================================================
# Input format firehose_rsem: Broad GDAC Firehose level-3 RSEM gene expression
# (e.g. *rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt).
# Config (input:): expression = file name under $DATA/<input.dir>/.
# Samples come from the TCGA barcode (group = participant, class from the sample-type code,
# type_code kept); features from the "SYMBOL|ENTREZ" ids. Values are on the raw RSEM scale
# (use transform: log2p1).
# Sourced by 01_prepare.R; see that script for the load_input() contract.
# ============================================================
source(file.path(snakemake@scriptdir, "inputs", "tcga.R"))

load_input <- function(icfg, files, classes) {
  path <- files$expression
  # Line 1: "Hybridization REF" + barcodes; line 2: "gene_id" + "normalized_count" labels.
  barcodes <- scan(path, what = "", sep = "\t", nlines = 1, quiet = TRUE)[-1]
  dt <- fread(path, skip = 2, header = FALSE, sep = "\t")
  ids <- dt[[1]]
  X <- t(as.matrix(dt[, -1]))
  dimnames(X) <- list(barcodes, ids)
  list(X = X, samples = tcga_sample_metadata(barcodes, classes), features = firehose_features(ids))
}
