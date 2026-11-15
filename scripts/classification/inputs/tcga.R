# ============================================================
# TCGA helpers, shared by inputs/firehose_rsem.R and 08_copy_number.R.
# ============================================================

# TCGA barcode: TCGA-<TSS>-<participant>-<sample type><vial>-...
# Sample-type codes 01-09 are tumors (01 primary, 06 metastatic), 10-19 normals
# (10 blood, 11 solid tissue normal). group = participant (TCGA-<TSS>-<participant>).
tcga_sample_metadata <- function(barcodes, classes) {
  parts <- tstrsplit(barcodes, "-", fixed = TRUE)
  type_code <- as.integer(substr(parts[[4]], 1, 2))
  data.table(
    sample_id = barcodes,
    group = paste(parts[[1]], parts[[2]], parts[[3]], sep = "-"),
    type_code = type_code,
    class = ifelse(type_code < 10, classes$positive, classes$negative)
  )
}

# Firehose gene ids "FHL1|2273" -> symbol "FHL1", entrez "2273". Unannotated genes
# ("?|100130426") get the Entrez id as their symbol.
firehose_features <- function(ids) {
  symbol <- sub("\\|.*$", "", ids)
  entrez <- sub("^.*\\|", "", ids)
  symbol[symbol == "?"] <- entrez[symbol == "?"]
  data.table(feature_id = ids, symbol = symbol, entrez = entrez)
}
