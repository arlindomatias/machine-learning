# ============================================================
# SCRIPT 00: shared setup (sourced by every step)
# Dataset-agnostic: every dataset-specific value comes from the project's config.yaml
# (`snakemake@config`). Scripts run through the Snakemake `script:` directive.
# ============================================================
if (!exists("snakemake")) {
  stop("Run this script through Snakemake (`snakemake --use-conda --cores 2`).")
}

# R messages and output go to the rule log.
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

cfg <- snakemake@config
out <- snakemake@output
inp <- snakemake@input

for (f in unlist(out)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

options(timeout = 600)
SEED <- cfg$seed %||% 123
set.seed(SEED)
if (!interactive()) pdf(NULL)

CLASS_LEVELS <- c(cfg$classes$negative, cfg$classes$positive)
CLASS_COLORS <- unlist(cfg$classes$colors)

theme_set(theme_bw(base_size = 11))

save_plot <- function(path, plot, ...) {
  tryCatch(ggsave(path, plot = plot, dpi = 300, bg = "white", ...),
           error = function(e) message("Plot ", path, " failed: ", conditionMessage(e)))
}

write_tsv <- function(x, path) fwrite(x, path, sep = "\t")

# --- Sample metadata ------------------------------------------------------------

# TCGA barcode: TCGA-<TSS>-<participant>-<sample type><vial>-...
# Sample-type codes 01-09 are tumors (01 primary, 06 metastatic), 10-19 normals
# (10 blood, 11 solid tissue normal).
sample_metadata <- function(barcodes) {
  parts <- tstrsplit(barcodes, "-", fixed = TRUE)
  type_code <- as.integer(substr(parts[[4]], 1, 2))
  data.table(
    barcode = barcodes,
    patient = paste(parts[[1]], parts[[2]], parts[[3]], sep = "-"),
    type_code = type_code,
    class = factor(ifelse(type_code < 10, cfg$classes$positive, cfg$classes$negative),
                   levels = CLASS_LEVELS)
  )
}

read_samples <- function(path) {
  s <- fread(path)
  s[, class := factor(class, levels = CLASS_LEVELS)]
  s
}

# "FHL1|2273" -> symbol "FHL1", entrez "2273". Unannotated genes ("?|100130426") get the
# Entrez ID as their symbol.
split_gene_ids <- function(ids) {
  symbol <- sub("\\|.*$", "", ids)
  entrez <- sub("^.*\\|", "", ids)
  symbol[symbol == "?"] <- entrez[symbol == "?"]
  data.table(gene_id = ids, symbol = symbol, entrez = entrez)
}

# --- Feature selection ----------------------------------------------------------

# Between-group / within-group sum-of-squares ratio per gene (Dudoit et al., 2002).
# Vectorized over genes; any number of classes. Constant genes are dropped.
bss_wss <- function(X, y) {
  y <- as.factor(y)
  overall <- colMeans(X)
  bss <- wss <- numeric(ncol(X))
  for (lv in levels(y)) {
    Xk <- X[y == lv, , drop = FALSE]
    mk <- colMeans(Xk)
    bss <- bss + nrow(Xk) * (mk - overall)^2
    wss <- wss + colSums(sweep(Xk, 2, mk)^2)
  }
  ratio <- bss / wss
  res <- data.table(gene_id = colnames(X), bss_wss = ratio)[is.finite(bss_wss)]
  setorder(res, -bss_wss)
  res[, rank := .I]
  res
}
