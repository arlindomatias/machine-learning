# ============================================================
# SCRIPT 07: Reactome over-representation of the top full-cohort genes
# (universe = every gene measured on the platform).
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages({
  library(ReactomePA)
  library(enrichplot)
})

ecfg <- cfg$enrichment
top <- fread(inp$top_genes)
universe <- split_gene_ids(colnames(readRDS(inp$expression)))$entrez

ora <- enrichPathway(gene = as.character(top$entrez), universe = universe, organism = "human",
                     pvalueCutoff = ecfg$padj_cutoff, pAdjustMethod = "BH", readable = TRUE)
res <- if (is.null(ora)) data.table() else as.data.table(as.data.frame(ora))
write_tsv(res, out$table)
message("Reactome pathways with adjusted p < ", ecfg$padj_cutoff, ": ", nrow(res))

# <figure name="reactome_dotplot">
if (nrow(res) >= 1) {
  save_plot(out$dotplot, dotplot(ora, showCategory = 20) + ggtitle("Reactome over-representation"),
            width = 8, height = 7)
}
# </figure>

# <figure name="reactome_emapplot">
if (nrow(res) >= 2) {
  ora <- pairwise_termsim(ora)
  save_plot(out$emapplot, emapplot(ora, showCategory = ecfg$show_category), width = 9, height = 9)
}
# </figure>

# Steps without results still create their declared outputs.
for (f in c(out$dotplot, out$emapplot)) if (!file.exists(f)) {
  png(f, width = 600, height = 200); plot.new(); text(0.5, 0.5, "No enriched pathways"); invisible(dev.off())
}
