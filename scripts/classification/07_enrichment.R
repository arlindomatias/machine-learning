# ============================================================
# SCRIPT 07 (optional module `enrichment`): Reactome over-representation of the top full-cohort
# features. Uses the Entrez ids in features.csv; universe = every measured feature with an Entrez id.
# Config modules.enrichment: organism (ReactomePA: human | rat | mouse | celegans | yeast | zebrafish
# | fly), padj_cutoff, show_category.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages({
  library(ReactomePA)
  library(enrichplot)
})

ecfg <- cfg$modules$enrichment
top <- fread(inp$top_genes, colClasses = list(character = "entrez"))
features <- read_features(inp$features)
universe <- unique(na.omit(features$entrez))
genes <- unique(na.omit(top$entrez[top$entrez != ""]))
message("Top features with an Entrez id: ", length(genes), " of ", nrow(top), " | universe: ", length(universe))
if (length(universe) == 0) stop("No Entrez ids in features.csv: disable modules.enrichment or annotate the features")

ora <- enrichPathway(gene = genes, universe = universe, organism = ecfg$organism %||% "human",
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
