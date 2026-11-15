# R session of the analysis environment (every package used by the pipeline).
source(file.path(snakemake@scriptdir, "00_common.R"))
suppressPackageStartupMessages({
  library(class); library(randomForest); library(cluster); library(ReactomePA); library(enrichplot)
})
writeLines(capture.output(sessionInfo()), out$session)
