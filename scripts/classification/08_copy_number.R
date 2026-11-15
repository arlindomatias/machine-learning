# ============================================================
# SCRIPT 08 (optional module `copy_number`, TCGA Firehose data only): somatic copy-number segments
# (Affymetrix SNP 6.0, germline CNVs removed). Segment_Mean is log2(copy number / 2).
# Config modules.copy_number: file (under $DATA/<input.dir>/), gain_loss_threshold, example_sample.
# ============================================================
source(file.path(snakemake@scriptdir, "00_common.R"))
source(file.path(snakemake@scriptdir, "inputs", "tcga.R"))

ccfg <- cfg$modules$copy_number
thr <- ccfg$gain_loss_threshold
cnv <- fread(inp$copy_number)
message("Segments x columns: ", paste(dim(cnv), collapse = " x "))
print(summary(cnv$Segment_Mean))

cnv[, length := End - Start + 1]
per_sample <- cnv[, .(
  n_segments = .N,
  n_positive = sum(Segment_Mean > 0),
  n_negative = sum(Segment_Mean < 0),
  mean_segment = mean(Segment_Mean),
  frac_genome_gained = sum(length[Segment_Mean > thr]) / sum(length),
  frac_genome_lost = sum(length[Segment_Mean < -thr]) / sum(length)
), by = .(sample_id = Sample)]
per_sample[, frac_genome_altered := frac_genome_gained + frac_genome_lost]
per_sample <- merge(tcga_sample_metadata(per_sample$sample_id, cfg$classes), per_sample, by = "sample_id")
per_sample[, class := factor(class, levels = CLASS_LEVELS)]
write_tsv(per_sample, out$per_sample)

message("\nExample sample ", ccfg$example_sample, ":")
print(per_sample[sample_id == ccfg$example_sample])
print(per_sample[, .(n = .N, median_fga = median(frac_genome_altered)), by = .(type_code, class)])

# <figure name="fraction_genome_altered">
per_sample[, tissue := fcase(type_code == 10, "Blood normal", type_code == 11, "Solid tissue normal",
                             type_code == 1, "Primary tumor", type_code == 6, "Metastasis",
                             default = paste("Type", type_code))]
p <- ggplot(per_sample, aes(tissue, frac_genome_altered, fill = class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.3) +
  scale_fill_manual(values = CLASS_COLORS, name = NULL) +
  labs(x = NULL, y = "Fraction of genome altered",
       title = "Somatic copy-number alteration burden",
       subtitle = sprintf("Segments with |log2 ratio| > %.1f", thr))
save_plot(out$fga_plot, p, width = 6.5, height = 4.5)
# </figure>
