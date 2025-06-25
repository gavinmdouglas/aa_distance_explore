rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ggbeeswarm)

intab <- read.table(file="/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/combined_output.tsv.gz",
                    header=TRUE, sep="\t", stringsAsFactors = FALSE)

metric_id_map <- read.table('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                            header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names=1)

intab$distance_type <- sub('-', '.', intab$distance_type)
intab$distance_type <- sub('_prob', '', intab$distance_type)
intab$distance_type <- sub('_maxscaled', '', intab$distance_type)

intab$distance_type <- metric_id_map[intab$distance_type, 'Clean']

miyata_grantham_i <- grep('Miyata|Grantham', intab$distance_type)

miyata_grantham_nonorig_i <- setdiff(miyata_grantham_i, grep("orig.", intab$distance_type))

intab <- intab[-miyata_grantham_nonorig_i, ]

# Restrict this analysis to subsamples 1-10.
focal_subsamples <- paste("subsample", 1:10, sep = '')
intab <- intab[intab$subsample %in% focal_subsamples, ]
intab$LnL_rank <- NA

intab[which(intab$lineage == 'Strep'), 'lineage'] <- 'Streptococcus'

# Get rank of each distance matrix by LnL, per lineage and subsample.
for (lineage in unique(intab$lineage)) {
  lineage_df <- intab[intab$lineage == lineage, ]
  for (subsample in unique(lineage_df$subsample)) {
    subsample_df <- lineage_df[lineage_df$subsample == subsample, ]
    subsample_df <- subsample_df[order(subsample_df$lnL, decreasing = TRUE), ]
    subsample_df$LnL_rank <- 1:nrow(subsample_df)
    intab[intab$lineage == lineage & intab$subsample == subsample, ] <- subsample_df
  }
}

intab$distance_type <- gsub(' \\(orig.\\)', '', intab$distance_type)

dist_means <- aggregate(LnL_rank ~ distance_type, data=intab, FUN=mean)
dist_means <- dist_means[order(dist_means$LnL_rank, decreasing = TRUE), ]

intab$distance_type <- factor(intab$distance_type, levels = dist_means$distance_type)

AADist_PAML_LnL_rank <- ggplot(data=intab, aes(x = LnL_rank, y = distance_type)) +
  geom_quasirandom(col='grey', size=0.8) +
  geom_boxplot(alpha=0.7, outlier.shape = NA) +
  facet_wrap(~lineage, nrow = 1) +
  theme_bw() +
  theme(position_quasirandom(orientation = 'x')) +
  xlab('Model rank (lower is better)') +
  ylab('Amino acid distance basis for model')

ggsave(plot = AADist_PAML_LnL_rank,
       filename = '~/AADist_PAML_LnL_rank_subsample1-10.pdf',
       width = 8, height = 6, useDingbats = FALSE)

# Get median a and b values per AA distance type.
median_vals_raw <- list()
for (disttype in unique(lineage_df$distance_type)) {
  subset_df <- intab[intab$distance_type == disttype, ]
  a_std <- (subset_df$a - mean(subset_df$a)) / sd(subset_df$a)
  b_std <- (subset_df$b - mean(subset_df$b)) / sd(subset_df$b)
  std_dist <- sqrt(a_std ** 2 + b_std ** 2)
  median_idx <- which.min(abs(std_dist - median(std_dist)))
  median_vals_raw[[disttype]] <- subset_df[median_idx, ]
}

median_vals <- do.call(rbind, median_vals_raw)
median_vals <- median_vals[, c("distance_type", "a", "b")]
rownames(median_vals) <- NULL

write.table(x = median_vals,
            file = "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/median_a_and_b_values.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
