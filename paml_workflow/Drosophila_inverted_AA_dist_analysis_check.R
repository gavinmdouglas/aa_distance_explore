rm(list = ls(all.names = TRUE))

# Compare PAML LnL ranks across AA distance metrics for those that are inverted vs. not
# (sanity check to make sure there's not something obviously wonky about some of the metrics and that
# maybe they were mistakingly inverted).

library(ggplot2)
library(ggbeeswarm)
library(ggtext)

intab <- read.table(file="/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/Drosophila_normal_and_inverted_check.tsv.gz",
                    header=TRUE, sep="\t", stringsAsFactors = FALSE)

normal_tab <- intab[intab$lineage == "Drosophila_normal", ]
inverted_tab <- intab[intab$lineage == "Drosophila_inverted", ]


intab$distance_type <- paste(intab$lineage, intab$distance_type, sep = ' ')

intab$LnL_rank <- NA

for (subsample in unique(intab$subsample)) {
  subsample_df <- intab[intab$subsample == subsample, ]
  subsample_df <- subsample_df[order(subsample_df$lnL, decreasing = TRUE), ]
  subsample_df$LnL_rank <- 1:nrow(subsample_df)
  intab[intab$subsample == subsample, ] <- subsample_df
}

median_rank_by_type <- aggregate(LnL_rank ~ distance_type, data=intab, FUN=median)
intab$distance_type <- factor(intab$distance_type,
                              levels = median_rank_by_type$distance_type[order(median_rank_by_type$LnL_rank, decreasing = TRUE)])

ggplot(data=intab, aes(x = LnL_rank, y = distance_type)) +
  geom_quasirandom(col='black', size=0.8, orientation='y') +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_markdown()) +
  xlab('Model rank (lower is better)') +
  ylab('Amino acid distance basis for model')


inverted_medians <- median_rank_by_type[grep('inverted', median_rank_by_type$distance_type), 'LnL_rank']
names(inverted_medians) <- median_rank_by_type[grep('inverted', median_rank_by_type$distance_type), 'distance_type']

normal_medians <- median_rank_by_type[grep('normal', median_rank_by_type$distance_type), 'LnL_rank']
names(normal_medians) <- median_rank_by_type[grep('normal', median_rank_by_type$distance_type), 'distance_type']
