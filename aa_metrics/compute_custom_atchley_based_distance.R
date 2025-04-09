rm(list = ls(all.names = TRUE))

# Compute distances based on min-max transformed Atchley factors. In this case, all factors are equally weighed.

source('~/Drive/research/aa_distance/aa_distance_explore/aa_metrics/aa_pairwise_diff_func.R')

atchley <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/Atchley2005_factors.tsv.gz',
                      header = TRUE, sep = '\t',
                      stringsAsFactors = FALSE, row.names = 1)
atchley_minmax <- pairwise_diff_normalized(tab = atchley,
                                           normtype = "minmax",
                                           sig_digits=4)
write.table(x = atchley_minmax,
            file = '~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/distances/atchley_minmax.tsv',
            sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE)