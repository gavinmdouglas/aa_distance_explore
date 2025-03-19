rm(list = ls(all.names = TRUE))

source('~/Drive/ncsu/aa_selection/github_repo/folding_selection/aa_metrics/aa_pairwise_diff_func.R')

# Compute distances based on min-max transformed Atchley factors.

atchley <- read.table('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/literature_data/Atchley2005_factors.tsv.gz',
                      header = TRUE, sep = '\t',
                      stringsAsFactors = FALSE, row.names = 1)

atchley_minmax <- pairwise_diff_normalized(tab = atchley,
                                           normtype = "minmax",
                                           sig_digits=4)

write.table(x = atchley_minmax,
            file = '/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/distances/atchley_minmax.tsv',
              sep = '\t',
              col.names = NA,
              row.names = TRUE,
              quote = FALSE)
