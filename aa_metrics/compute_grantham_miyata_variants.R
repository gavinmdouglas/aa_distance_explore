rm(list = ls(all.names = TRUE))

source('~/Drive/ncsu/aa_selection/github_repo/folding_selection/aa_metrics/aa_pairwise_diff_func.R')

# Compute variants of Grantham and miyata distances.

# First of all, re-compute based on original values in table.
# There were some minor errors in the original Grantham for instance, so worth
# doing, even with mean normalization approach.
# Get re-computed miyata too (for SD approach), in case there are differences
# with that one too.

# Then use updated polarity and volume vectors to get the Grantham and Miyata-like values, 
# based on mean and SD normalization, respectively, but also with min-max distances.

grantham_properties <- read.table('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/literature_data/Grantham1974_Table1_properties.tsv',
                                  header = TRUE, sep = '\t',
                                  stringsAsFactors = FALSE, row.names = 1)

grantham_properties <- grantham_properties[, 3:5]

grantham_mean_norm <- pairwise_diff_normalized(tab = grantham_properties,
                                               normtype = "mean_diff",
                                               set_mean_100 = TRUE)

grantham_minmax <- pairwise_diff_normalized(tab = grantham_properties,
                                            normtype = "minmax",
                                            set_mean_100 = TRUE)

miyata_sd_norm <- pairwise_diff_normalized(tab = grantham_properties[, c('polarity', 'volume')],
                                           normtype = "sd_diff")

miyata_minmax <- pairwise_diff_normalized(tab = grantham_properties[, c('polarity', 'volume')],
                                          normtype = "minmax")

# Then read in updated variables.
polarity <- read.table('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/literature_data/polarity_mean_minmax.tsv',
                       header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

volume <- read.table('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/literature_data/Tsai2002_Table3_volumes.tsv',
                       header = TRUE, sep = '\t', stringsAsFactors = FALSE)
volume <- volume[which(! is.na(volume$AA)), ]
rownames(volume) <- volume$AA

updated_properties <- data.frame(c = grantham_properties$composition,
                                 p = polarity[rownames(grantham_properties), 'mean_minmax'],
                                 v = volume[rownames(grantham_properties), 'ProtOr_volume'])
rownames(updated_properties) <- rownames(grantham_properties)

custom_grantham_mean_norm <- pairwise_diff_normalized(tab = updated_properties,
                                                      normtype = "mean_diff",
                                                      set_mean_100 = TRUE)

custom_grantham_minmax <- pairwise_diff_normalized(tab = updated_properties,
                                                  normtype = "minmax",
                                                  set_mean_100 = TRUE)

custom_miyata_sd_norm <- pairwise_diff_normalized(tab = updated_properties[, c('p', 'v')],
                                                  normtype = "sd_diff")

custom_miyata_minmax <- pairwise_diff_normalized(tab = updated_properties[, c('p', 'v')],
                                                normtype = "minmax")

# Add all these output tables to a single list, with names as keys.
output_tabs <- list(
  'grantham_mean_norm' = grantham_mean_norm,
  'grantham_minmax' = grantham_minmax,
  'miyata_sd_norm' = miyata_sd_norm,
  'miyata_minmax' = miyata_minmax,
  'custom_grantham_mean_norm' = custom_grantham_mean_norm,
  'custom_grantham_minmax' = custom_grantham_minmax,
  'custom_miyata_sd_norm' = custom_miyata_sd_norm,
  'custom_miyata_minmax' = custom_miyata_minmax)

for (key in names(output_tabs)) {
  write.table(x = output_tabs[[key]],
              file = paste0('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/distances/',
                            key, '.tsv'),
              sep = '\t',
              col.names = NA,
              row.names = TRUE,
              quote = FALSE)
}
