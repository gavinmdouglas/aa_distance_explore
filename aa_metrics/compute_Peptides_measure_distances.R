rm(list = ls(all.names = TRUE))

# Compute AA distances based on additional metrics in Peptides database, using the same equation used for Grantham and Miyata metrics.

source('~/Drive/research/aa_distance/aa_distance_explore/aa_metrics/aa_pairwise_diff_func.R')

feature_table_files <- list.files('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/Peptides_R_package_data/',
                             full.names = TRUE, pattern = '.tsv.gz$')

for (feature_table_file in feature_table_files) {
  feature_type <- gsub('.tsv.gz', '', basename(feature_table_file))
  feature_tab <- read.table(feature_table_file, header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)
  feature_dist <- pairwise_diff_normalized(tab = feature_tab, normtype = "minmax")

  write.table(x = feature_dist,
              file = gzfile(paste0('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/distances/', feature_type, '.tsv.gz')),
              sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE)
}
