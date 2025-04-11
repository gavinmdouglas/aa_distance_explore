rm(list = ls(all.names = TRUE))

PC_var_explained_files <- list.files('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/AAOntology_key_files/PCs_by_category',
                                     full.names = TRUE, pattern = '_var_explained.tsv.gz$')

for (PC_var_explained_file in PC_var_explained_files) {
  feature_type <- gsub('_var_explained.tsv.gz', '', basename(PC_var_explained_file))
  PC_file <- sub('_var_explained.tsv.gz$', '.tsv.gz', PC_var_explained_file)
  PC_tab <- read.table(PC_file, header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)
  PC_var_explained <- read.table(PC_var_explained_file, header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

  PC_dist <- matrix(NA, nrow = nrow(PC_tab), ncol = nrow(PC_tab))
  rownames(PC_dist) <- rownames(PC_tab)
  colnames(PC_dist) <- rownames(PC_tab)
  for (i in 1:nrow(PC_tab)) {
    for (j in 1:nrow(PC_tab)) {
      PC_dist[i, j] <- sqrt(sum(((PC_tab[i, ] - PC_tab[j, ]) ^ 2) * PC_var_explained$prop_variance_explained))
    }
  }
  
  write.table(x = data.frame(PC_dist),
              file = gzfile(paste0('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/distances/', 'AAOntology-', feature_type, '.tsv.gz')),
              sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE)
}
