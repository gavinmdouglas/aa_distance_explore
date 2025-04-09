rm(list = ls(all.names = TRUE))

# Parse scales in AAOntology databae (mainly from the AAindex database),
# and run PCA on scales within each category to get key components.

in_scales <- read.table('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/AAOntology_key_files/AAOntology_subset_Supp_Table1_normalized_scales.tsv',
                        header=TRUE, sep = '\t', row.names=1)

scale_info <- read.table('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/AAOntology_key_files/AAOntology_subset_Supp_Table3_scale_info.tsv',
                         header=TRUE, sep = '\t', row.names=1)

# Confirm that scale IDs match exactly.
identical(sort(colnames(in_scales)), sort(rownames(scale_info)))

scale_info$category <- sub("/", "_", scale_info$category)
scale_info$category <- sub("-", "_", scale_info$category)

scale_categories <- unique(scale_info$category)

for (scale_category in scale_categories) {
  category_scales <- rownames(scale_info[scale_info$category == scale_category, ])
  category_scale_subset <- in_scales[, category_scales]
  
  write.table(x = category_scale_subset,
              file=gzfile(paste0('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/AAOntology_key_files/scales_by_category/', scale_category, '.tsv.gz')),
              sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE)

  category_pca <- prcomp(category_scale_subset, center = TRUE, scale. = TRUE)
  
  category_PCs <- data.frame(category_pca$x)
  
  eigenvalues <- category_pca$sdev ** 2
  
  PC_variance_explained <- data.frame(PC=1:length(eigenvalues),
                                      eigenvalue=eigenvalues / sum(eigenvalues),
                                      prop_variance_explained=eigenvalues / sum(eigenvalues))

  write.table(x = category_PCs,
              file=gzfile(paste0('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/AAOntology_key_files/PCs_by_category/PCs_', scale_category, '.tsv.gz')),
              sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE)
  
  write.table(x = PC_variance_explained,
              file=gzfile(paste0('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/AAOntology_key_files/PCs_by_category/PCs_', scale_category, '_var_explained.tsv.gz')),
              sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
}
