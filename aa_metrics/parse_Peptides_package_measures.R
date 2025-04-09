rm(list = ls(all.names = TRUE))

# The Peptides R package handily comes with several amino acid characteristics.
# Went through and wrote these to tables for consistent processing.

library("Peptides")

data(AAdata)
?AAdata

parse_aa_vectors_to_df <- function(in_list) {
  exp_aa_order <- names(in_list[[1]])
  if (length(exp_aa_order) != 20) { stop("Expected 20 amino acids!") }
  prepped <- data.frame(matrix(NA, nrow=20, ncol=0))
  prepped$amino_acid <- exp_aa_order
  for (feature in names(in_list)) {
    if (length(in_list[[feature]]) != 20) { stop("Expected 20 amino acids!") }
    if (! identical(names(in_list[[feature]]), exp_aa_order)) {
      stop("Amino acid order is not the same!")
    }
    prepped[, feature] <- in_list[[feature]]
  }

  return(prepped)
}

Cruciani2004 <- parse_aa_vectors_to_df(AAdata$crucianiProperties)
Kidera2008 <- parse_aa_vectors_to_df(AAdata$kideraFactors)
zScales <- parse_aa_vectors_to_df(AAdata$zScales)
FASAGI <- parse_aa_vectors_to_df(AAdata$FASGAI)
VHSE <- parse_aa_vectors_to_df(AAdata$VHSE)

# Hydrophobicity was notably left out.
# Decided not to analyze these, as they are reflected in other metrics,
# esepcially with AAIndex clustersings.

write.table(x = Cruciani2004,
            file=gzfile('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/Peptides_R_package_data/Cruciani2004.tsv.gz'),
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(x = Kidera2008,
            file=gzfile('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/Peptides_R_package_data/Kidera2008.tsv.gz'),
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(x = zScales,
            file=gzfile('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/Peptides_R_package_data/zScales.tsv.gz'),
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(x = FASAGI,
            file=gzfile('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/Peptides_R_package_data/FASAGI.tsv.gz'),
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(x = VHSE,
            file=gzfile('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/Peptides_R_package_data/VHSE.tsv.gz'),
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
