rm(list = ls(all.names = TRUE))

# Read in components from Venkatarajan et al. 2001, and get distance measure by averaging distances (weighted by eigenvalues).

eigenvalues <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/Venkatarajan2001_eigenvalues.tsv",
                          header=TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

eigenvalues$prop <- eigenvalues$Eigenvalue / sum(eigenvalues$Eigenvalue)

PC_tab <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/Venkatarajan2001.tsv",
                      header=TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

PC_dist <- matrix(NA, nrow = nrow(PC_tab), ncol = nrow(PC_tab))
rownames(PC_dist) <- rownames(PC_tab)
colnames(PC_dist) <- rownames(PC_tab)

for (i in 1:nrow(PC_tab)) {
  for (j in 1:nrow(PC_tab)) {
    PC_dist[i, j] <- sqrt(sum(((PC_tab[i, ] - PC_tab[j, ]) ^ 2) * eigenvalues$prop))
  }
}

PC_dist$amino_a

write.table(x = PC_dist,
            file=gzfile('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/distances/Venkatarajan2001.tsv.gz'),
            sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE)
