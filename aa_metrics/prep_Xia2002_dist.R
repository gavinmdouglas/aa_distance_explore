rm(list = ls(all.names = TRUE))

# Back-calculate Xia 2002 distance from their Table 10, which is the mean of their distance with Miyata.
table10 <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/Xia2002_Table10.txt",
                     header=TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

miyata <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/distances/miyata_orig.tsv.gz",
                     header=TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

aa_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/aa_map.tsv.gz',
                     header=TRUE, sep = '\t', row.names = 2, stringsAsFactors = FALSE)

colnames(table10) <- aa_map[colnames(table10), 'single']
rownames(table10) <- aa_map[rownames(table10), 'single']
table10['A', ] <- NA

working <- data.frame(matrix(NA, nrow = nrow(miyata), ncol = nrow(miyata)))
rownames(working) <- rownames(miyata)
colnames(working) <- rownames(miyata)

for (aa1 in rownames(working)) {
  for (aa2 in colnames(working)) {
    if (aa1 == aa2) {
      working[aa1, aa2] <- 0.0
    } else if (! is.na(table10[aa1, aa2])) {
      working[aa1, aa2] <- table10[aa1, aa2]
    } else if (! is.na(table10[aa2, aa1])) {
      working[aa1, aa2] <- table10[aa2, aa1]
    } else {
      stop(paste("No entry for", aa1, aa2))
    }
  }
}

# Remove Miyata
working <- (working * 2) - miyata

write.table(x = working,
            file=gzfile('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/distances/Xia2002.tsv.gz'),
            sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE)
