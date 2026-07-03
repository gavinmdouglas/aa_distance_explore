rm(list = ls(all.names = TRUE))

aa_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/aa_map.tsv.gz',
                     sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 3)
aa_map$long_w_letter <- stringr::str_to_title(paste0(aa_map$full, ' (', rownames(aa_map), ')'))
aa_map$long_w_letter <- sub(' Acid ', ' acid ', aa_map$long_w_letter)

metrics_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

tab <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/codon_exchangeabilities.tsv.gz',
                  sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

tab <- tab[, -which(colnames(tab) == 'num_nonsyn_single')]


tab$amino_acid <- aa_map[tab$amino_acid, 'long_w_letter']

for (coln in colnames(tab)) {
  if (coln %in% rownames(metrics_map)) {
    colnames(tab)[which(colnames(tab) == coln)] <- metrics_map[coln, 'Clean']
  }
}

prefs_to_rm <- c("Grantham (update, recalc.)", "Grantham (update, min-max)", "Grantham (recalc.)", "Grantham (min-max)",
                 "Miyata (recalc.)", "Miyata (min-max)", "Miyata (update, min-max)", "Miyata (update)", "Graur stability index",
                 "DMS-EX (mean alternative)")

prefs_to_rm <- c(prefs_to_rm, grep("\\+", colnames(tab), value=TRUE))

tab <- tab[, which(! colnames(tab) %in% prefs_to_rm)]

colnames(tab)[which(colnames(tab) == "Grantham (orig.)")] <- 'Grantham'
colnames(tab)[which(colnames(tab) == "Miyata (orig.)")] <- 'Miyata'

by_aa <- aggregate(. ~ amino_acid, data = tab, FUN = mean)
rownames(by_aa) <- by_aa$amino_acid
by_aa <- by_aa[, -1]

by_aa_scaled <- t(scale(by_aa, center = TRUE, scale = TRUE))

# Included stability index as a test originally.
# by_aa_scaled["Graur stability index", ] <- by_aa_scaled["Graur stability index", ] * -1

write.table(x = by_aa_scaled,
            file=gzfile('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/codon_exchangeabilities_prepped_scaled.tsv.gz'),
            sep = '\t', col.names = NA, row.names = TRUE, quote = FALSE)
