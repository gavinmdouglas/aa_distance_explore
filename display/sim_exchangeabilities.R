rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(grid)
library(ggplot2)

aa_map <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/aa_map.tsv.gz',
                     sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 3)

metrics_map <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  
tab <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/codon_exchangeabilities.tsv.gz',
                  sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

tab <- tab[, -which(colnames(tab) == 'num_nonsyn_single')]

tab$amino_acid <- aa_map[tab$amino_acid, 'short']

for (coln in colnames(tab)) {
  if (coln %in% rownames(metrics_map)) {
    colnames(tab)[which(colnames(tab) == coln)] <- metrics_map[coln, 'Clean']
  }
}

prefs_to_rm <- c("Grantham (update, recalc.)", "Grantham (update, min-max)", "Grantham (recalc.)", "Grantham (min-max)",
                 "Miyata (recalc.)", "Miyata (min-max)", "Miyata (update, min-max)", "Miyata (update)")
tab <- tab[, which(! colnames(tab) %in% prefs_to_rm)]

colnames(tab)[which(colnames(tab) == "Grantham (orig.)")] <- 'Grantham'
colnames(tab)[which(colnames(tab) == "Miyata (orig.)")] <- 'Miyata'

by_aa <- aggregate(. ~ amino_acid, data = tab, FUN = mean)
rownames(by_aa) <- by_aa$amino_acid
by_aa <- by_aa[, -1]

by_aa_scaled <- t(scale(by_aa))

by_aa_scaled["Graur index (inverted)", ] <- by_aa_scaled["Graur index (inverted)", ] * -1

raw_heatmap <- ComplexHeatmap::Heatmap(matrix = by_aa_scaled,
                             cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        column_names_rot = 45,
                        name = 'Scaled\nmean\nexchange-\nabilities')

# Save as PNG, with ggsave.
ggsave(plot = grid.grabExpr(draw(raw_heatmap)),
       filename = '~/Drive/ncsu/aa_selection/aa_manuscript/figures/Supp_sim_exchangeabilities.png',
       dpi = 300,
       width = 9,
       height = 5)
