rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(grid)
library(cowplot)
library(ggplot2)
library(ggrepel)

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
                 "Miyata (recalc.)", "Miyata (min-max)", "Miyata (update, min-max)", "Miyata (update)", "Graur index (inverted)")
tab <- tab[, which(! colnames(tab) %in% prefs_to_rm)]

colnames(tab)[which(colnames(tab) == "Grantham (orig.)")] <- 'Grantham (Inverted Graur index)'
colnames(tab)[which(colnames(tab) == "Miyata (orig.)")] <- 'Miyata'

by_aa <- aggregate(. ~ amino_acid, data = tab, FUN = mean)
rownames(by_aa) <- by_aa$amino_acid
by_aa <- by_aa[, -1]

by_aa_scaled <- t(scale(by_aa, center = TRUE, scale = TRUE))

#by_aa_scaled["Graur index (inverted)", ] <- by_aa_scaled["Graur index (inverted)", ] * -1

rownames(by_aa_scaled)[which(rownames(by_aa_scaled) == "Grantham (Inverted Graur index)")] <- 'Grantham'

raw_heatmap <- ComplexHeatmap::Heatmap(matrix = by_aa_scaled,
                                       cluster_rows = TRUE,
                                       cluster_columns = TRUE,
                                       clustering_method_rows = 'complete',
                                       clustering_method_columns = 'complete',
                                       clustering_distance_rows = 'euclidean',
                                       clustering_distance_columns = 'euclidean',
                                       column_names_rot = 45,
                                       name = 'Scaled\nmean\nexchange-\nabilities')

# Pairwise distance matrix.
# by_aa_scaled_dist <- dist(by_aa_scaled, method = "euclidean")
# Heatmap(as.matrix(by_aa_scaled_dist), rect_gp = gpar(type = "none"),
#         show_column_dend = FALSE,
#         row_labels = rownames(by_aa_scaled),
#         column_labels = rownames(by_aa_scaled),
#         column_names_rot = 45,
#         cell_fun = function(j, i, x, y, w, h, fill) {
#           if(as.integer(x) <= 1 - as.integer(y)) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#           }
#         })

# Run PCA
PCA_metrics <- prcomp(by_aa_scaled, center = TRUE, scale. = TRUE)
eigenvalues <- PCA_metrics$sdev ** 2
percent_explained <- sprintf('%.1f', eigenvalues / sum(eigenvalues) * 100)

PC1_vs_PC2 <- ggplot(PCA_metrics$x, aes(x = PC1, y = PC2, label = rownames(PCA_metrics$x))) +
  geom_point(color = "black", size = 1.0) +
  geom_text_repel(box.padding = 0.5,
                  point.padding = 0.3,
                  min.segment.length = 0,
                  segment.color = "grey50",
                  segment.alpha = 0.8,
                  max.overlaps = Inf,
                  size = 3) +
  theme_bw() +
  labs(x = paste0("Principal Component 1 (", percent_explained[1], '% variance explained)'),
       y = paste0("Principal Component 2 (", percent_explained[2], '% variance explained)')) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

combined_plot <- plot_grid(grid.grabExpr(draw(raw_heatmap, padding = unit(c(2, 20, 2, 2), "mm"))),
          PC1_vs_PC2,
          labels = c('a', 'b'), rel_widths = c(1, 0.8))

ggsave(plot = combined_plot,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Main_metric_similarity.pdf',
       dpi = 600,
       width = 14,
       height = 6)
