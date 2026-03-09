# Run PCA
PCA_metrics <- prcomp(by_aa_scaled, center = TRUE, scale. = TRUE)
eigenvalues <- PCA_metrics$sdev ** 2
percent_explained <- sprintf('%.1f', eigenvalues / sum(eigenvalues) * 100)

PC1_vs_PC2 <- ggplot(PCA_metrics$x, aes(x = PC1, y = PC2, label = rownames(PCA_metrics$x))) +
  geom_point(color = "grey80", size = 0.5) +
  geom_text_repel(box.padding = 0.5,
                  point.padding = 0.3,
                  min.segment.length = 0,
                  segment.color = "grey50",
                  segment.alpha = 0.8,
                  max.overlaps = Inf,
                  size = 3) +
  theme_bw() +
  labs(x = paste0("Principal Component 1 (", percent_explained[1], '% variance explained)'),
       y = paste0("Principal Component 2 (", percent_explained[2], '% variance explained)'))

prop_explained <- eigenvalues / sum(eigenvalues)
weighted_PCA_dist <- matrix(NA, nrow(by_aa_scaled), nrow(by_aa_scaled))
for (i in 1:nrow(by_aa_scaled)) {
  for (j in 1:nrow(by_aa_scaled)) {
    weighted_PCA_dist[i, j] <- sqrt(sum((PCA_metrics$x[i, ] - PCA_metrics$x[j, ])^2 * prop_explained))
  }
}

MDS_metrics_PC_space <- cmdscale(weighted_PCA_dist, k = nrow(PCA_metrics$x) - 1, eig = TRUE)

MDS_metrics_PC_space_points <- data.frame(MDS_metrics_PC_space$points)
colnames(MDS_metrics_PC_space_points) <- paste0('PC', 1:ncol(MDS_metrics_PC_space_points))

PC1_vs_PC2 <- ggplot(MDS_metrics_PC_space_points, aes(x = PC1, y = PC2, label = rownames(PCA_metrics$x))) +
  geom_point(color = "grey80", size = 1) +
  geom_text_repel(box.padding = 0.5,
                  point.padding = 0.3,
                  min.segment.length = 0,
                  segment.color = "grey50",
                  segment.alpha = 0.8,
                  max.overlaps = Inf,
                  size = 3) +
  theme_bw() +
  labs(x = paste0("Principal Component 1 (", percent_explained[1], '% variance explained)'),
       y = paste0("Principal Component 2 (", percent_explained[2], '% variance explained)'))


library(Rtsne)

set.seed(42)

tsne_result <- Rtsne(X = by_aa_scaled,
                     dims = 2,
                     perplexity = 5,
                     theta = 0,
                     check_duplicates = TRUE,
                     pca = TRUE,
                     max_iter = 100000,
                     verbose = FALSE)

tsne_df <- data.frame(tSNE1 = tsne_result$Y[, 1], tSNE2 = tsne_result$Y[, 2])

ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, label = rownames(PCA_metrics$x))) +
  geom_point(color = "grey80", size = 0.5) +
  geom_text_repel(box.padding = 0.5,
                  point.padding = 0.3,
                  min.segment.length = 0,
                  segment.color = "grey50",
                  segment.alpha = 0.8,
                  max.overlaps = Inf,
                  size = 3) +
  theme_bw() +
  labs(x = 'tSNE Dimension 1', y = 'tSNE Dimension 2')



prop_explained <- eigenvalues / sum(eigenvalues)
weighted_PCA_dist <- matrix(NA, nrow(by_aa_scaled), nrow(by_aa_scaled))
for (i in 1:nrow(by_aa_scaled)) {
  for (j in 1:nrow(by_aa_scaled)) {
    weighted_PCA_dist[i, j] <- sqrt(sum((PCA_metrics$x[i, ] - PCA_metrics$x[j, ])^2 * prop_explained))
  }
}

Heatmap(weighted_PCA_dist, rect_gp = gpar(type = "none"),
        show_column_dend = FALSE,
        row_labels = rownames(by_aa_scaled),
        column_labels = rownames(by_aa_scaled),
        column_names_rot = 45,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(as.integer(x) <= 1 - as.integer(y)) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        })