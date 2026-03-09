rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(ggtext)

by_aa_scaled <- read.table('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/codon_exchangeabilities_prepped_scaled.tsv.gz',
                           sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

by_aa_scaled <- by_aa_scaled[-which(rownames(by_aa_scaled) == "DEX"), ]
   
clean_to_group <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                             sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 2)

row_labels <- rownames(by_aa_scaled)

row_group <- clean_to_group[rownames(by_aa_scaled), 'Grouping']

group_colours <- c("Physicochemical/structural" = "black",
                   "Experimental"              = "#0072B2",
                   #"Deep learning"             = "#CC79A7",
                   "Substitution-based"        = "#009E73",
                   "Radical vs. conservative"  = "#D55E00")
label_colours <- group_colours[row_group]

raw_heatmap <- ComplexHeatmap::Heatmap(matrix = as.matrix(by_aa_scaled),
                                       cluster_rows = TRUE,
                                       cluster_columns = TRUE,
                                       clustering_method_rows = 'complete',
                                       clustering_method_columns = 'complete',
                                       clustering_distance_rows = 'euclidean',
                                       clustering_distance_columns = 'euclidean',
                                       column_names_rot = 45,
                                       row_labels = row_labels,
                                       row_names_gp = gpar(col = label_colours),
                                       name = 'Scaled mean\nexchangeabilities')


pca_workflow <- function(in_scaled_matrix) {
 
  PCA_metrics <- prcomp(in_scaled_matrix, center = FALSE, scale. = FALSE)
  eigenvalues <- PCA_metrics$sdev ** 2
  percent_explained <- sprintf('%.1f', eigenvalues / sum(eigenvalues) * 100)
  
  pca_df <- as.data.frame(PCA_metrics$x)
  pca_df$label <- rownames(PCA_metrics$x)
  pca_df$group <- clean_to_group[rownames(pca_df), 'Grouping']

  PC1_vs_PC2 <- ggplot(pca_df, aes(x = PC1, y = PC2, label = label, colour=group)) +
    geom_point(size = 1.0) +
    geom_text_repel(box.padding = 0.5,
                    point.padding = 0.2,
                    min.segment.length = 0,
                    force = 1,
                    force_pull = 0.1,
                    segment.color = "grey90",
                    segment.alpha = 0.8,
                    max.overlaps = Inf,
                    size = 3) +
    theme_bw() +
    scale_colour_manual(values = group_colours, guide = "none") +
    labs(x = paste0("PC1 (", percent_explained[1], '%)'),
         y = paste0("PC2 (", percent_explained[2], '%)')) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())

  loadings_df <- as.data.frame(PCA_metrics$rotation[, 1:2])
  loadings_df$aa <- gsub(".*\\((.*)\\)", "\\1", rownames(loadings_df))
  
  loadings_plot <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = aa)) +
                          geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                                       arrow = arrow(length = unit(0.2, "cm")), colour = "steelblue") +
                          geom_text_repel(size = 3, box.padding = 0.4, max.overlaps = Inf, segment.color = NA) +
                        theme_bw() +
                        labs(x = paste0("PC1 loadings (", percent_explained[1], "%)"),
                             y = paste0("PC2 loadings (", percent_explained[2], "%)")) +
                        theme(panel.grid = element_blank(),
                              plot.margin = margin(1, 1, 1, 1, "cm"))
  
  
  
  return(list(PCs=PC1_vs_PC2, loadings=loadings_plot))
}

# First, calculate full PCA with all measures tested.
PC1_vs_PC2_full_out <- pca_workflow(by_aa_scaled)

# Decided to remove these too as they are major outliers that dominate the PCA, despite being minor metrics:
by_aa_scaled <- by_aa_scaled[! rownames(by_aa_scaled) %in% c("CSW", "EMPAR", "Conf. (AA-Ont.)", "Comp. (AA-Ont.)",
                                                             "RvC (charge)", "RvC (Zhang)", "Venkatarajan"), ]

PC1_vs_PC2_subset_out <- pca_workflow(by_aa_scaled)
PC1_vs_PC2_subset <- PC1_vs_PC2_subset_out$PCs
PC1_vs_PC2_subset_loadings <- PC1_vs_PC2_subset_out$loadings

group_legend <- Legend(
  labels = names(group_colours),
  title = "Measure grouping",
  type = "points",
  pch = NA,
  background = "white",
  border = NA,
  grid_width = unit(0, "mm"),
  labels_gp = gpar(col = group_colours)
)

draw(raw_heatmap, annotation_legend_list = list(group_legend))

heatmap_row <- plot_grid(grid.grabExpr(draw(raw_heatmap,
                                            padding = unit(c(2, 20, 2, 2), "mm"),
                                            merge_legend = TRUE,
                                            annotation_legend_list = list(group_legend),
                                            heatmap_legend_side = "right",
                                            annotation_legend_side = "right")),
                         labels = 'a')
pca_row <- plot_grid(PC1_vs_PC2_subset, PC1_vs_PC2_subset_loadings,
                     labels = c('b', 'c'))
combined_plot <- plot_grid(heatmap_row, pca_row, ncol = 1, rel_heights = c(1.5, 1))
combined_plot
ggsave(plot = combined_plot,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Main_measure_similarity.pdf',
       dpi = 600,
       width = 11.5,
       height = 13)

full_PCA_plot <- plot_grid(PC1_vs_PC2_full_out$PCs, PC1_vs_PC2_full_out$loadings, labels = c('a', 'b'))
ggsave(plot = full_PCA_plot,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Supp_measure_PCA_full.pdf',
       dpi = 600,
       width = 14,
       height = 7)
