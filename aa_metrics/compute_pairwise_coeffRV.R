rm(list = ls(all.names = TRUE))

library(FactoMineR)
library(ggrepel)
library(ggplot2)

aa_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/aa_map.tsv.gz',
                     sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 3)

aa_sim_files <- list.files('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/',
                           pattern = 'tsv.gz', full.names = TRUE)

aa_sim_files <- grep("test_", aa_sim_files, invert = TRUE, value = TRUE)

specific_str_to_keep <- c("grantham_orig.tsv.gz",
                          "miyata_orig.tsv.gz",
                          "proteinGym_rsa_median_custom.tsv.gz",
                          "proteinGym_rsa_mean_custom.tsv.gz",
                          "proteinGym_mean_standard_score_rsa_weighted.tsv.gz",
                          "proteinGym_median_robust_score_rsa_weighted.tsv.gz")

specific_file_to_keep <- character()
for (str in specific_str_to_keep) {
  specific_file_to_keep <- c(specific_file_to_keep, grep(str, aa_sim_files, value = TRUE))
}

aa_sim_files <- grep("grantham", aa_sim_files, invert = TRUE, value = TRUE)
aa_sim_files <- grep("miyata", aa_sim_files, invert = TRUE, value = TRUE)
aa_sim_files <- grep("proteinGym", aa_sim_files, invert = TRUE, value = TRUE)

aa_sim_files <- c(aa_sim_files, specific_file_to_keep)

metrics_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

aas <- sort(rownames(aa_map))

dist_dummy <- data.frame(matrix(NA, nrow=20, ncol=20))
rownames(dist_dummy) <- aas
colnames(dist_dummy) <- aas

dist_tabs <- list()
for (simfile in aa_sim_files) {
  tab <- read.table(simfile, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  metric_name <- gsub('.tsv.gz', '', basename(simfile))
  dist_tabs[[metric_name]] <- dist_dummy
  for (aa1 in aas) {
    for (aa2 in aas) {
      if (aa1 == aa2) {
        dist_tabs[[metric_name]][aa1, aa2] <- 0.0
        next 
      }
      dist_tabs[[metric_name]][aa1, aa2] <- 1.0 - tab[which(tab$aa1 == aa1 & tab$aa2 == aa2), 3]
    }
  }
  
  # Throw error if any NA's.
  if (any(is.na(dist_tabs[[metric_name]]))) {
    stop(paste0("NA values found in metric: ", metric_name, ". Please check the input data."))
  }
}

metric_id_map <- read.table('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                            header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names=1)

names(dist_tabs) <- sub('-', '.', names(dist_tabs))
names(dist_tabs) <- sub('.prob', '', names(dist_tabs))
names(dist_tabs) <- sub('_prob', '', names(dist_tabs))
names(dist_tabs) <- sub('_maxscaled', '', names(dist_tabs))

names(dist_tabs) <- metric_id_map[names(dist_tabs), 'Clean']
names(dist_tabs) <- gsub(' \\(orig.\\)', '', names(dist_tabs))


metric_coeffRV <- data.frame(matrix(NA, nrow = length(dist_tabs), ncol = length(dist_tabs)))
rownames(metric_coeffRV) <- names(dist_tabs)
colnames(metric_coeffRV) <- names(dist_tabs)

for (i in seq_along(dist_tabs)) {
  for (j in seq_along(dist_tabs)) {
    if (i == j) {
      metric_coeffRV[i, j] <- 1.0
    } else {
      coeffRV <- coeffRV(dist_tabs[[i]], dist_tabs[[j]])
      metric_coeffRV[i, j] <- coeffRV$rv
    }
  }
}

metric_coeffRV_dist <- as.dist(1 - metric_coeffRV)

PCA_metrics <- prcomp(metric_coeffRV_dist, center = TRUE, scale. = TRUE)
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


ggsave(plot = PC1_vs_PC2,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Additional_custom_metric_similarity.pdf',
       dpi = 600,
       width = 6,
       height = 6)


# Also make the PCA without the extra custom metrics:
#"proteinGym_rsa_mean_custom.tsv.gz",
#"proteinGym_mean_standard_score_rsa_weighted.tsv.gz",
#"proteinGym_median_robust_score_rsa_weighted.tsv.gz"

metric_coeffRV_trim <- metric_coeffRV[!grepl("coarse", rownames(metric_coeffRV)), 
                                      !grepl("coarse", colnames(metric_coeffRV))]

metric_coeffRV_trim <- metric_coeffRV_trim[!grepl("This study \\(mean\\)", rownames(metric_coeffRV_trim)), 
                                           !grepl("This study \\(mean\\)", colnames(metric_coeffRV_trim))]

metric_coeffRV_dist_trim <- as.dist(1 - metric_coeffRV_trim)
PCA_metrics_trim <- prcomp(metric_coeffRV_dist_trim, center = TRUE, scale. = TRUE)

eigenvalues_trim <- PCA_metrics_trim$sdev ** 2
percent_explained_trim <- sprintf('%.1f', eigenvalues_trim / sum(eigenvalues_trim) * 100)

PC1_vs_PC2_trim <- ggplot(PCA_metrics_trim$x, aes(x = PC1, y = PC2, label = rownames(PCA_metrics_trim$x))) +
  geom_point(color = "black", size = 1.0) +
  geom_text_repel(box.padding = 0.5,
                  point.padding = 0.3,
                  min.segment.length = 0,
                  segment.color = "grey50",
                  segment.alpha = 0.8,
                  max.overlaps = Inf,
                  size = 3) +
  theme_bw() +
  labs(x = paste0("Principal Component 1 (", percent_explained_trim[1], '% variance explained)'),
       y = paste0("Principal Component 2 (", percent_explained_trim[2], '% variance explained)')) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

ggsave(plot = PC1_vs_PC2_trim,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Main_metric_similarity.pdf',
       dpi = 600,
       width = 6,
       height = 6)