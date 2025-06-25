rm(list = ls(all.names = TRUE))

# Explore correspondance between EX exchangeability and distance from other metrics.
# Are there obvious outliers in EX that are not found by other metrics?
# Is there a best way to re-scale AA distances to get them closer to EX exchangeabilities?

library(ggplot2)

aa_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/aa_map.tsv.gz',
                     sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 3)
aa_map$long_w_letter <- stringr::str_to_title(paste0(aa_map$full, ' (', rownames(aa_map), ')'))
aa_map$long_w_letter <- sub(' Acid ', ' acid ', aa_map$long_w_letter)

metrics_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  

metric_sim_files <- list.files('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity/',
                              pattern = '', full.names = TRUE)

aas <- sort(rownames(aa_map))

aa1_set <- character()
aa2_set <- character()
for (aa1 in aas[-length(aas)]) {
  for (aa2 in aas[(which(aas == aa1) + 1 ):length(aas)]) {
    aa1_set <- c(aa1_set, aa1)
    aa2_set <- c(aa2_set, aa2)
  }
}

metric_sim <- list()

for (metric_sim_file in metric_sim_files) {
  metric_name <- sub('.tsv.gz', '', basename(metric_sim_file))
  raw_tab <- read.table(metric_sim_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

  metric_vec <- numeric()
  
  for (aa1 in aas[-length(aas)]) {
    for (aa2 in aas[(which(aas == aa1) + 1 ):length(aas)]) {
        metric_vec <- c(metric_vec, raw_tab[which(raw_tab$aa1 == aa1 & raw_tab$aa2 == aa2), 3])
    }
  }
  
  metric_sim[[metric_name]] <- metric_vec
}

plot(metric_sim$ex.maxscaled, metric_sim$demask)


plot(metric_sim$ex.maxscaled, metric_sim$zScales.maxscaled.invert)

plot(metric_sim$ex.maxscaled, metric_sim$`AAOntology-PCs_Others.maxscaled.invert`)


focal_metrics <- data.frame('VHSE'=metric_sim$VHSE.maxscaled.invert,
                            'zScales'=metric_sim$zScales.maxscaled.invert,
                            'Miyata'=metric_sim$miyata_orig.maxscaled.invert,
                            'Other_AAOnt'=metric_sim$`AAOntology-PCs_Others.maxscaled.invert`,
                            'EMPAR'=metric_sim$empar.maxscaled,
                            'Sneath'=metric_sim$sneath.maxscaled)

focal_metrics_w_EX <- data.frame('VHSE'=metric_sim$VHSE.maxscaled.invert,
                            'zScales'=metric_sim$zScales.maxscaled.invert,
                            'Miyata'=metric_sim$miyata_orig.maxscaled.invert,
                            'Other_AAOnt'=metric_sim$`AAOntology-PCs_Others.maxscaled.invert`,
                            'EMPAR'=metric_sim$empar.maxscaled,
                            'Sneath'=metric_sim$sneath.maxscaled,
                            'EX'=metric_sim$ex.maxscaled)

mean_sim_focal_metrics <- rowMeans(focal_metrics)
mean_sim_focal_metrics_w_EX <- rowMeans(focal_metrics_w_EX)



focal_metrics_scaled <- data.frame(scale(focal_metrics, center = TRUE, scale = TRUE))
colnames(focal_metrics_scaled) <- colnames(focal_metrics)

ex_scaled <- scale(metric_sim$ex.maxscaled, center = TRUE, scale = TRUE)

focal_metrics_scaled_diff <- focal_metrics_scaled
for (metric in colnames(focal_metrics_scaled)) {
  focal_metrics_scaled_diff[, metric] <- abs(focal_metrics_scaled[, metric] - ex_scaled)
}

focal_metrics_diff <- data.frame(focal_metrics)
for (metric in colnames(focal_metrics)) {
  focal_metrics_diff[, metric] <- abs(focal_metrics[, metric] - metric_sim$ex.maxscaled)
}
