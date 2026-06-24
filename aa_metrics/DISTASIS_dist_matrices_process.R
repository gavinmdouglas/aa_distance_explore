rm(list = ls(all.names = TRUE))

# DISTASIS - combined distance matrices to create, based on PAML analysis on first 10 alignments,
# and relationship of metrics in PCA based on 1 - RVcoeff.

# These are the focal metrics to explore:
# ○ Custom median
# ○ EX
# ○ VHSE
# ○ Miyata
# ○ Other (AA-Ont)
# 
# Every individual with custom (median approach), along with growing list in that order
# (Also a growing list without EX too!)
# 
# So, 9 combinations: Custom +…
# EX
# VHSE
# Miyata
# Other (AA-Ont)
# EX + VHSE
# EX + V + M
# EX + V + M + O
# VHSE + M
# VHSE + M + O

# Also decided to try combinations with EX and others, excluding custom median, just for comparison.

# And decided to include DeMaSk too for an additional subset:
# EX + DeMaSk
# Custom median + DeMaSk
# Custom median + EX+ DeMaSk

library('DistatisR')

check_symmetrical <- function(mat) {
  if (! all(mat == t(mat))) {
    FALSE
  } else {
    TRUE
  }
}

aa_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/aa_map.tsv.gz',
                     sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 3)

metric_id_map <- read.table('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                            header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names=1)

focal_sim_files <- c("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/proteinGym_rsa_median_custom.tsv.gz",
                     "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/ex.tsv.gz",
                     "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/VHSE.tsv.gz",
                     "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/miyata_orig.tsv.gz",
                     "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/AAOntology-PCs_Others.tsv.gz",
                     "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/demask.tsv.gz")

aas <- sort(rownames(aa_map))

dist_dummy <- data.frame(matrix(NA, nrow=20, ncol=20))
rownames(dist_dummy) <- aas
colnames(dist_dummy) <- aas

dist_tabs <- list()
for (simfile in focal_sim_files) {
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
  
  # Check if symmetrical.
  if (!check_symmetrical(dist_tabs[[metric_name]])) {
    message(paste0("Metric ", metric_name, " is not symmetrical. Will make symmetrical now."))
    dist_tabs[[metric_name]] <- (dist_tabs[[metric_name]] + t(dist_tabs[[metric_name]])) / 2
  }
}

names(dist_tabs) <- sub('-', '.', names(dist_tabs))
names(dist_tabs) <- sub('.prob', '', names(dist_tabs))
names(dist_tabs) <- sub('_prob', '', names(dist_tabs))
names(dist_tabs) <- sub('_maxscaled', '', names(dist_tabs))
names(dist_tabs) <- metric_id_map[names(dist_tabs), 'Clean']
names(dist_tabs) <- gsub(' \\(orig.\\)', '', names(dist_tabs))
names(dist_tabs)[which(names(dist_tabs) == "DMS-EX")] <- "Custom"
names(dist_tabs)[which(names(dist_tabs) == "Other (AA-Ont.)")] <- "OtherAAOnt"

focal_combinations <- list(
  c("Custom", "EX"),
  c("Custom", "VHSE"),
  c("Custom", "Miyata"),
  c("Custom", "OtherAAOnt"),
  c("Custom", "EX", "VHSE"),
  c("Custom", "EX", "VHSE", "Miyata"),
  c("Custom", "EX", "VHSE", "Miyata", "OtherAAOnt"),
  c("Custom", "VHSE", "Miyata"),
  c("Custom", "VHSE", "Miyata", "OtherAAOnt"),
  c('EX', 'VHSE'),
  c('EX', 'VHSE', 'Miyata'),
  c('EX', 'VHSE', 'Miyata', "OtherAAOnt"),
  c('Custom', 'DeMaSk'),
  c('EX', 'DeMaSk'),
  c('Custom', 'EX', 'DeMaSk')
)

transform_to_range <- function(numbers, min_val = 0.01, max_val = 0.99) {
  input_min <- min(numbers)
  input_max <- max(numbers)
  return(((numbers - input_min) / (input_max - input_min)) * (max_val - min_val) + min_val)
}

for (focal_combo_i in seq_along(focal_combinations)) {
  focal_combo <- focal_combinations[[focal_combo_i]]
  
  if (! all(focal_combo %in% names(dist_tabs))) {
    stop(paste("Missing metrics:", paste(focal_combo[!focal_combo %in% names(dist_tabs)], collapse=", ")))
  }

  focal_combo_name <- paste(focal_combo, collapse = "_")
  
  all_dist_matrices <- array(0, dim = c(20, 20, length(focal_combo)))
  
  for (i in seq_along(focal_combo)) {
    all_dist_matrices[, , i] <- as.matrix(dist_tabs[[focal_combo[i]]])
  }
  
  distasis_output <- DistatisR::distatis(all_dist_matrices, nfact2keep = 20)
  
  eigenvalues <- distasis_output$res4Splus$eigValues
  
  write.table(x = eigenvalues,
              file = paste0("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/DISTATIS_working/eigenvalues/", focal_combo_name, ".txt"),
              quote = FALSE, col.names = FALSE, row.names = FALSE)

  factor_scores <- data.frame(distasis_output$res4Splus$F)
  orig_col <- colnames(factor_scores)
  factor_scores$amino_acids <- aas
  factor_scores <- factor_scores[, c("amino_acids", orig_col)]
  write.table(x = factor_scores,
              file = paste0("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/DISTATIS_working/factor_scores/", focal_combo_name, ".tsv"),
              quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  
  # Keep factors for positive eigenvalues.
  pos_eigenval_i <- which(eigenvalues > 0)
  eigenvalues <- eigenvalues[pos_eigenval_i]
  
  factors <- distasis_output$res4Splus$F[, pos_eigenval_i]
  
  # Initially thought I should normalized for eigenvalues, but then realized
  # DISTATIS already did this implicitly, so the below commands commented
  # out were *not* used.
  #percent_explained <- (eigenvalues / sum(eigenvalues)) * 100
  #weighted_factors <- sweep(factors, 2, sqrt(percent_explained / 100), "*")
  #pairwise_dist <- data.frame(as.matrix(dist(weighted_factors)))
  
  pairwise_dist <- data.frame(as.matrix(dist(factors)))
  colnames(pairwise_dist) <- aas
  pairwise_dist$aa1 <- aas
  
  out_tab <- reshape2::melt(pairwise_dist, id.vars = "aa1")
  colnames(out_tab) <- c('aa1', 'aa2', 'similarity')
  
  out_tab <- out_tab[out_tab$aa1 != out_tab$aa2, ]
  
  out_tab$similarity <- 1 - transform_to_range(out_tab$similarity)
  
  write.table(x = out_tab,
              file = gzfile(paste0("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/", "DISTATIS_", focal_combo_name, ".tsv.gz")),
              quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

}
