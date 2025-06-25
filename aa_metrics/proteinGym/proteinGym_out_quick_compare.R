rm(list = ls(all.names = TRUE))

ex <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/ex.tsv.gz',
                 header=TRUE, sep = '\t', stringsAsFactors = FALSE)

median_custom <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/proteinGym_rsa_median_custom.tsv.gz',
                          header=TRUE, sep = '\t', stringsAsFactors = FALSE)

mean_custom <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/proteinGym_rsa_mean_custom.tsv.gz',
                          header=TRUE, sep = '\t', stringsAsFactors = FALSE)

earlier_mean <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/proteinGym_mean_standard_score_rsa_weighted.tsv.gz',
                           header=TRUE, sep = '\t', stringsAsFactors = FALSE)

beta_median <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/proteinGym_rsa_median_ex2005_style_beta.tsv.gz',
                 header=TRUE, sep = '\t', stringsAsFactors = FALSE)

beta_mean <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/proteinGym_rsa_mean_ex2005_style_beta.tsv.gz',
                          header=TRUE, sep = '\t', stringsAsFactors = FALSE)

rsa_latest <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/proteinGym_rsa_median_ex2005_style_nontrimmed.tsv.gz',
                            header=TRUE, sep = '\t', stringsAsFactors = FALSE)

ex_method_out <- read.table("~/test.tsv",
                            header=TRUE, sep = '\t', stringsAsFactors = FALSE)

proteingym_scores <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_processed.tsv',
                                header=TRUE, sep = '\t', stringsAsFactors = FALSE)

proteingym_scores_rsa <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_processed_rsa_weighted.tsv',
                                    header=TRUE, sep = '\t', stringsAsFactors = FALSE)

proteingym_scores_directed <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_processed_directed.tsv',
                                         header=TRUE, sep = '\t', stringsAsFactors = FALSE)

ex <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/ex.tsv.gz',
                 header=TRUE, sep = '\t', stringsAsFactors = FALSE)

proteingym_scores$weighted_median_robust_score <- NA
proteingym_scores$weighted_median_robust_score_min <- NA
proteingym_scores$weighted_median_robust_score_max <- NA

proteingym_scores$weighted_mean_standard_score <- NA
proteingym_scores$weighted_mean_standard_score_min <- NA
proteingym_scores$weighted_mean_standard_score_max <- NA

proteingym_scores$mean_standard_score_directed <- NA
proteingym_scores$median_robust_score_directed <- NA
proteingym_scores$mean_standard_score_directed_min <- NA
proteingym_scores$median_robust_score_directed_min <- NA
proteingym_scores$mean_standard_score_directed_max <- NA
proteingym_scores$median_robust_score_directed_max <- NA

proteingym_scores$ex <- NA

proteingym_scores$ex_style <- NA

for (i in 1:nrow(proteingym_scores)) {
  aa1 <- proteingym_scores$aa1[i]
  aa2 <- proteingym_scores$aa2[i]

  proteingym_scores[i, 'ex'] <- ex[ex$aa1 == aa1 & ex$aa2 == aa2, 'similarity']
  
  rsa_standard_scores <- c(proteingym_scores_rsa[proteingym_scores_rsa$ref_aa == aa1 & proteingym_scores_rsa$mut_aa == aa2, 'weighted_median_standard_score'],
                           proteingym_scores_rsa[proteingym_scores_rsa$ref_aa == aa2 & proteingym_scores_rsa$mut_aa == aa1, 'weighted_median_standard_score'])
  
  rsa_robust_scores <- c(proteingym_scores_rsa[proteingym_scores_rsa$ref_aa == aa1 & proteingym_scores_rsa$mut_aa == aa2, 'weighted_median_robust_score'],
                         proteingym_scores_rsa[proteingym_scores_rsa$ref_aa == aa2 & proteingym_scores_rsa$mut_aa == aa1, 'weighted_median_robust_score'])
  
  directed_standard_scores <- c(proteingym_scores_directed[proteingym_scores_directed$aa1 == aa1 & proteingym_scores_directed$aa2 == aa2, 'mean_standard_score'],
                                 proteingym_scores_directed[proteingym_scores_directed$aa1 == aa2 & proteingym_scores_directed$aa2 == aa1, 'mean_standard_score'])
  
  directed_robust_scores <- c(proteingym_scores_directed[proteingym_scores_directed$aa1 == aa1 & proteingym_scores_directed$aa2 == aa2, 'median_robust_score'],
                               proteingym_scores_directed[proteingym_scores_directed$aa1 == aa2 & proteingym_scores_directed$aa2 == aa1, 'median_robust_score'])
  
  proteingym_scores[i, 'weighted_median_robust_score'] <- mean(rsa_robust_scores)
  proteingym_scores[i, 'weighted_median_robust_score_min'] <- min(rsa_robust_scores)
  proteingym_scores[i, 'weighted_median_robust_score_max'] <- max(rsa_robust_scores)
  
  proteingym_scores[i, 'weighted_mean_standard_score'] <- mean(rsa_standard_scores)
  proteingym_scores[i, 'weighted_mean_standard_score_min'] <- min(rsa_standard_scores)
  proteingym_scores[i, 'weighted_mean_standard_score_max'] <- max(rsa_standard_scores)
  
  proteingym_scores[i, 'mean_standard_score_directed'] <- mean(directed_standard_scores)
  proteingym_scores[i, 'median_robust_score_directed'] <- mean(directed_robust_scores)
  proteingym_scores[i, 'mean_standard_score_directed_min'] <- min(directed_standard_scores)
  proteingym_scores[i, 'median_robust_score_directed_min'] <- min(directed_robust_scores)
  proteingym_scores[i, 'mean_standard_score_directed_max'] <- max(directed_standard_scores)
  proteingym_scores[i, 'median_robust_score_directed_max'] <- max(directed_robust_scores)
  
  exstyle_aa1_val <- ex_method_out[ex_method_out$ref_aa == aa1 & ex_method_out$mut_aa == aa2, 'weighted_mean']
  exstyle_aa1_N <- ex_method_out[ex_method_out$ref_aa == aa1 & ex_method_out$mut_aa == aa2, 'N_exposed'] + ex_method_out[ex_method_out$ref_aa == aa1 & ex_method_out$mut_aa == aa2, 'N_buried']
  
  exstyle_aa2_val <- ex_method_out[ex_method_out$ref_aa == aa2 & ex_method_out$mut_aa == aa1, 'weighted_mean']
  exstyle_aa2_N <- ex_method_out[ex_method_out$ref_aa == aa2 & ex_method_out$mut_aa == aa1, 'N_exposed'] + ex_method_out[ex_method_out$ref_aa == aa2 & ex_method_out$mut_aa == aa1, 'N_buried']
  
  ex_style_weighted_mean <- (exstyle_aa1_val * exstyle_aa1_N + exstyle_aa2_val * exstyle_aa2_N) / (exstyle_aa1_N + exstyle_aa2_N)
  proteingym_scores$ex_style[i] <- ex_style_weighted_mean
  
}

proteingym_scores <- proteingym_scores[, -which(colnames(proteingym_scores) == "num_genes_w_obs")]

aa_combos <- proteingym_scores[, c(1, 2)]

proteingym_scores <- proteingym_scores[, -c(1, 2)]

ex_vec <- proteingym_scores$ex

proteingym_scores <- proteingym_scores[, -which(colnames(proteingym_scores) == 'ex')]

spearman_r <- numeric()
for (metric in colnames(proteingym_scores)) {
  spearman_r <- c(spearman_r, cor.test(proteingym_scores[[metric]], ex_vec, method = 'spearman', exact = FALSE)$estimate)
}

names(spearman_r) <- colnames(proteingym_scores)

proteingym_scores_subset <- proteingym_scores

proteingym_scores_subset <- proteingym_scores_subset[, c("weighted_median_robust_score_min",
                                                         "weighted_median_robust_score_max",
                                                         "weighted_mean_standard_score_min",
                                                         "weighted_mean_standard_score_max",
                                                         "mean_standard_score_directed",
                                                         "median_robust_score_directed",
                                                         "mean_standard_score_directed_min",
                                                         "median_robust_score_directed_min",
                                                         "mean_standard_score_directed_max",
                                                         "median_robust_score_directed_max")]


unique_aas <- sort(unique(c(aa_combos$aa1, aa_combos$aa2)))

for (metric in colnames(proteingym_scores_subset)) {

  outfile <- paste0('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/test_proteingym_', metric, '.tsv.gz')
  
  tab <- data.frame(matrix(0, nrow=20, ncol=20))
  rownames(tab) <- unique_aas
  colnames(tab) <- unique_aas
  
  for (i in 1:nrow(proteingym_scores_subset)) {
    aa1 <- aa_combos$aa1[i]
    aa2 <- aa_combos$aa2[i]
    tab[aa1, aa2] <- proteingym_scores_subset[i, metric]
    tab[aa2, aa1] <- proteingym_scores_subset[i, metric]
  }

  write.table(x = tab, file = gzfile(outfile),
              sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)
}
