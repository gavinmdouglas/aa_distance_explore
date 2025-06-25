rm(list = ls(all.names = TRUE))

unique_aas <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

prep_proteingym_sim_tab_w_rsa <- function(infile, focal_col, N_col1, N_col2) {
  
  intab <- read.table(infile, header=TRUE, sep='\t', stringsAsFactors = FALSE)
  
  prepped <- data.frame(matrix(0, nrow=20, ncol=20))
  rownames(prepped) <- unique_aas
  colnames(prepped) <- unique_aas
  
  for (aa1 in unique_aas[-length(unique_aas)]) {

    aa1_index <- which(unique_aas == aa1)

    for (aa2 in unique_aas[(aa1_index + 1):length(unique_aas)]) {

      forward_sub_row <- intab[intab$ref_aa == aa1 & intab$mut_aa == aa2, ]
      backward_sub_row <- intab[intab$ref_aa == aa2 & intab$mut_aa == aa1, ]
      
      forward_N <- forward_sub_row[, N_col1] + forward_sub_row[, N_col2]
      backward_N <- backward_sub_row[, N_col1] + backward_sub_row[, N_col2]
      
      weighted_mean_val <- ((forward_sub_row[, focal_col] * forward_N) +
                            (backward_sub_row[, focal_col] * backward_N)) / (forward_N + backward_N)

      prepped[aa1, aa2] <- weighted_mean_val
      prepped[aa2, aa1] <- weighted_mean_val

    }
  }
  
  return(prepped)
}

weighted_mean_ex2005_style <- prep_proteingym_sim_tab_w_rsa('~/beta_test.tsv', 'weighted_mean', N_col1="N_buried", N_col2="N_exposed")
weighted_median_ex2005_style <- prep_proteingym_sim_tab_w_rsa('~/beta_test.tsv', 'weighted_median', N_col1="N_buried", N_col2="N_exposed")

write.table(weighted_mean_ex2005_style,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_rsa_mean_ex2005_style_beta.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)

write.table(weighted_median_ex2005_style,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_rsa_median_ex2005_style_beta.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)

weighted_mean_custom <- prep_proteingym_sim_tab_w_rsa('~/test.tsv', 'weighted_mean_standard_score', N_col1="N_exposed_mean_standard_score", N_col2="N_buried_mean_standard_score")

weighted_median_custom <- prep_proteingym_sim_tab_w_rsa('~/test.tsv', 'weighted_median_robust_score', N_col1="N_exposed_median_robust_score", N_col2="N_buried_median_robust_score")

write.table(weighted_mean_custom,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_rsa_mean_custom.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)

write.table(weighted_median_custom,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_rsa_median_custom.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)

proteingym_scores <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_processed.tsv',
                                header=TRUE, sep = '\t', stringsAsFactors = FALSE)

proteingym_ex2005 <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/',
                                    header=TRUE, sep = '\t', stringsAsFactors = FALSE)

# And do the same, but for RSA category weighted values too.
# Note that there are different values per ref and mut AAs, so the mean is taken.
proteingym_scores_rsa <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_processed_rsa_weighted.tsv',
                                    header=TRUE, sep = '\t', stringsAsFactors = FALSE)

median_tab_rsa <- data.frame(matrix(0, nrow=20, ncol=20))
mean_tab_rsa <- data.frame(matrix(0, nrow=20, ncol=20))

median_score_rsa <-  data.frame(matrix(0, nrow=20, ncol=20))

unique_aas_rsa <- sort(unique(c(proteingym_scores_rsa$ref_aa, proteingym_scores_rsa$mut_aa)))

rownames(median_tab_rsa) <- unique_aas_rsa
rownames(mean_tab_rsa) <- unique_aas_rsa
colnames(median_tab_rsa) <- unique_aas_rsa
colnames(mean_tab_rsa) <- unique_aas_rsa
rownames(median_score_rsa) <- unique_aas_rsa
colnames(median_score_rsa) <- unique_aas_rsa

# Also for the directional overall score computations.
proteingym_scores_directed <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_processed_directed.tsv',
                                         header=TRUE, sep = '\t', stringsAsFactors = FALSE)

mean_tab_directed <- data.frame(matrix(0, nrow=20, ncol=20))
median_tab_directed <- data.frame(matrix(0, nrow=20, ncol=20))

mean_tab_directed_min <- data.frame(matrix(0, nrow=20, ncol=20))
median_tab_directed_min <- data.frame(matrix(0, nrow=20, ncol=20))

mean_tab_directed_max <- data.frame(matrix(0, nrow=20, ncol=20))
median_tab_directed_max <- data.frame(matrix(0, nrow=20, ncol=20))

unique_aas_directed <- sort(unique(c(proteingym_scores_directed$aa1, proteingym_scores_directed$aa2)))
rownames(mean_tab_directed) <- unique_aas_directed
rownames(median_tab_directed) <- unique_aas_directed
colnames(mean_tab_directed) <- unique_aas_directed
colnames(median_tab_directed) <- unique_aas_directed
rownames(mean_tab_directed_min) <- unique_aas_directed
rownames(median_tab_directed_min) <- unique_aas_directed
colnames(mean_tab_directed_min) <- unique_aas_directed
colnames(median_tab_directed_min) <- unique_aas_directed
rownames(mean_tab_directed_max) <- unique_aas_directed
rownames(median_tab_directed_max) <- unique_aas_directed
colnames(mean_tab_directed_max) <- unique_aas_directed
colnames(median_tab_directed_max) <- unique_aas_directed

for (i in 1:nrow(proteingym_scores)) {
  aa1 <- proteingym_scores$aa1[i]
  aa2 <- proteingym_scores$aa2[i]
  mean_tab[aa1, aa2] <- proteingym_scores$mean_standard_score[i]
  mean_tab[aa2, aa1] <- proteingym_scores$mean_standard_score[i]
  median_tab[aa1, aa2] <- proteingym_scores$median_robust_score[i]
  median_tab[aa2, aa1] <- proteingym_scores$median_robust_score[i]
  
  # For the RSA category weighted values
  mean_weighted_mean <- (proteingym_scores_rsa[proteingym_scores_rsa$ref_aa == aa1 & proteingym_scores_rsa$mut_aa == aa2, 'weighted_mean_standard_score'] +
                           proteingym_scores_rsa[proteingym_scores_rsa$ref_aa == aa2 & proteingym_scores_rsa$mut_aa == aa1, 'weighted_mean_standard_score']) / 2
  
  mean_weighted_median <- (proteingym_scores_rsa[proteingym_scores_rsa$ref_aa == aa1 & proteingym_scores_rsa$mut_aa == aa2, 'weighted_median_robust_score'] +
                           proteingym_scores_rsa[proteingym_scores_rsa$ref_aa == aa2 & proteingym_scores_rsa$mut_aa == aa1, 'weighted_median_robust_score']) / 2

  mean_weighted_median_score <- (proteingym_scores_rsa[proteingym_scores_rsa$ref_aa == aa1 & proteingym_scores_rsa$mut_aa == aa2, 'weighted_median_standard_score'] +
                                 proteingym_scores_rsa[proteingym_scores_rsa$ref_aa == aa2 & proteingym_scores_rsa$mut_aa == aa1, 'weighted_median_standard_score']) / 2
  
  mean_tab_rsa[aa1, aa2] <- mean_weighted_mean
  mean_tab_rsa[aa2, aa1] <- mean_weighted_mean
  
  median_tab_rsa[aa1, aa2] <- mean_weighted_median
  median_tab_rsa[aa2, aa1] <- mean_weighted_median
  
  median_score_rsa[aa1, aa2] <- mean_weighted_median_score
  median_score_rsa[aa2, aa1] <- mean_weighted_median_score
  
  # For the directed tables
  mean_mean_val <- (proteingym_scores_directed[proteingym_scores_directed$aa1 == aa1 & proteingym_scores_directed$aa2 == aa2, 'mean_standard_score'] +
                    proteingym_scores_directed[proteingym_scores_directed$aa1 == aa2 & proteingym_scores_directed$aa2 == aa1, 'mean_standard_score']) / 2
  mean_tab_directed[aa1, aa2] <- mean_mean_val
  mean_tab_directed[aa2, aa1] <- mean_mean_val
  
  mean_median_val <- (proteingym_scores_directed[proteingym_scores_directed$aa1 == aa1 & proteingym_scores_directed$aa2 == aa2, 'median_robust_score'] +
                      proteingym_scores_directed[proteingym_scores_directed$aa1 == aa2 & proteingym_scores_directed$aa2 == aa1, 'median_robust_score']) / 2
  median_tab_directed[aa1, aa2] <- mean_median_val
  median_tab_directed[aa2, aa1] <- mean_median_val
  
  # For the min and max values
  min_mean_val <- min(proteingym_scores_directed[proteingym_scores_directed$aa1 == aa1 & proteingym_scores_directed$aa2 == aa2, 'mean_standard_score'],
                      proteingym_scores_directed[proteingym_scores_directed$aa1 == aa2 & proteingym_scores_directed$aa2 == aa1, 'mean_standard_score'])
  mean_tab_directed_min[aa1, aa2] <- min_mean_val
  mean_tab_directed_min[aa2, aa1] <- min_mean_val
  
  min_median_val <- min(proteingym_scores_directed[proteingym_scores_directed$aa1 == aa1 & proteingym_scores_directed$aa2 == aa2, 'median_robust_score'],
                        proteingym_scores_directed[proteingym_scores_directed$aa1 == aa2 & proteingym_scores_directed$aa2 == aa1, 'median_robust_score'])
  median_tab_directed_min[aa1, aa2] <- min_median_val
  median_tab_directed_min[aa2, aa1] <- min_median_val
  
  max_mean_val <- max(proteingym_scores_directed[proteingym_scores_directed$aa1 == aa1 & proteingym_scores_directed$aa2 == aa2, 'mean_standard_score'],
                      proteingym_scores_directed[proteingym_scores_directed$aa1 == aa2 & proteingym_scores_directed$aa2 == aa1, 'mean_standard_score'])
  mean_tab_directed_max[aa1, aa2] <- max_mean_val
  mean_tab_directed_max[aa2, aa1] <- max_mean_val
  
  
  
}


write.table(mean_tab_rsa,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_mean_standard_score_rsa_weighted.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)

write.table(median_tab_rsa,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_median_robust_score_rsa_weighted.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)

write.table(median_score_rsa,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_median_standard_score_rsa_weighted.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)
