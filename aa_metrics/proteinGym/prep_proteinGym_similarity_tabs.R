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

# First prepare the EX2005-style weighted mean and median similarity tables.
EX2005_parsed_filepath <- "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/parsed_EX2005_applied_to_proteinGym.tsv"
weighted_mean_ex2005_style <- prep_proteingym_sim_tab_w_rsa(EX2005_parsed_filepath, 'weighted_mean', N_col1="N_buried", N_col2="N_exposed")
weighted_median_ex2005_style <- prep_proteingym_sim_tab_w_rsa(EX2005_parsed_filepath, 'weighted_median', N_col1="N_buried", N_col2="N_exposed")

write.table(weighted_mean_ex2005_style,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_rsa_mean_ex2005_style.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)

write.table(weighted_median_ex2005_style,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_rsa_median_ex2005_style.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)


# Then prep the updated DMS-based values using an approach that to me is more robust and straight-forward.
custom_filepath <- "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/parsed_dms_enrichments_to_aa_w_rsa.tsv"

weighted_mean_custom <- prep_proteingym_sim_tab_w_rsa(custom_filepath, 'weighted_mean_standard_score', N_col1="N_exposed_mean_standard_score", N_col2="N_buried_mean_standard_score")

weighted_median_custom <- prep_proteingym_sim_tab_w_rsa(custom_filepath, 'weighted_median_robust_score', N_col1="N_exposed_median_robust_score", N_col2="N_buried_median_robust_score")

write.table(weighted_mean_custom,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_rsa_mean_custom.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)

write.table(weighted_median_custom,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_rsa_median_custom.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)
