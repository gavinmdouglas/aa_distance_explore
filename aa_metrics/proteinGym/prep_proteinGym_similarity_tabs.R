rm(list = ls(all.names = TRUE))

unique_aas <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

prep_symmetric_proteingym_sim_tab_w_rsa <- function(infile, focal_col, N_col1=NULL, N_col2=NULL, weighted_col=TRUE) {
  
  if (weighted_col) {
    if (is.null(N_col1) || is.null(N_col2)) {
      stop("If weighted_col is TRUE, both N_col1 and N_col2 must be provided.")
    }
  } else {
    if (!is.null(N_col1) || !is.null(N_col2)) {
      stop("If weighted_col is FALSE, N_col1 and N_col2 should not be provided.")
    }
  }
  
  intab <- read.table(infile, header=TRUE, sep='\t', stringsAsFactors = FALSE)
  
  prepped <- data.frame(matrix(0, nrow=20, ncol=20))
  rownames(prepped) <- unique_aas
  colnames(prepped) <- unique_aas
  
  for (aa1 in unique_aas[-length(unique_aas)]) {

    aa1_index <- which(unique_aas == aa1)

    for (aa2 in unique_aas[(aa1_index + 1):length(unique_aas)]) {

      forward_sub_row <- intab[intab$ref_aa == aa1 & intab$mut_aa == aa2, ]
      backward_sub_row <- intab[intab$ref_aa == aa2 & intab$mut_aa == aa1, ]
      
      if (weighted_col) {
        forward_N <- forward_sub_row[, N_col1] + forward_sub_row[, N_col2]
        backward_N <- backward_sub_row[, N_col1] + backward_sub_row[, N_col2]
        
        weighted_mean_val <- ((forward_sub_row[, focal_col] * forward_N) +
                              (backward_sub_row[, focal_col] * backward_N)) / (forward_N + backward_N)
  
        mean_val <- weighted_mean_val

      } else {
        mean_val <- mean(c(forward_sub_row[, focal_col], backward_sub_row[, focal_col])) 
      }
      
      prepped[aa1, aa2] <- mean_val
      prepped[aa2, aa1] <- mean_val

    }
  }
  
  return(prepped)
}

transform_to_range <- function(numbers, min_val = 0.01, max_val = 0.99) {
  input_min <- min(numbers)
  input_max <- max(numbers)
  return(((numbers - input_min) / (input_max - input_min)) * (max_val - min_val) + min_val)
}

# Prep the updated DMS-based values using our approach.
custom_filepath <- "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/parsed_dms_enrichments_to_aa_w_rsa.tsv"

weighted_median_custom_symmetric <- prep_symmetric_proteingym_sim_tab_w_rsa(custom_filepath, 'weighted_median_robust_score', N_col1="N_exposed_median_robust_score", N_col2="N_buried_median_robust_score")

write.table(weighted_median_custom_symmetric,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/proteinGym_rsa_median_custom.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)


# Also parse out asymmetric versions for future usage.
intab <- read.table(custom_filepath, header=TRUE, sep='\t', stringsAsFactors=FALSE)

subset_tab <- intab[, c("ref_aa", "mut_aa", "weighted_median_robust_score", "exposed_median_robust_score", "buried_median_robust_score")]

write.table(subset_tab,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/proteinGym_asymmetric_scores.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)


subset_tab$weighted <- transform_to_range(subset_tab$weighted_median_robust_score)
subset_tab$exposed <- transform_to_range(subset_tab$exposed_median_robust_score)
subset_tab$buried <- transform_to_range(subset_tab$buried_median_robust_score)

subset_tab <- subset_tab[, c("ref_aa", "mut_aa", "weighted", "exposed", "buried")]

write.table(subset_tab,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/proteinGym_asymmetric_similarity.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)
