rm(list = ls(all.names = TRUE))

pairwise_diff_normalized <- function(tab, normtype="mean_diff", set_mean_100=FALSE, sig_digits=2) {
  
  # General function for computing Miyata and Grantham-like amino acid pairwise distances. 
  # Allows different properties to be input, and also different ways of normalizing pairwise differences.
  
  # "tab" needs to be a dataframe with single-letter amino acid codes as rownames.
  # The columns must be the properties of interest that will be compared.
  
  # "norm_type is the type of normalization to be used. Options are:
  # 1) "mean_diff" - The mean difference between all pairs of amino acids (based on each column individually) is used to normalize the differences per column.
  # 2) "sd_diff" - The standard deviation of the difference between all pairs of amino acids (based on each column individually) is used to normalize the differences per column.
  # 3) "standard" - Prior to taking differences, convert all variables to standard scores (mean centre and divide by SD). No normalization of absolute differences is done.
  # 4) "minmax" - Prior to taking differences convert all variables to min-max scaled values (subtract min and divide by range). No normalization of absolute differences is done.
  
  # "set_mean_100" to True if you want the mean distance to be 100 (all values will be scaled to make this the case).

  # "sig_digits" is the number of significant digits to round the output to.
  
  amino_acids <- c("F", "L", "S", "Y", "C", "W",
                   "P", "H", "Q", "R", "I", "M",
                   "T", "N", "K", "V", "A", "D",
                   "E", "G")
  aa_pairs <- combn(amino_acids, 2, simplify = FALSE)

  # Input parameter checks.
  if (! all(rownames(tab) %in% amino_acids)) { stop("Not all amino acids are present in the rownames.") }
  if (nrow(tab) != 20) { stop("The number of rows is not equal to 20.") }
  if (any(is.na(tab))) { stop("There are NA values in the input dataframe.") }
  if (ncol(tab) < 1) { stop("There needs to be at least one column.") }
  if (! (normtype %in% c("mean_diff", "sd_diff", "minmax", "standard"))) { stop("normtype must be one of 'mean_diff', 'sd_diff', 'minmax', or 'standard'.") }
  if (!all(sapply(tab, is.numeric))) { stop("All columns must be numeric.") }

  if (normtype %in% c("mean_diff", "sd_diff")) {
    mean_diffs <- numeric()
    sd_diffs <- numeric()
  
    for (category in colnames(tab)) {
      category_pairwise_abs_diff <- sapply(aa_pairs, function(x) { abs(tab[x[1], category] -  tab[x[2], category]) })
      mean_diffs <- c(mean_diffs, mean(category_pairwise_abs_diff))
      sd_diffs <- c(sd_diffs, sd(category_pairwise_abs_diff))
    }
  } else if (normtype == "minmax") {
    tab <- as.data.frame(apply(tab, 2, function(x) { (x - min(x)) / (max(x) - min(x)) }))
  } else if (normtype == "standard") {
    tab <- as.data.frame(apply(tab, 2, function(x) { (x - mean(x)) / sd(x) }))
  }
  
  out_dist <- data.frame(matrix(NA, nrow = 20, ncol = 20))
  rownames(out_dist) <- amino_acids
  colnames(out_dist) <- amino_acids
  
  # Set the diagonal to 0.
  diag(out_dist) <- 0
  
  for (i in 1:length(aa_pairs)) {
    aa1 <- aa_pairs[[i]][1]
    aa2 <- aa_pairs[[i]][2]
    
    delta <- as.numeric(tab[aa1, ] - tab[aa2, ])

    if (normtype == "mean_diff") {
      combined_dist <- sqrt(sum((delta / mean_diffs) ** 2))
    } else if (normtype == "sd_diff") {
      combined_dist <- sqrt(sum((delta / sd_diffs) ** 2))
    } else if (normtype == "minmax" || normtype == "standard") {
      combined_dist <- sqrt(sum(delta ** 2))
    } else {
     stop('Error: normtype should have matched one of those categories...') 
    }
  
   out_dist[aa1, aa2] <- combined_dist
   out_dist[aa2, aa1] <- combined_dist
    
  }
  
  if (set_mean_100) {
    out_dist <- out_dist * (100 / mean(out_dist[upper.tri(out_dist)]))
  }
  
  return(signif(out_dist, digits = sig_digits))

}
