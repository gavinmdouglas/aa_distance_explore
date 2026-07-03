
fisher_tests_by_maf_and_aa_pair <- function(intab,
                                            maf_cutoffs = c(0.0001, 0.001, 0.01)) {
  
  raw <- list()
  aa_combos <- unique(intab$aa_combo)
  for (maf_cutoff in maf_cutoffs) {
    
    working <- intab
    working$freq_class <- NA
    working$freq_class[which(working$AF < maf_cutoff)] <- 'Rare'
    working$freq_class[which(working$AF >= maf_cutoff)] <- 'Higher'
    
    for (aa_comb in aa_combos) {
        working$combo_group <- NA
        working$combo_group[working$aa_combo == aa_comb] <- 'Focal'
        working$combo_group[working$aa_combo != aa_comb] <- 'Other'
        
        rare_focal_count <- length(which(working$freq_class == 'Rare' & working$combo_group == 'Focal'))
        rare_other_count <- length(which(working$freq_class == 'Rare' & working$combo_group == 'Other'))
        higher_focal_count <- length(which(working$freq_class == 'Higher' & working$combo_group == 'Focal'))
        higher_other_count <- length(which(working$freq_class == 'Higher' & working$combo_group == 'Other'))
        
        working_matrix <- matrix(c(rare_focal_count, higher_focal_count, rare_other_count, higher_other_count), nrow = 2, ncol = 2)
        
        # Add 1 to all cells if there are any 0's.
        if (any(working_matrix == 0)) {
          working_matrix <- working_matrix + 1
        }
        
        fisher_test <- fisher.test(working_matrix)
        
        raw[[paste(maf_cutoff, aa_comb, sep = '_')]] <- data.frame(maf_cutoff = maf_cutoff,
                                                                        aa_combo = aa_comb,
                                                                        rare_subset_count = rare_focal_count,
                                                                        higher_subset_count = higher_focal_count,
                                                                        rare_other_count = rare_other_count,
                                                                        higher_other_count = higher_other_count,
                                                                        fisher_OR = fisher_test$estimate,
                                                                        fisher_lower_95_OR = fisher_test$conf.int[1],
                                                                        fisher_upper_95_OR = fisher_test$conf.int[2],
                                                                        fisher_p = fisher_test$p.value)
        
        
    }
    
  }
  
  fisher_out <- do.call(rbind, raw)
  rownames(fisher_out) <- NULL
  
  return(fisher_out)
}



fisher_tests_by_maf_and_metric_quantile <- function(intab,
                                                    prefs,
                                                    maf_cutoffs = c(0.0001, 0.001, 0.01),
                                                    quantile_cutoffs = c(0.001, 0.01, 0.1, 0.25)) {
  
  raw <- list()
  for (maf_cutoff in maf_cutoffs) {
    
    working <- intab
    working$freq_class <- NA
    working$freq_class[which(working$AF < maf_cutoff)] <- 'Rare'
    working$freq_class[which(working$AF >= maf_cutoff)] <- 'Higher'
    
    for (cutoff in quantile_cutoffs) {
      for (pref in prefs) {
        pref_quantile <- quantile(working[, pref], probs = cutoff)
        
        working$quantile_group <- NA
        working$quantile_group[which(working[, pref] <= pref_quantile)] <- 'Focal'
        working$quantile_group[which(working[, pref] > pref_quantile)] <- 'Other'
        
        rare_focal_count <- length(which(working$freq_class == 'Rare' & working$quantile_group == 'Focal'))
        rare_other_count <- length(which(working$freq_class == 'Rare' & working$quantile_group == 'Other'))
        higher_focal_count <- length(which(working$freq_class == 'Higher' & working$quantile_group == 'Focal'))
        higher_other_count <- length(which(working$freq_class == 'Higher' & working$quantile_group == 'Other'))
        
        working_matrix <- matrix(c(rare_focal_count, higher_focal_count, rare_other_count, higher_other_count), nrow = 2, ncol = 2)
        
        # Add 1 to all cells if there are any 0's.
        if (any(working_matrix == 0)) {
          working_matrix <- working_matrix + 1
        }
        
        fisher_test <- fisher.test(working_matrix)
        
        raw[[paste(maf_cutoff, cutoff, pref, sep = '_')]] <- data.frame(maf_cutoff = maf_cutoff,
                                                                        pref_cutoff = cutoff,
                                                                        pref = pref,
                                                                        rare_subset_count = rare_focal_count,
                                                                        higher_subset_count = higher_focal_count,
                                                                        rare_other_count = rare_other_count,
                                                                        higher_other_count = higher_other_count,
                                                                        fisher_OR = fisher_test$estimate,
                                                                        fisher_lower_95_OR = fisher_test$conf.int[1],
                                                                        fisher_upper_95_OR = fisher_test$conf.int[2],
                                                                        fisher_p = fisher_test$p.value)
        
        
      }
    }
    
  }
  
  fisher_out <- do.call(rbind, raw)
  rownames(fisher_out) <- NULL
  
  return(fisher_out)
}

