library(exact2x2)

fisher_tests_by_maf_cutoff_aa_pair <- function(intab,
                                               count_cutoffs = 2,
                                               maf_cutoffs = c(0.0001, 0.001, 0.01),
                                               aa_combo_colname = "aa_combo",
                                               AF_colname = 'AF',
                                               count_colname = 'AC') {

  cutoffs <- c(count_cutoffs, maf_cutoffs)
  cutoff_types <- c(rep('count', length(count_cutoffs)), rep('maf', length(maf_cutoffs)))
  
  raw <- list()
  for (cutoff_i in 1:length(cutoffs)) {
    
    cutoff <- cutoffs[cutoff_i]
    cutoff_type <- cutoff_types[cutoff_i]
    
    working <- intab
    working$freq_class <- NA
    
    if (cutoff_type == 'maf') {
      working$freq_class <- ifelse(working[[AF_colname]] < cutoff, 'Rare', 'Higher')
    } else if (cutoff_type == 'count') {
      # Note for count it's <=, not <.
      working$freq_class <- ifelse(working[[count_colname]] <= cutoff, 'Rare', 'Higher')
    } else {
      stop('Invalid cutoff type:', cutoff_type) 
    }
    
    if (any(is.na(working$freq_class))) {
      stop('NA freq_class found after classification!')
    }
    
    unique_aa_combos <- sort(unique(working[[aa_combo_colname]]))
    
    for (aa_combo in unique_aa_combos) {
      num_focal_combo_rare <- sum(working[[aa_combo_colname]] == aa_combo & working$freq_class == 'Rare')
      num_focal_combo_higher <- sum(working[[aa_combo_colname]] == aa_combo & working$freq_class == 'Higher')
      num_other_combo_rare <- sum(working[[aa_combo_colname]] != aa_combo & working$freq_class == 'Rare')
      num_other_combo_higher <- sum(working[[aa_combo_colname]] != aa_combo & working$freq_class == 'Higher')
      
      working_matrix <- matrix(c(num_focal_combo_rare, num_focal_combo_higher,
                                 num_other_combo_rare, num_other_combo_higher), nrow = 2, ncol = 2)
      
      fisher_out <- exact2x2::exact2x2(working_matrix, alternative = "two.sided", conf.level = 0.95)
      
      raw[[paste(cutoff_type, cutoff, aa_combo, sep = '_')]] <- data.frame(cutoff_type = cutoff_type,
                                                                              cutoff = cutoff,
                                                                              aa_combo = aa_combo,
                                                                              focal_rare_count = num_focal_combo_rare,
                                                                              focal_higher_count = num_focal_combo_higher,
                                                                              other_rare_count = num_other_combo_rare,
                                                                              other_higher_count = num_other_combo_higher,
                                                                              fisher_OR = fisher_out$estimate,
                                                                              fisher_lower_95_OR = fisher_out$conf.int[1],
                                                                              fisher_upper_95_OR = fisher_out$conf.int[2],
                                                                              fisher_p = fisher_out$p.value)
    }
  
  }
  
  combined_fisher_out <- do.call(rbind, raw)
  rownames(combined_fisher_out) <- NULL
  
  return(combined_fisher_out)
  
}
