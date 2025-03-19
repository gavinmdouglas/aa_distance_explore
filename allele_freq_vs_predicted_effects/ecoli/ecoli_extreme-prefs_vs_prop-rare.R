rm(list = ls(all.names = TRUE))

library(grid)
library(ggplot2)
library(exact2x2)

# Rather than classifying rare and frequent, which doesn't make any sense realistically,
# what does make (more) sense is being able to classify substitutions with extremely low
# probabilities as rare (or invariant, as can be done in future analyses...).
# Look at proportion of substitutions that are rare based on the bottom X% of the scores for each preference.
# A higher proportion here is better!
metrics_map <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

sub_tab <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_w_prefs.tsv.gz',
                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

prefs <- colnames(sub_tab)[8:ncol(sub_tab)]

# Only consider sites with at least 60,000 unambiguous sites.
sub_tab_filt <- sub_tab[which(sub_tab$N >= 60000), ]

# Ignore singletons
sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$AC > 1), ]

# Classify remaining subs as "Rare" (AF < 0.0001) and "Higher" (>= 0.0001)
sub_tab_filt$freq_class <- NA
sub_tab_filt$freq_class[which(sub_tab_filt$AF < 0.0001)] <- 'Rare'
sub_tab_filt$freq_class[which(sub_tab_filt$AF >= 0.0001)] <- 'Higher'

rare_overall <- length(which(sub_tab_filt$freq_class == 'Rare'))
higher_overall <- length(which(sub_tab_filt$freq_class == 'Higher'))

raw <- list()
for (cutoff in c(0.001, 0.01, 0.1)) {
  for (pref in prefs) {
    pref_quantile <- quantile(sub_tab_filt[, pref], probs = cutoff)

    tab_subset <- sub_tab_filt[which(sub_tab_filt[, pref] <= pref_quantile), ]
    tab_other <- sub_tab_filt[which(sub_tab_filt[, pref] > pref_quantile), ]

    rare_subset_count <- length(which(tab_subset$freq_class == 'Rare'))
    higher_subset_count <- length(which(tab_subset$freq_class == 'Higher'))
    
    rare_other_count <- length(which(tab_other$freq_class == 'Rare'))
    higher_other_count <- length(which(tab_other$freq_class == 'Higher'))

    fisher_test <- fisher.test(x = matrix(c(rare_subset_count, higher_subset_count, rare_other_count, higher_other_count),
                                         nrow = 2, ncol = 2))
    
    raw[[paste(cutoff, pref, sep = '_')]] <- data.frame(cutoff = cutoff,
                                                        pref = pref,
                                                        rare_subset_count = rare_subset_count,
                                                        higher_subset_count = higher_subset_count,
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

gzfile_outfile <- gzfile('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz', 'w')
write.table(fisher_out,
            file = gzfile_outfile,
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
close(gzfile_outfile)
