rm(list = ls(all.names = TRUE))

library(grid)
library(ggplot2)
library(exact2x2)

source('~/Drive/research/aa_distance/aa_distance_explore/allele_freq_vs_predicted_effects/function_compute_maf_vs_quantile_breakdown.R')

# Rather than classifying rare and frequent, which doesn't make any sense realistically,
# what does make (more) sense is being able to classify substitutions with extremely low
# probabilities as rare (or invariant, as can be done in future analyses...).
# Look at proportion of substitutions that are rare based on the bottom X% of the scores for each preference.
# A higher proportion here is better!
metrics_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

sub_tab <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_w_prefs.tsv.gz',
                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

prefs <- colnames(sub_tab)[8:ncol(sub_tab)]

# Only consider sites with at least 60,000 unambiguous sites.
sub_tab_filt <- sub_tab[which(sub_tab$N >= 60000), ]

# Ignore singletons
sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$AC > 1), ]

all_output <- fisher_tests_by_maf_and_metric_quantile(sub_tab_filt)

all_output$pref <- metrics_map[all_output$pref, 'Clean']

write.table(all_output,
            file =  gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
