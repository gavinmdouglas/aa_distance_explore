rm(list = ls(all.names = TRUE))

source("~/Drive/research/aa_distance/aa_distance_explore/allele_freq_vs_predicted_effects/functions_maf_vs_predicted_impact_fisher_exact.R")

# Compute odds ratios of (predicted) most extreme substitutions among rare vs freq substitutions.
# Do for VespaG and then for each AA pair separately for the focal new DISTASIS measure.
# Meant to show the (pretty trivial) point that VespaG and related tools convey a lot more information than
# simply knowing what the AA pairs are.

ecoli_tab <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_w_vespag.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)
colnames(ecoli_tab)[which(colnames(ecoli_tab) == 'VespaG')] <- 'vespag'

# Only consider sites with at least 60,000 unambiguous sites.
sub_tab_filt <- ecoli_tab[which(ecoli_tab$N >= 60000), ]

# Ignore singletons
sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$AC > 1), ]

# Create sorted AA combinations
sub_tab_filt$aa_combo <- apply(sub_tab_filt[, c('consensus_aa', 'seg_aa')], 1, function(x) {
  paste(sort(x), collapse = '_')
})

aa_combo_tallies <- table(sub_tab_filt$aa_combo)

atleast100_obs <- names(aa_combo_tallies)[which(aa_combo_tallies >= 100)]

sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$aa_combo %in% atleast100_obs), ]

# Get vespag standardized by protein too.
sub_tab_filt$vespag_standard <- NA
for (prot in unique(sub_tab_filt$protein)) {
  prot_inds <- which(sub_tab_filt$protein == prot)
  sub_tab_filt$vespag_standard[prot_inds] <- (sub_tab_filt$vespag[prot_inds] - mean(sub_tab_filt$vespag[prot_inds])) / sd(sub_tab_filt$vespag[prot_inds])
}

sub_tab_filt <- sub_tab_filt[rowSums(is.na(sub_tab_filt)) == 0, ]

# Get odds ratio and run Fisher's exact test for number of AA pair subs among rare and frequent subs.
fisher_summary <- fisher_tests_by_maf_and_metric_quantile(intab = sub_tab_filt,
                                                          prefs = c("vespag", "vespag_standard"),
                                                          maf_cutoffs = c(0.0001, 0.001, 0.01))

write.table(fisher_summary,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)


# Also get enrichment for each AA pair separately.
aa_comb_fisher_summary <- fisher_tests_by_maf_and_aa_pair(intab = sub_tab_filt,
                                                          maf_cutoffs = c(0.0001, 0.001, 0.01))

write.table(aa_comb_fisher_summary,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_aa_pairs_by_freq_fisher_OR.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
