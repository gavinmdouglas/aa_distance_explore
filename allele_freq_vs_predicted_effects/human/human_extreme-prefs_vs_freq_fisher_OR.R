rm(list = ls(all.names = TRUE))

source("~/Drive/research/aa_distance/aa_distance_explore/allele_freq_vs_predicted_effects/functions_maf_vs_predicted_impact_fisher_exact.R")

# Rather than classifying rare and frequent, which doesn't make any sense realistically,
# what does make (more) sense is being able to classify substitutions with extremely low
# probabilities as rare (or invariant, as can be done in future analyses...).
# Look at proportion of substitutions that are rare based on the bottom X% of the scores.
# A higher proportion here is better!

# Note that this is only done for VespaG and RaSP, which have per-site preferences.

# To get an idea of the enrichment explained by AA pairs alone, I also ran analogous enrichments based on AA pairs.

sub_tab <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_snvs_w_prefs.tsv.gz',
                    header = TRUE, sep = '\t', stringsAsFactors = FALSE)

sub_tab$rasp <- sub_tab$rasp * -1

# Ignore singletons.
sub_tab_filt <- sub_tab[which(sub_tab$variant_count > 1), ]

sub_tab_filt$AF <- sub_tab_filt$variant_count / sub_tab_filt$sample_count

# Only consider segregating sites with AF < 0.5 (as the reference allele was used as ancestral).
# Very few SNPs are thrown out by this anyway.
sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$AF < 0.5), ]

# Create sorted AA combinations
sub_tab_filt$aa_combo <- apply(sub_tab_filt[, c('ref_aa', 'alt_aa')], 1, function(x) {
  paste(sort(x), collapse = '_')
})

aa_combo_tallies <- table(sub_tab_filt$aa_combo)
atleast100_obs <- names(aa_combo_tallies)[which(aa_combo_tallies >= 100)]
sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$aa_combo %in% atleast100_obs), ]

# Get vespag and rasp standardized by protein too.
sub_tab_filt$vespag_standard <- NA
sub_tab_filt$rasp_standard <- NA
for (prot in unique(sub_tab_filt$protein_id)) {
  prot_inds <- which(sub_tab_filt$protein_id == prot)
  sub_tab_filt$vespag_standard[prot_inds] <- (sub_tab_filt$vespag[prot_inds] - mean(sub_tab_filt$vespag[prot_inds])) / sd(sub_tab_filt$vespag[prot_inds])
  sub_tab_filt$rasp_standard[prot_inds] <- (sub_tab_filt$rasp[prot_inds] - mean(sub_tab_filt$rasp[prot_inds])) / sd(sub_tab_filt$rasp[prot_inds])
}

sub_tab_filt <- sub_tab_filt[rowSums(is.na(sub_tab_filt)) == 0, ]

all_output <- fisher_tests_by_maf_and_metric_quantile(intab = sub_tab_filt,
                                                      prefs = c("vespag", "vespag_standard", "rasp", "rasp_standard"),
                                                      maf_cutoffs = c(0.0001, 0.001, 0.01))

write.table(all_output,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

# Also get enrichment for each AA pair separately.
aa_comb_fisher_summary <- fisher_tests_by_maf_and_aa_pair(intab = sub_tab_filt,
                                                          maf_cutoffs = c(0.0001, 0.001, 0.01))

write.table(aa_comb_fisher_summary,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_seg_subs_aa_pairs_by_freq_fisher_OR.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
