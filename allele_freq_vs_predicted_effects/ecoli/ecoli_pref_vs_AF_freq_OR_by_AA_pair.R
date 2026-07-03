rm(list = ls(all.names = TRUE))

source("~/Drive/research/aa_distance/aa_distance_explore/allele_freq_vs_predicted_effects/functions_maf_vs_predicted_impact_fisher_exact.R")

# Compute odds ratios of (predicted) most extreme substitutions among rare vs freq substitutions.
# Do for VespaG and then for each AA pair separately for the focal new DISTASIS measure.
# Meant to show the (pretty trivial) point that VespaG and related tools convey a lot more information than
# simply knowing what the AA pairs are.

ecoli_tab <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/parsed_allelome_variants_w_prefs_dssp_w_anc.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)

ecoli_tab_forward <- ecoli_tab[which(ecoli_tab$consensus_aa == ecoli_tab$Likely_Anc_AA), ]
ecoli_tab_forward$anc_aa <- ecoli_tab_forward$consensus_aa
ecoli_tab_forward$anc_count <- ecoli_tab_forward$consensus_count
ecoli_tab_forward$derived_aa <- ecoli_tab_forward$alt_aa
ecoli_tab_forward$derived_count <- ecoli_tab_forward$alt_count

ecoli_tab_backward <- ecoli_tab[which(ecoli_tab$alt_aa == ecoli_tab$Likely_Anc_AA), ]
ecoli_tab_backward$anc_aa <- ecoli_tab_backward$alt_aa
ecoli_tab_backward$anc_count <- ecoli_tab_backward$alt_count
ecoli_tab_backward$derived_aa <- ecoli_tab_backward$consensus_aa
ecoli_tab_backward$derived_count <- ecoli_tab_backward$consensus_count

ecoli_tab_backward$AF <- 1 - ecoli_tab_backward$AF

ecoli_tab <- rbind(ecoli_tab_forward, ecoli_tab_backward)

ecoli_tab$count <- ecoli_tab$anc_count + ecoli_tab$derived_count


aa_combo_tallies <- table(sub_tab_filt$aa_combo)

atleast50_obs <- names(aa_combo_tallies)[which(aa_combo_tallies >= 50)]

sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$aa_combo %in% atleast50_obs), ]


# Get odds ratio and run Fisher's exact test for number of AA pair subs among rare and frequent subs.
fisher_summary <- fisher_tests_by_maf_and_metric_quantile(intab = sub_tab_filt,
                                                          prefs = c("thermoMPNN", "thermoMPNN_sum_scaled", "vespag", "vespag_sum_scaled"),
                                                          maf_cutoffs = c(0.001, 0.01))

write.table(fisher_summary,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)


# Also get enrichment for each AA pair separately.
aa_comb_fisher_summary <- fisher_tests_by_maf_and_aa_pair(intab = sub_tab_filt,
                                                          maf_cutoffs = c(0.001, 0.01))

write.table(aa_comb_fisher_summary,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_aa_pairs_by_freq_fisher_OR.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
