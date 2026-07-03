rm(list = ls(all.names = TRUE))

# Compute correlations mean preferences vs. minor allele frequencies for each amino acid pair.
metrics_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

ecoli_tab <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_w_vespag.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# Get vespag standardized by protein too.
ecoli_tab$vespag_standard <- NA
for (prot in unique(ecoli_tab$protein)) {
  prot_inds <- which(ecoli_tab$protein == prot)
  ecoli_tab$vespag_standard[prot_inds] <- (ecoli_tab$vespag[prot_inds] - mean(ecoli_tab$vespag[prot_inds])) / sd(ecoli_tab$vespag[prot_inds])
}

aa_sim_metrics <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/combined_symmetric_prepped_similarity_metrics.tsv.gz",
                             header=TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

# Only consider sites with at least 60,000 unambiguous sites.
sub_tab_filt <- ecoli_tab[which(ecoli_tab$N >= 60000), ]

# Ignore mutations in start codon position.
sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$pos != 0), ]

# Ignore singletons
sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$AC > 1), ]

# Create sorted AA combinations
sub_tab_filt$aa_combo <- apply(sub_tab_filt[, c('consensus_aa', 'seg_aa')], 1, function(x) {
  paste(sort(x), collapse = '_')
})

aa_combo_tallies <- table(sub_tab_filt$aa_combo)

atleast100_obs <- names(aa_combo_tallies)[which(aa_combo_tallies >= 100)]

sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$aa_combo %in% atleast100_obs), ]

mean_by_aa_combo <- aggregate(cbind(AF, vespag, vespag_standard) ~ aa_combo,
                              data = sub_tab_filt,
                              FUN = mean)

rownames(mean_by_aa_combo) <- mean_by_aa_combo$aa_combo
mean_by_aa_combo <- mean_by_aa_combo[, -1, drop = FALSE]
for (aa_metric in colnames(aa_sim_metrics)) {
  mean_by_aa_combo[, aa_metric] <- aa_sim_metrics[rownames(mean_by_aa_combo), aa_metric]
}

tested_names <- character()
spearman_rho <- numeric()
spearman_p <- numeric()
for (coln in colnames(mean_by_aa_combo)) {
 if (coln == "AF") { next }
   test_res <- cor.test(mean_by_aa_combo[, coln], mean_by_aa_combo$AF, method = 'spearman', exact=FALSE)
   tested_names <- c(tested_names, coln)
   spearman_rho <- c(spearman_rho, test_res$estimate)
   spearman_p <- c(spearman_p, test_res$p.value)
}

spearman_tab <- data.frame(Metric = tested_names,
                           Spearman_rho = spearman_rho,
                           Spearman_p = spearman_p)

spearman_tab <- spearman_tab[order(spearman_tab$Spearman_rho), ]

write.table(x=spearman_tab,
            file=gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_spearman_pref_vs_maf.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
