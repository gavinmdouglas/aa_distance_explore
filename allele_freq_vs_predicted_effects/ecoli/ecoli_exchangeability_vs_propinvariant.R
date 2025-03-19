rm(list = ls(all.names = TRUE))

library(reshape2)

# Similar to the analysis looking at quantiles of extreme preferences vs.
# proportion that were rare, look at most extreme exchangeabilities by site
# vs. proportion that are invariant compared to variant sites (excluding singletons).

metrics_map <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

exchange_tab <- read.table('/Users/gavin/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_per_codon_exchangeability_invariant_vs_freq.tsv.gz',
                           header = TRUE, sep = '\t', stringsAsFactors = FALSE)

colnames(exchange_tab)[which(colnames(exchange_tab) == 'VespaG')] <- 'vespag'

prefs <- colnames(exchange_tab)[7:ncol(exchange_tab)]

# Only consider sites with at least 60,000 unambiguous sites.
# Should have been filtered already but can always check!
exchange_tab <- exchange_tab[which(exchange_tab$N >= 60000), ]
exchange_tab <- exchange_tab[-which(exchange_tab$Max_AC == 1), ]

exchange_tab$Class <- NA
exchange_tab$Class[which(exchange_tab$Max_AC == 0)] <- 'Invariant'
exchange_tab$Class[which(exchange_tab$Max_AF < 0.0001 & exchange_tab$Max_AF > 0)] <- 'Rare'
exchange_tab$Class[which(exchange_tab$Max_AF >= 0.0001)] <- 'Higher'

if (length(which(is.na(exchange_tab$Class))) > 0) { stop('ERROR - not classifiable sites?') }

for (pref in prefs) {
  exchange_tab[, pref] <- as.numeric(scale(exchange_tab[, pref], center=TRUE, scale=TRUE))
}

exchange_tab_long <- melt(data = exchange_tab,
                          id.vars = 'Class',
                          measure.vars = prefs,
                          variable.name = 'Pref',
                          value.name = 'exchange')
exchange_tab_long$Class <- factor(exchange_tab_long$Class, rev(c('Invariant', 'Rare', 'Higher')))

ggplot(data = exchange_tab_long, aes(x = exchange, y = Pref, fill = Class)) +
  geom_boxplot(outlier.shape = NA) +
  xlab('Scaled exchangeability') +
  ylab('Preference type') +
  theme_bw() +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values=c('#1f78b4', '#ff7f00', '#33a02c')) +
  labs(fill = 'Site\nclass') +
  coord_cartesian(xlim = c(-3, 3))

# Get proportion variant sites by AA.
variant_prop_by_aa <- data.frame(AA = sort(unique(exchange_tab$Consensus_AA)),
                                 prop_variant = NA)
rownames(variant_prop_by_aa) <- variant_prop_by_aa$AA
for (aa in unique(variant_prop_by_aa$AA)) {
  variant_prop_by_aa[aa, 'prop_variant'] <- length(which(exchange_tab$Consensus_AA == aa & exchange_tab$Class == 'Variant')) / length(which(exchange_tab$Consensus_AA == aa))
}

mean_invariant_vs_variant <- data.frame(pref = prefs,
                                       mean_invariant = NA,
                                       mean_variant = NA)
rownames(mean_invariant_vs_variant) <- mean_invariant_vs_variant$pref

rho <- numeric()
pval <- numeric()
for (pref in prefs) {
  mean_invariant_vs_variant[pref, 'mean_invariant'] <- mean(exchange_tab[which(exchange_tab$Class == 'Invariant'), pref])
  mean_invariant_vs_variant[pref, 'mean_variant'] <- mean(exchange_tab[which(exchange_tab$Class == 'Variant'), pref])
  exchange_tab$tmp <- exchange_tab[, pref]
  pref_subset <- aggregate(x = tmp ~ Consensus_AA, data = exchange_tab, FUN = mean)
  rownames(pref_subset) <- pref_subset$Consensus_AA
  variant_prop_by_aa[, pref] <- pref_subset[variant_prop_by_aa$AA, 'tmp']
  spearman_out <- cor.test(variant_prop_by_aa[, pref], variant_prop_by_aa$prop_variant, method = 'spearman', exact = FALSE)
  rho <- c(rho, spearman_out$estimate)
  pval <- c(pval, spearman_out$p.value)
}

mean_invariant_vs_variant$ratio <- mean_invariant_vs_variant$mean_variant / mean_invariant_vs_variant$mean_invariant

spearman_propvariant_vs_exchange <- data.frame(pref = prefs, rho = rho, p = pval)
spearman_propvariant_vs_exchange$BH <- p.adjust(spearman_propvariant_vs_exchange$p, method = 'BH')

spearman_propvariant_vs_exchange <- spearman_propvariant_vs_exchange[order(spearman_propvariant_vs_exchange$rho, decreasing = TRUE), ]

spearman_propvariant_vs_exchange$mean_variant_exchange <- mean_invariant_vs_variant[spearman_propvariant_vs_exchange$pref, 'mean_variant']
spearman_propvariant_vs_exchange$mean_invariant_exchange <- mean_invariant_vs_variant[spearman_propvariant_vs_exchange$pref, 'mean_invariant']

spearman_propvariant_vs_exchange$pref <- metrics_map[spearman_propvariant_vs_exchange$pref, 'Clean']

gz_outfile <- gzfile('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_exchange_vs_prop_variant.tsv.gz', 'w')
write.table(spearman_propvariant_vs_exchange,
            file = gz_outfile,
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
close(gz_outfile)
