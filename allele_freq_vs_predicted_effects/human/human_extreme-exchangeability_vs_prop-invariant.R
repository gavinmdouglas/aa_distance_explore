rm(list = ls(all.names = TRUE))

library(grid)
library(ggplot2)

# Similar to the analysis looking at quantiles of extreme preferences vs.
# proportion that were rare, look at most extreme exchangeabilities by site
# vs. proportion that are invariant compared to frequent (at least one sub >= 0.01 AF).

aa_map <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/aa_map.tsv.gz',
                     sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 3)

metrics_map <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

exchange_tab <- read.table('/Users/gavin/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/human_variants/human_per_codon_exchangeability_invariant_vs_freq.tsv.gz',
                           header = TRUE, sep = '\t', stringsAsFactors = FALSE)

exchange_tab$Class <- 'Variant'
exchange_tab[which(exchange_tab$Max_AC == 0), 'Class'] <- 'Invariant'

colnames(exchange_tab)[which(colnames(exchange_tab) == 'VespaG')] <- 'vespag'

prefs <- colnames(exchange_tab)[6:ncol(exchange_tab)]

prefs <- prefs[which(! prefs %in% 'Class')]

raw <- list()

for (cutoff in c(0.001, 0.01, 0.1)) {
  
  prop_invariant <- c()
  
  for (pref in prefs) {
    pref_quantile <- quantile(exchange_tab[, pref], probs = cutoff)
    tab_subset <- exchange_tab[which(exchange_tab[, pref] <= pref_quantile), ]
    prop_invariant <- c(prop_invariant, length(which(tab_subset$Class == 'Invariant')) / nrow(tab_subset))
  }
  
  prep_prop_invariant_at_extreme <- data.frame(cutoff=as.character(cutoff), pref = prefs, prop_invariant = prop_invariant)
  prep_prop_invariant_at_extreme <- prep_prop_invariant_at_extreme[order(prep_prop_invariant_at_extreme$prop_invariant, decreasing = TRUE), ]
  
  prep_prop_invariant_at_extreme$pref <- metrics_map[prep_prop_invariant_at_extreme$pref, 'Clean']
  
  raw[[as.character(cutoff)]] <- prep_prop_invariant_at_extreme
  
}

combined_prop_invariant <- do.call(rbind, raw)

combined_prop_invariant$pref <- as.character(combined_prop_invariant$pref)

combined_prop_invariant$background_overall_invariant_prop <- length(which(exchange_tab$Class == 'Invariant')) / nrow(exchange_tab)

write.table(combined_prop_invariant,
            file = gzfile('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/human_variants/sites_extreme_exchangeabilities_prop_invariant.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
