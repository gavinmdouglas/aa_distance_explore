rm(list = ls(all.names = TRUE))

library(ggplot2)
library(reshape2)

# Similar to the analysis looking at quantiles of extreme preferences vs.
# proportion that were rare, look at most extreme exchangeabilities by site
# vs. proportion that are invariant compared to variant sites (excluding singletons).

metrics_map <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

ecoli_tab <- read.table('/Users/gavin/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_per_codon_exchangeability_invariant_vs_freq.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)
colnames(ecoli_tab)[which(colnames(ecoli_tab) == 'VespaG')] <- 'vespag'

prefs <- colnames(ecoli_tab)[7:ncol(ecoli_tab)]

# Only consider sites with at least 60,000 unambiguous sites.
# Should have been filtered already but can always check!
ecoli_tab <- ecoli_tab[which(ecoli_tab$N >= 60000), ]
ecoli_tab <- ecoli_tab[-which(ecoli_tab$Max_AC == 1), ]

ecoli_tab$Class <- NA
ecoli_tab$Class[which(ecoli_tab$Max_AC == 0)] <- 'Invariant'
ecoli_tab$Class[which(ecoli_tab$Max_AF < 0.0001 & ecoli_tab$Max_AF > 0)] <- 'Rare'
ecoli_tab$Class[which(ecoli_tab$Max_AF >= 0.0001)] <- 'Higher'

if (length(which(is.na(ecoli_tab$Class))) > 0) { stop('ERROR - not classifiable sites?') }

for (pref in prefs) {
  ecoli_tab[, pref] <- as.numeric(scale(ecoli_tab[, pref], center=TRUE, scale=TRUE))
}

ecoli_long <- melt(data = ecoli_tab,
                   id.vars = 'Class',
                   measure.vars = prefs,
                   variable.name = 'Pref',
                   value.name = 'exchange')
ecoli_long$Species <- 'E. coli'

# Read in human data.
prefs <- c(prefs, 'rasp')

human_tab <- read.table('/Users/gavin/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/human_variants/human_per_codon_exchangeability_invariant_vs_freq.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)

colnames(human_tab)[which(colnames(human_tab) == 'VespaG')] <- 'vespag'

human_tab <- human_tab[-which(human_tab$Max_AC == 1), ]

human_tab$Class <- NA
human_tab$Class[which(human_tab$Max_AC == 0)] <- 'Invariant'
human_tab$Class[which(human_tab$Max_AF < 0.0001 & human_tab$Max_AF > 0)] <- 'Rare'
human_tab$Class[which(human_tab$Max_AF >= 0.0001)] <- 'Higher'

if (length(which(is.na(human_tab$Class))) > 0) { stop('ERROR - not classifiable sites?') }

for (pref in prefs) {
  human_tab[, pref] <- as.numeric(scale(human_tab[, pref], center=TRUE, scale=TRUE))
}

human_long <- melt(data = human_tab,
                   id.vars = 'Class',
                   measure.vars = prefs,
                   variable.name = 'Pref',
                   value.name = 'exchange')
human_long$Species <- 'Human'

combined_data <- rbind(ecoli_long, human_long)
combined_data$Class <- factor(combined_data$Class, rev(c('Invariant', 'Rare', 'Higher')))

combined_data$Pref <- as.character(combined_data$Pref)
combined_data$Pref <- metrics_map[combined_data$Pref, 'Clean']

# Keep just the Grantham and Miyata original models, for simplicity.
prefs_to_rm <- c("Grantham (update, recalc.)", "Grantham (update, min-max)", "Grantham (recalc.)", "Grantham (min-max)",
                 "Miyata (recalc.)", "Miyata (min-max)", "Miyata (update, min-max)", "Miyata (update)")
combined_data <- combined_data[which(! combined_data$Pref %in% prefs_to_rm), ]
combined_data$Pref[which(combined_data$Pref == "Grantham (orig.)")] <- 'Grantham'
combined_data$Pref[which(combined_data$Pref == "Miyata (orig.)")] <- 'Miyata'

combined_data$Species[which(combined_data$Species == 'E. coli')] <- "italic('E. coli')"
combined_data$Species <- factor(combined_data$Species, levels = c("italic('E. coli')", "Human"))
combined_data$Pref <- factor(combined_data$Pref, levels = rev(sort(unique(combined_data$Pref))))

# Subset data to avoid memory limit issue.
combined_data <- combined_data[sample(1:nrow(combined_data), 1000000), ]

exchange_by_freq_boxplots <- ggplot(data = combined_data, aes(x = exchange, y = Pref, fill = Class)) +
                                  geom_boxplot(outlier.shape = NA) +
                                  xlab('Scaled exchangeability per site') +
                                  ylab('Preference type') +
                                  theme_bw() +
                                  guides(fill = guide_legend(reverse = TRUE)) +
                                  scale_fill_manual(values=c('#1f78b4', '#ff7f00', '#33a02c')) +
                                  labs(fill = 'Site class') +
                                  coord_cartesian(xlim = c(-3, 3)) +
                                  facet_wrap(Species ~ ., labeller = label_parsed)

ggsave(plot = exchange_by_freq_boxplots,
       filename = '~/Drive/ncsu/aa_selection/aa_manuscript/figures/Supp_exchangeability_by_freq.png',
       dpi = 300,
       width = 6,
       height = 7,
       units = 'in',
       device = 'png')
