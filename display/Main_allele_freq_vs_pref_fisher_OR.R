rm(list = ls(all.names = TRUE))

library(grid)
library(ggplot2)
library(cowplot)

# Figure summarizing how well observed non-syn allele frequencies across E. coli and human
# agree with site AA preferences from a variety of approaches.

aa_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/aa_map.tsv.gz',
                     sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 3)

metrics_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)


ecoli_fisher_OR <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz',
                                       sep = '\t', stringsAsFactors = FALSE, header = TRUE)
ecoli_fisher_OR$Species <- 'E. coli'


human_fisher_OR <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz',
                              sep = '\t', stringsAsFactors = FALSE, header = TRUE)
human_fisher_OR$Species <- 'Human'

fisher_OR <- rbind(ecoli_fisher_OR, human_fisher_OR)

fisher_OR <- fisher_OR[which(fisher_OR$maf_cutoff == 0.01 & fisher_OR$pref_cutoff == 0.01), ]

fisher_OR$pref <- metrics_map[fisher_OR$pref, 'Clean']

# Only keep the original Grantham and Miyata as it is information overflow otherwise.
prefs_to_rm <- c("Grantham (update, recalc.)", "Grantham (update, min-max)", "Grantham (recalc.)", "Grantham (min-max)",
                 "Miyata (recalc.)", "Miyata (min-max)", "Miyata (update, min-max)", "Miyata (update)")
fisher_OR <- fisher_OR[which(! fisher_OR$pref %in% prefs_to_rm), ]
fisher_OR$pref[which(fisher_OR$pref == "Grantham (orig.)")] <- 'Grantham'
fisher_OR$pref[which(fisher_OR$pref == "Miyata (orig.)")] <- 'Miyata'

fisher_OR$Significance <- ifelse(fisher_OR$fisher_p < 0.05, 'P < 0.05', 'Non-sig.')
fisher_OR$Significance <- factor(fisher_OR$Significance, levels = c('Non-sig.', 'P < 0.05'))

fisher_OR$Species[which(fisher_OR$Species == 'E. coli')] <- "italic('Escherichia coli')"
fisher_OR$Species <- factor(fisher_OR$Species, levels = c("italic('Escherichia coli')", "Human"))

pref_means <- aggregate(data = fisher_OR, fisher_OR ~ pref, FUN = mean)
pref_means <- pref_means[order(pref_means$fisher_OR, decreasing = FALSE), ]

pref_means_ML <- pref_means[which(pref_means$pref %in% c('RaSP', 'VespaG')), ]
pref_means_metrics <- pref_means[which(! pref_means$pref %in% c('RaSP', 'VespaG')), ]
pref_means <- rbind(pref_means_metrics, pref_means_ML)

fisher_OR$pref <- factor(fisher_OR$pref, levels = pref_means$pref)

OR_summary_plot <- ggplot(data = fisher_OR,
                          aes(x = log2(fisher_OR),
                              y = pref,
                              colour = Significance)) +
  geom_errorbarh(aes(xmin = log2(fisher_lower_95_OR), xmax = log2(fisher_upper_95_OR)), height = 0.2, colour='black') +
  geom_point(size = 3, alpha=0.8) +
  theme_bw() +
  xlab(expression(log[2]~"Odds ratio (enrichment of predicted divergent substitutions among rare mutations)")) +
  ylab("Preference type") +
  geom_vline(xintercept = 0.0, lty = 2) +
  facet_wrap(Species ~ ., labeller = label_parsed, scales = 'free_x') +
  scale_colour_manual(values=c('grey80', 'cornflowerblue')) +
  theme(axis.text.y = element_text(color = ifelse(levels(factor(fisher_OR$pref)) %in% c('RaSP', 'VespaG'), "#4CAF50", "black")))

ggsave(plot = OR_summary_plot,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Main_allele_freq_by_extreme_prefs_Fisher.pdf',
       dpi = 600,
       width = 7,
       height = 5,
       units = 'in',
       device = 'pdf')
