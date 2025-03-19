rm(list = ls(all.names = TRUE))

library(grid)
library(ggplot2)
library(cowplot)

# Figure summarizing how well observed non-syn allele frequencies across E. coli and human
# agree with site AA preferences from a variety of approaches.

aa_map <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/aa_map.tsv.gz',
                     sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 3)

metrics_map <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# Only keep the original Grantham and Miyata as it is information overflow otherwise.
prefs_to_rm <- c("Grantham (update, recalc.)", "Grantham (update, min-max)", "Grantham (recalc.)", "Grantham (min-max)",
                 "Miyata (recalc.)", "Miyata (min-max)", "Miyata (update, min-max)", "Miyata (update)")

ecoli_fisher_OR <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz',
                                       sep = '\t', stringsAsFactors = FALSE, header = TRUE)
ecoli_fisher_OR$Species <- 'E. coli'


human_fisher_OR <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/human_variants/human_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz',
                              sep = '\t', stringsAsFactors = FALSE, header = TRUE)
human_fisher_OR$Species <- 'Human'

fisher_OR <- rbind(ecoli_fisher_OR, human_fisher_OR)
fisher_OR$pref <- metrics_map[fisher_OR$pref, 'Clean']
fisher_OR$`Quantile (%)` <- NA
fisher_OR$`Quantile (%)`[which(fisher_OR$cutoff == '0.001')] <- '0.1'
fisher_OR$`Quantile (%)`[which(fisher_OR$cutoff == '0.01')] <- '1'
fisher_OR$`Quantile (%)`[which(fisher_OR$cutoff == '0.1')] <- '10'
fisher_OR$`Quantile (%)` <- factor(fisher_OR$`Quantile (%)`, levels = c('10', '1', '0.1'))

fisher_OR <- fisher_OR[which(! fisher_OR$pref %in% prefs_to_rm), ]
fisher_OR$pref[which(fisher_OR$pref == "Grantham (orig.)")] <- 'Grantham'
fisher_OR$pref[which(fisher_OR$pref == "Miyata (orig.)")] <- 'Miyata'

fisher_OR$Sig. <- ifelse(fisher_OR$fisher_p < 0.05, 'P < 0.05', 'Non-sig.')
fisher_OR$Sig. <- factor(fisher_OR$Sig., levels = c('Non-sig.', 'P < 0.05'))

fisher_OR$Species[which(fisher_OR$Species == 'E. coli')] <- "italic('E. coli')"
fisher_OR$Species <- factor(fisher_OR$Species, levels = c("italic('E. coli')", "Human"))

pref_means <- aggregate(data = fisher_OR, fisher_OR ~ pref, FUN = mean)
pref_means <- pref_means[order(pref_means$fisher_OR, decreasing = FALSE), ]

fisher_OR$pref <- factor(fisher_OR$pref, levels = pref_means$pref)

OR_summary_plot <- ggplot(data = fisher_OR,display_items/Main_allele_freq_vs_pref_fisher_OR.R
                          aes(x = fisher_OR,
                              y = pref,
                              colour = Sig.,
                              shape = `Quantile (%)`)) +
  geom_point(size = 3) +
 # xlim(c(0.5, 1.75)) +
  theme_bw() +
  xlab("Odd's ratio (enrichment in rare mutations)") +
  ylab("Preference type") +
  geom_vline(xintercept = 1, lty = 2) +
  facet_wrap(Species ~ ., labeller = label_parsed) +
  scale_colour_manual(values=c('grey80', 'cornflowerblue'))

ggsave(plot = OR_summary_plot,
       filename = '~/Drive/ncsu/aa_selection/aa_manuscript/figures/Main_allele_freq_by_extreme_prefs_Fisher.png',
       dpi = 300,
       width = 7,
       height = 5,
       units = 'in',
       device = 'png')

ggsave(plot = OR_summary_plot,
       filename = '~/Drive/ncsu/aa_selection/aa_manuscript/figures/Main_allele_freq_by_extreme_prefs_Fisher.pdf',
       dpi = 600,
       width = 7,
       height = 5,
       units = 'in',
       device = 'pdf')
