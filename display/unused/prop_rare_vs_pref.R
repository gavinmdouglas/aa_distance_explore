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

ecoli_combined_prop_rare <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/seg_subs_extreme_prefs_prop_rare.tsv.gz',
                                       sep = '\t', stringsAsFactors = FALSE, header = TRUE)

# combined_prop_rare$background_overall_rare_prop just included as column for convenience, should be one duplicated value.
if (length(unique(ecoli_combined_prop_rare$background_overall_rare_prop)) > 1) {
  stop('background_overall_rare_prop is not a single value!')
}

ecoli_overall_rare_prop <- ecoli_combined_prop_rare$background_overall_rare_prop[1]

# First look at E. coli observed substitutions vs. most extreme preferences.
ecoli_combined_prop_rare$cutoff <- factor(ecoli_combined_prop_rare$cutoff, levels = c('0.1', '0.01', '0.001'))

ecoli_combined_prop_rare <- ecoli_combined_prop_rare[which(! ecoli_combined_prop_rare$pref %in% prefs_to_rm), ]
ecoli_combined_prop_rare$pref[which(ecoli_combined_prop_rare$pref == "Grantham (orig.)")] <- 'Grantham'
ecoli_combined_prop_rare$pref[which(ecoli_combined_prop_rare$pref == "Miyata (orig.)")] <- 'Miyata'

# Sort by max and then by mean to get pref plotting order.
ecoli_prop_rare_prefs <- sort(unique(ecoli_combined_prop_rare$pref))
ecoli_combined_prop_rare_stats <- data.frame(pref=ecoli_prop_rare_prefs,
                                       max=sapply(ecoli_prop_rare_prefs, function(x) { max(ecoli_combined_prop_rare[which(ecoli_combined_prop_rare$pref == x), 'prop_rare']) }),
                                       mean=sapply(ecoli_prop_rare_prefs, function(x) { mean(ecoli_combined_prop_rare[which(ecoli_combined_prop_rare$pref == x), 'prop_rare']) }))
ecoli_combined_prop_rare_stats <- ecoli_combined_prop_rare_stats[order(ecoli_combined_prop_rare_stats$max, ecoli_combined_prop_rare_stats$mean, decreasing = FALSE), ]

ecoli_extremeprefs_rareprop <- ggplot(data = ecoli_combined_prop_rare, aes(x = prop_rare, y = pref, fill = cutoff)) +
                                    geom_bar(stat = 'identity', position = 'dodge') +
                                    theme_bw() +
                                    xlab('Proportion of sub. in cut-off quantile of lowest preferences that are rare\n(restricted to rare and frequent sub. only)') +
                                    ylab('Preference') +
                                    scale_y_discrete(limits = ecoli_combined_prop_rare_stats$pref) +
                                    scale_fill_manual(values = c('0.001' = 'grey30', '0.01' = 'grey60', '0.1' = 'grey90')) +
                                    coord_cartesian(xlim = c(0.85, 1)) +
                                    labs(fill = "Pref.\nquantile\ncut-off") +
                                    geom_vline(xintercept = ecoli_overall_rare_prop, linetype = 'dashed', color = 'blue', size = 0.5) +
                                    guides(fill = guide_legend(reverse = TRUE)) +
                                    ggtitle(expression(paste(italic("E. coli"), ' most negative sub. predictions vs. observed frequency'))) +
                                    theme(plot.title = element_text(hjust = 0.5),
                                          legend.position = 'none')


# Then similar idea, but invariant sites vs. frequent, based on mean codon exchangeabilities.
ecoli_combined_prop_invariant <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/sites_extreme_exchangeabilities_prop_invariant.tsv.gz',
                                            sep = '\t', stringsAsFactors = FALSE, header = TRUE)

# combined_prop_invariant$background_overall_rare_prop just included as column for convenience, should be one duplicated value.
if (length(unique(ecoli_combined_prop_invariant$background_overall_rare_prop)) > 1) {
  stop('background_overall_rare_prop is not a single value!')
}

ecoli_overall_invariant_prop <- ecoli_combined_prop_invariant$background_overall_rare_prop[1]

ecoli_combined_prop_invariant$cutoff <- factor(ecoli_combined_prop_invariant$cutoff, levels = c('0.1', '0.01', '0.001'))

# Only keep the original Grantham and Miyata as it is information overflow otherwise.
ecoli_combined_prop_invariant <- ecoli_combined_prop_invariant[which(! ecoli_combined_prop_invariant$pref %in% prefs_to_rm), ]
ecoli_combined_prop_invariant$pref[which(ecoli_combined_prop_invariant$pref == "Grantham (orig.)")] <- 'Grantham'
ecoli_combined_prop_invariant$pref[which(ecoli_combined_prop_invariant$pref == "Miyata (orig.)")] <- 'Miyata'

ecoli_prop_invariant_prefs <- sort(unique(ecoli_combined_prop_invariant$pref))
ecoli_combined_prop_invariant_stats <- data.frame(pref=ecoli_prop_invariant_prefs,
                                       max=sapply(ecoli_prop_invariant_prefs, function(x) { max(ecoli_combined_prop_invariant[which(ecoli_combined_prop_invariant$pref == x), 'prop_invariant']) }),
                                       mean=sapply(ecoli_prop_invariant_prefs, function(x) { mean(ecoli_combined_prop_invariant[which(ecoli_combined_prop_invariant$pref == x), 'prop_invariant']) }))
ecoli_combined_prop_invariant_stats <- ecoli_combined_prop_invariant_stats[order(ecoli_combined_prop_invariant_stats$max, ecoli_combined_prop_invariant_stats$mean, decreasing = FALSE), ]

ecoli_extremeprefs_invariantprop <- ggplot(data = ecoli_combined_prop_invariant, aes(x = prop_invariant, y = pref, fill = cutoff)) +
                                        geom_bar(stat = 'identity', position = 'dodge') +
                                        theme_bw() +
                                        xlab('Proportion of sites in cut-off quantile of lowest preferences that are invariant\n(restricted to sites that are invariant or have a frequent sub.)') +
                                        ylab('Preference') +
                                        scale_y_discrete(limits = ecoli_combined_prop_invariant_stats$pref) +
                                        scale_fill_manual(values = c('0.001' = 'grey30', '0.01' = 'grey60', '0.1' = 'grey90')) +
                                        coord_cartesian(xlim = c(0.85, 1)) +
                                        labs(fill = "Pref.\nquantile\ncut-off") +
                                        geom_vline(xintercept = ecoli_overall_invariant_prop, linetype = 'dashed', color = 'blue', size = 0.5) +
                                        guides(fill = guide_legend(reverse = TRUE)) +
                                        ggtitle(expression(paste(italic("E. coli"), ' least exchangeable codons vs. invariant status'))) +
                                        theme(plot.title = element_text(hjust = 0.5),
                                              legend.position = 'none')



human_combined_prop_rare <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/human_variants/seg_subs_extreme_prefs_prop_rare.tsv.gz',
                                       sep = '\t', stringsAsFactors = FALSE, header = TRUE)

# combined_prop_rare$background_overall_rare_prop just included as column for convenience, should be one duplicated value.
if (length(unique(human_combined_prop_rare$background_overall_rare_prop)) > 1) {
  stop('background_overall_rare_prop is not a single value!')
}

human_overall_rare_prop <- human_combined_prop_rare$background_overall_rare_prop[1]

# First look at E. coli observed substitutions vs. most extreme preferences.
human_combined_prop_rare$cutoff <- factor(human_combined_prop_rare$cutoff, levels = c('0.1', '0.01', '0.001'))

human_combined_prop_rare <- human_combined_prop_rare[which(! human_combined_prop_rare$pref %in% prefs_to_rm), ]
human_combined_prop_rare$pref[which(human_combined_prop_rare$pref == "Grantham (orig.)")] <- 'Grantham'
human_combined_prop_rare$pref[which(human_combined_prop_rare$pref == "Miyata (orig.)")] <- 'Miyata'

# Sort by max and then by mean to get pref plotting order.
human_prop_rare_prefs <- sort(unique(human_combined_prop_rare$pref))
human_combined_prop_rare_stats <- data.frame(pref=human_prop_rare_prefs,
                                             max=sapply(human_prop_rare_prefs, function(x) { max(human_combined_prop_rare[which(human_combined_prop_rare$pref == x), 'prop_rare']) }),
                                             mean=sapply(human_prop_rare_prefs, function(x) { mean(human_combined_prop_rare[which(human_combined_prop_rare$pref == x), 'prop_rare']) }))
human_combined_prop_rare_stats <- human_combined_prop_rare_stats[order(human_combined_prop_rare_stats$max, human_combined_prop_rare_stats$mean, decreasing = FALSE), ]

human_extremeprefs_rareprop <- ggplot(data = human_combined_prop_rare, aes(x = prop_rare, y = pref, fill = cutoff)) +
                                      geom_bar(stat = 'identity', position = 'dodge') +
                                      theme_bw() +
                                      xlab('Proportion of sub. in cut-off quantile of lowest preferences that are rare\n(restricted to rare and frequent sub. only)') +
                                      ylab('Preference') +
                                      scale_y_discrete(limits = human_combined_prop_rare_stats$pref) +
                                      scale_fill_manual(values = c('0.001' = 'grey30', '0.01' = 'grey60', '0.1' = 'grey90')) +
                                      coord_cartesian(xlim = c(0.99, 1)) +
                                      labs(fill = "Pref.\nquantile\ncut-off") +
                                      geom_vline(xintercept = human_overall_rare_prop, linetype = 'dashed', color = 'blue', size = 0.5) +
                                      guides(fill = guide_legend(reverse = TRUE)) +
                                      ggtitle('Human most negative sub. predictions vs. observed frequency') +
                                      theme(plot.title = element_text(hjust = 0.5),
                                            legend.position = 'none')


# Then similar idea, but invariant sites vs. frequent, based on mean codon exchangeabilities.
human_combined_prop_invariant <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/human_variants/sites_extreme_exchangeabilities_prop_invariant.tsv.gz',
                                            sep = '\t', stringsAsFactors = FALSE, header = TRUE)

# combined_prop_invariant$background_overall_rare_prop just included as column for convenience, should be one duplicated value.
if (length(unique(human_combined_prop_invariant$background_overall_rare_prop)) > 1) {
  stop('background_overall_rare_prop is not a single value!')
}

human_overall_invariant_prop <- human_combined_prop_invariant$background_overall_rare_prop[1]

human_combined_prop_invariant$cutoff <- factor(human_combined_prop_invariant$cutoff, levels = c('0.1', '0.01', '0.001'))

# Only keep the original Grantham and Miyata as it is information overflow otherwise.
human_combined_prop_invariant <- human_combined_prop_invariant[which(! human_combined_prop_invariant$pref %in% prefs_to_rm), ]
human_combined_prop_invariant$pref[which(human_combined_prop_invariant$pref == "Grantham (orig.)")] <- 'Grantham'
human_combined_prop_invariant$pref[which(human_combined_prop_invariant$pref == "Miyata (orig.)")] <- 'Miyata'

human_prop_invariant_prefs <- sort(unique(human_combined_prop_invariant$pref))
human_combined_prop_invariant_stats <- data.frame(pref=human_prop_invariant_prefs,
                                                  max=sapply(human_prop_invariant_prefs, function(x) { max(human_combined_prop_invariant[which(human_combined_prop_invariant$pref == x), 'prop_invariant']) }),
                                                  mean=sapply(human_prop_invariant_prefs, function(x) { mean(human_combined_prop_invariant[which(human_combined_prop_invariant$pref == x), 'prop_invariant']) }))
human_combined_prop_invariant_stats <- human_combined_prop_invariant_stats[order(human_combined_prop_invariant_stats$max, human_combined_prop_invariant_stats$mean, decreasing = FALSE), ]

human_extremeprefs_invariantprop <- ggplot(data = human_combined_prop_invariant, aes(x = prop_invariant, y = pref, fill = cutoff)) +
                                          geom_bar(stat = 'identity', position = 'dodge') +
                                          theme_bw() +
                                          xlab('Proportion of sites in cut-off quantile of lowest preferences that are invariant\n(restricted to sites that are invariant or have a frequent sub.)') +
                                          ylab('Preference') +
                                          scale_y_discrete(limits = human_combined_prop_invariant_stats$pref) +
                                          scale_fill_manual(values = c('0.001' = 'grey30', '0.01' = 'grey60', '0.1' = 'grey90')) +
                                          coord_cartesian(xlim = c(0.96, 1)) +
                                          labs(fill = "Pref.\nquantile\ncut-off") +
                                          geom_vline(xintercept = human_overall_invariant_prop, linetype = 'dashed', color = 'blue', size = 0.5) +
                                          guides(fill = guide_legend(reverse = TRUE)) +
                                          ggtitle('Human least exchangeable codons vs. invariant status') +
                                          theme(plot.title = element_text(hjust = 0.5),
                                                legend.position = 'none') +
                                          # Write value of Sneath at 0.001 (which is off the scale)
                                          annotate("text", x = 0.9583, y = 1.25, hjust = 0, vjust = 0, size = 3,
                                                   label = as.character(round(human_combined_prop_invariant[which(human_combined_prop_invariant$pref == 'Sneath' & human_combined_prop_invariant$cutoff == '0.001'), 'prop_invariant'], 2))) +
                                          annotate("text", x = 0.9583, y = 6.25, hjust = 0, vjust = 0, size = 3,
                                                   label = as.character(round(human_combined_prop_invariant[which(human_combined_prop_invariant$pref == 'RvC (charge)' & human_combined_prop_invariant$cutoff == '0.001'), 'prop_invariant'], 2)))
                                        

# Also parse out legend only for one plot.
grobs <- ggplotGrob(human_extremeprefs_invariantprop + theme(legend.position = 'right'))$grobs
legend_only <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

main_plot <- plot_grid(ecoli_extremeprefs_rareprop, ecoli_extremeprefs_invariantprop,
                       human_extremeprefs_rareprop, human_extremeprefs_invariantprop, 
                       labels = c('a', 'b', 'c', 'd'))

main_plot_w_legend <- plot_grid(main_plot, legend_only, ncol = 2, rel_widths = c(1, 0.1))

ggsave(plot = main_plot_w_legend,
       filename = '~/Drive/ncsu/aa_selection/aa_manuscript/figures/Main_allele_freq_by_extreme_prefs.png',
       dpi = 300,
       width = 14,
       height = 8,
       units = 'in',
       device = 'png')

ggsave(plot = main_plot_w_legend,
       filename = '~/Drive/ncsu/aa_selection/aa_manuscript/figures/Main_allele_freq_by_extreme_prefs.pdf',
       dpi = 600,
       width = 14,
       height = 8,
       units = 'in',
       device = 'pdf')
