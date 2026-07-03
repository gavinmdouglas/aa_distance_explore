rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ggbeeswarm)
library(cowplot)

# This package is needed to get italicized genera names as facet labels:
library(ggtext) 

clean_to_group <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                             sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 2)

clean_to_group["Standard M0", 'Grouping'] <- 'No amino acid distance'
group_colours <- c("Physicochemical/structural" = "black",
                   "Experimental"              = "#0072B2",
                   "Substitution-based"        = "#009E73",
                   "Radical vs. conservative"  = "#D55E00",
                   "No amino acid distance" = "grey50")

subset_1_10 <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/focal_out_subsamples1_10.tsv.gz",
                          header=TRUE, sep="\t", stringsAsFactors = FALSE)

subset_11_20 <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/focal_out_subsamples11_20.tsv.gz",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE)

subset_1_10[subset_1_10$distance_type == "Standard M0 (no AAdist)", 'distance_type'] <- "Standard M0"
subset_11_20[subset_11_20$distance_type == "Standard M0 (no AAdist)", 'distance_type'] <- "Standard M0"

combined_only_11_20 <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/other_out_combined_metrics_subsamples11_20.tsv.gz",
                                  header=TRUE, sep="\t", stringsAsFactors = FALSE)

subset_miyata_grantham_1_10 <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/other_out_miyata_grantham_subsamples1_10.tsv.gz",
                                          header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Set factors of distance metrics by mean rank, and set lineage factors.
factors_and_sort_by_mean_rank <- function(in_df) {
  dist_means <- aggregate(Model_rank ~ distance_type, data=in_df, FUN=mean)
  dist_means <- dist_means[order(dist_means$Model_rank, decreasing = TRUE), ]
  in_df$distance_type <- factor(in_df$distance_type, levels = dist_means$distance_type)
  
  # Italicize genera names.
  in_df$lineage <- gsub("Drosophila", "*Drosophila*", in_df$lineage)
  in_df$lineage <- gsub("Streptococcus", "*Streptococcus*", in_df$lineage)
  in_df$lineage <- factor(in_df$lineage, levels=c("*Streptococcus*", "*Drosophila*", "Mammal"))
  
  return(in_df)
}

subset_1_10 <- factors_and_sort_by_mean_rank(subset_1_10)
subset_miyata_grantham_1_10 <- factors_and_sort_by_mean_rank(subset_miyata_grantham_1_10)
combined_only_11_20 <- factors_and_sort_by_mean_rank(combined_only_11_20)
subset_11_20 <- factors_and_sort_by_mean_rank(subset_11_20)

colour_y_labels <- function(in_df) {
  levels_order <- levels(in_df$distance_type)
  groups <- clean_to_group[levels_order, 'Grouping']
  in_df$group <- clean_to_group[as.character(in_df$distance_type), 'Grouping']
  in_df$group <- factor(in_df$group, levels = names(group_colours))
  return(in_df)
}

subset_1_10 <- colour_y_labels(subset_1_10)
subset_11_20 <- colour_y_labels(subset_11_20)

y_colours_1_10 <- group_colours[clean_to_group[levels(subset_1_10$distance_type), 'Grouping']]
y_colours_11_20 <- group_colours[clean_to_group[levels(subset_11_20$distance_type), 'Grouping']]

# Plot supplementary plot with subsamples 1-10 for focal categories.
subset_1_10_rank_plot <- ggplot(data=subset_1_10, aes(x = Model_rank, y = distance_type)) +
  geom_quasirandom(aes(colour = group), size=0.8, orientation='y') +
  geom_boxplot(aes(colour = group), alpha=0.5, outlier.shape = NA) +
  facet_wrap(~lineage, nrow = 1) +
  theme_bw() +
  scale_colour_manual(values = group_colours, name = "Measure grouping") +
  theme(strip.background = element_blank(),
        strip.text = element_markdown(),
        axis.text.y = element_text(colour = y_colours_1_10)) +
  xlab('Model rank (lower is better)') +
  ylab('Amino acid distance basis for model')

# Along with another supplementary plot with Grantham and Miyata variants for first 10 subsamples.
grantham_miyata_rank_plot <- ggplot(data=subset_miyata_grantham_1_10, aes(x = Model_rank, y = distance_type)) +
  geom_quasirandom(col='black', size=0.8, orientation='y') +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  facet_wrap(~lineage, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_markdown()) +
  xlab('Model rank (lower is better)') +
  ylab('Amino acid distance basis for model') +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0.0, NA), breaks= seq(0, 10, by=2))

# Also get plot of "Combined" separately.
combined_only_11_20_rank_plot <- ggplot(data=combined_only_11_20, aes(x = Model_rank, y = distance_type)) +
  geom_quasirandom(col='black', size=0.8, orientation='y') +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  facet_wrap(~lineage, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_markdown()) +
  xlab('Model rank (lower is better)') +
  ylab('Amino acid distance basis for model') +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0.0, NA), breaks= seq(0, 12, by=2))

# Then get plot of subsamples 11-20 for focal categories.
subset_11_20_rank_plot <- ggplot(data=subset_11_20, aes(x = Model_rank, y = distance_type)) +
  geom_quasirandom(aes(colour = group), size=0.8, orientation='y') +
  geom_boxplot(aes(colour = group), alpha=0.5, outlier.shape = NA) +
  facet_wrap(~lineage, nrow = 1) +
  theme_bw() +
  scale_colour_manual(values = group_colours, name = "Measure grouping") +
  theme(strip.background = element_blank(),
        strip.text = element_markdown(),
        axis.text.y = element_text(colour = y_colours_11_20))+
  xlab('Model rank (lower is better)') +
  ylab('Amino acid distance basis for model')

# Also get a and b values for subsets 11-20 in same order.
a_vs_b_scatterplot <- ggplot(data = subset_11_20[subset_11_20$distance_type == "DEX", ], aes(x=a, y=b, colour=lineage)) +
                              geom_point() +
                              theme_bw() +
                              theme(strip.background = element_blank(),
                                    strip.text = element_markdown(),
                                    legend.text = element_markdown(),
                                    legend.position = "inside",
                                    legend.position.inside = c(0.7, 0.3),
                                    legend.background = element_rect(colour = "black")) +
                              scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0.0, NA)) +
                              scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0.0, NA)) +
                              scale_colour_manual(name = "Lineage",
                                                  values = c("#E69F00", "#56B4E9", "#009E73")) +
                              guides(colour = guide_legend(override.aes = list(size = 3))) +
                              xlab(expression(paste("Overall nonsynonymous acceptance rate (", italic(a), ")"))) +
                              ylab(expression(paste("Amino acid distance penalty strength (", italic(b), ")")))

# Write out plots.
ggsave(plot = subset_1_10_rank_plot,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Main_model_rank_subsample1-10.pdf',
       width = 8, height = 6, useDingbats = FALSE)

ggsave(plot = grantham_miyata_rank_plot,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Supp_Grantham_Miyata_model_rank_subsample1-10.pdf',
       width = 8, height = 6, useDingbats = FALSE)

ggsave(plot = combined_only_11_20_rank_plot,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Supp_combined_metrics_model_rank_subsample11-20.pdf',
       width = 8, height = 6, useDingbats = FALSE)

ggsave(plot = subset_11_20_rank_plot,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Main_metrics_model_rank_subsample11-20.pdf',
       width = 8, height = 6, useDingbats = FALSE)

ggsave(plot = a_vs_b_scatterplot,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Main_a_vs_b_custom_metric.pdf',
       width = 3.5, height = 3.5, useDingbats = FALSE)
