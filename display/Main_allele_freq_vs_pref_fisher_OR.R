rm(list = ls(all.names = TRUE))

library(grid)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)

ecoli_ML_fisher_OR <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz',
                                 sep = '\t', stringsAsFactors = FALSE, header = TRUE)
ecoli_ML_fisher_OR <- ecoli_ML_fisher_OR[grep("standard", ecoli_ML_fisher_OR$pref, invert = TRUE), ]
ecoli_ML_fisher_OR <- ecoli_ML_fisher_OR[grep("sum_scaled", ecoli_ML_fisher_OR$pref, invert = TRUE), ]

ecoli_ML_fisher_OR <- ecoli_ML_fisher_OR[ecoli_ML_fisher_OR$pref_cutoff == 0.01, ]

human_ML_fisher_OR <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz',
                                 sep = '\t', stringsAsFactors = FALSE, header = TRUE)
human_ML_fisher_OR <- human_ML_fisher_OR[grep("standard", human_ML_fisher_OR$pref, invert = TRUE), ]
human_ML_fisher_OR <- human_ML_fisher_OR[human_ML_fisher_OR$pref_cutoff == 0.01, ]

ecoli_aa_pairs_fisher_OR <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_aa_pairs_by_freq_fisher_OR.tsv.gz',
                                       sep = '\t', stringsAsFactors = FALSE, header = TRUE)

# Filter out AA combos with low sample size:
ecoli_aa_pairs_fisher_OR <- ecoli_aa_pairs_fisher_OR[ecoli_aa_pairs_fisher_OR$higher_subset_count >= 10, ]

human_aa_pairs_fisher_OR <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_seg_subs_aa_pairs_by_freq_fisher_OR.tsv.gz',
                                       sep = '\t', stringsAsFactors = FALSE, header = TRUE)

ecoli_aa_pairs_fisher_OR$is_ML <- FALSE
ecoli_aa_pairs_fisher_OR$species <- "Escherichia coli"

ecoli_ML_points <- data.frame(maf_cutoff = ecoli_ML_fisher_OR$maf_cutoff,
                              aa_combo = ecoli_ML_fisher_OR$pref,
                              fisher_OR = ecoli_ML_fisher_OR$fisher_OR,
                              fisher_p = ecoli_ML_fisher_OR$fisher_p,
                              is_ML = TRUE,
                              species = "Escherichia coli",
                              pref = ecoli_ML_fisher_OR$pref)

ecoli_combined <- rbind(cbind(ecoli_aa_pairs_fisher_OR[, c("maf_cutoff", "aa_combo", "fisher_OR", "fisher_p", "is_ML", "species")], pref = NA),
                        ecoli_ML_points[, c("maf_cutoff", "aa_combo", "fisher_OR", "fisher_p", "is_ML", "species", "pref")])

human_aa_pairs_fisher_OR$is_ML <- FALSE
human_aa_pairs_fisher_OR$species <- "Human"

human_ML_points <- data.frame(maf_cutoff = human_ML_fisher_OR$maf_cutoff,
                              aa_combo = human_ML_fisher_OR$pref,
                              fisher_OR = human_ML_fisher_OR$fisher_OR,
                              fisher_p = human_ML_fisher_OR$fisher_p,
                              is_ML = TRUE,
                              species = "Human",
                              pref = human_ML_fisher_OR$pref)

human_combined <- rbind(cbind(human_aa_pairs_fisher_OR[, c("maf_cutoff", "aa_combo", "fisher_OR", "fisher_p", "is_ML", "species")], pref = NA),
                        human_ML_points[, c("maf_cutoff", "aa_combo", "fisher_OR", "fisher_p", "is_ML", "species", "pref")])

combined_data <- rbind(ecoli_combined, human_combined)

combined_data$is_significant <- p.adjust(combined_data$fisher_p, method = "BH") < 0.05

combined_data$colour_cat <- ifelse(combined_data$is_ML, combined_data$pref, "aa_pair")

combined_data$maf_cutoff_str <- factor(combined_data$maf_cutoff, 
                                       levels = sort(unique(combined_data$maf_cutoff)),
                                       labels = c("0.001", "0.01"))

combined_or_plot <- ggplot(combined_data, aes(x = maf_cutoff_str, y = log2(fisher_OR))) +
                      geom_quasirandom(aes(colour = colour_cat, alpha = is_significant, size = is_ML)) +
                      scale_colour_manual(values = c("aa_pair" = "black",
                                                     "vespag" = "red", 
                                                     "thermoMPNN" = "purple",
                                                     "rasp" = "blue"),
                                          labels = c("aa_pair" = "Amino acid pairs", "vespag" = "VespaG", "rasp" = "RaSP"),
                                          guide = guide_legend(override.aes = list(
                                            size = c(1, 2, 2),
                                            alpha = c(1, 1, 1)
                                          ))) +
                      scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) +
                      scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
                    
                      facet_wrap(~ species, labeller = labeller(species = c("Escherichia coli" = expression(italic("Escherichia coli")), "Human" = "Human"))) +
                      labs(x = "Minor allele frequency cut-off", y = expression(log[2]("Odds ratio"))) +
                      theme_bw() +
                      theme(legend.position = "right", 
                            strip.text = element_text(face = "italic"),
                            strip.background = element_blank(),
                            axis.line = element_line())

ggsave("~/Drive/research/aa_distance/aa_distance_ms/display/RAW_Main_allele_freq_fisher_OR_comparison.pdf", 
       width = 8, height = 4, units = "in", dpi = 300)

# Note that this plot needs to be cleaned up manually (in Affinity Designer) to format legend properly
# and un-italicize "Human".
