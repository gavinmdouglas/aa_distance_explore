rm(list = ls(all.names = TRUE))

library(grid)
library(ggplot2)

# Rather than classifying rare and frequent, which doesn't make any sense realistically,
# what does make (more) sense is being able to classify substitutions with extremely low
# probabilities as rare (or invariant, as can be done in future analyses...).
# Look at proportion of substitutions that are rare based on the bottom X% of the scores for each preference.
# A higher proportion here is better!

metrics_map <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

sub_tab <- read.table('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/human_variants/human_snvs_w_prefs.tsv.gz',
                    header = TRUE, sep = '\t', stringsAsFactors = FALSE)

sub_tab$rasp <- sub_tab$rasp * -1

prefs <- colnames(sub_tab)[7:ncol(sub_tab)]

# Ignore singletons
sub_tab_filt <- sub_tab[which(sub_tab$variant_count > 1), ]

sub_tab_filt$AF <- sub_tab_filt$variant_count / sub_tab_filt$sample_count

# Only consider segregating sites with AF < 0.5 (as the reference allele was used as ancestral).
# Very few SNPs are thrown out by this anyway.
sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$AF < 0.5), ]

# Classify remaining subs as "rare" (AF < 0.0001) and "Higher" (>= 0.0001).
sub_tab_filt$freq_class <- NA
sub_tab_filt$freq_class[which(sub_tab_filt$AF < 0.0001)] <- 'Rare'
sub_tab_filt$freq_class[which(sub_tab_filt$AF >= 0.0001)] <- 'Higher'

raw <- list()
for (cutoff in c(0.001, 0.01, 0.1)) {
  for (pref in prefs) {
    pref_quantile <- quantile(sub_tab_filt[, pref], probs = cutoff)
    
    tab_subset <- sub_tab_filt[which(sub_tab_filt[, pref] <= pref_quantile), ]
    tab_other <- sub_tab_filt[which(sub_tab_filt[, pref] > pref_quantile), ]
    
    rare_subset_count <- length(which(tab_subset$freq_class == 'Rare'))
    higher_subset_count <- length(which(tab_subset$freq_class == 'Higher'))
    
    rare_other_count <- length(which(tab_other$freq_class == 'Rare'))
    higher_other_count <- length(which(tab_other$freq_class == 'Higher'))
    
    fisher_test <- fisher.test(x = matrix(c(rare_subset_count, higher_subset_count, rare_other_count, higher_other_count),
                                          nrow = 2, ncol = 2))
    
    raw[[paste(cutoff, pref, sep = '_')]] <- data.frame(cutoff = cutoff,
                                                        pref = pref,
                                                        rare_subset_count = rare_subset_count,
                                                        higher_subset_count = higher_subset_count,
                                                        rare_other_count = rare_other_count,
                                                        higher_other_count = higher_other_count,
                                                        fisher_OR = fisher_test$estimate,
                                                        fisher_lower_95_Orm(list = ls(all.names = TRUE))
                                                        
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
                                                         R = fisher_test$conf.int[1],
                                                        fisher_upper_95_OR = fisher_test$conf.int[2],
                                                        fisher_p = fisher_test$p.value)
  }
}

fisher_out <- do.call(rbind, raw)
rownames(fisher_out) <- NULL

gzfile_outfile <- gzfile('~/Drive/ncsu/aa_selection/aa_selection_zenodo/allele_freq_vs_predicted_effects/human_variants/human_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz', 'w')
write.table(fisher_out,
            file = gzfile_outfile,
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
close(gzfile_outfile)
