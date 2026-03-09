rm(list = ls(all.names = TRUE))

ecoli_spearman <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_spearman_pref_vs_maf.tsv.gz',
                             header = TRUE, sep = '\t', stringsAsFactors = FALSE)




human_spearman <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_spearman_pref_vs_maf.tsv.gz',
                             header = TRUE, sep = '\t', stringsAsFactors = FALSE)

ecoli_spearman$Species <- 'E. coli'
human_spearman$Species <- 'Human'

combined <- rbind(ecoli_spearman, human_spearman)

extra_to_rm <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/nonfocal_metrics.txt',
                           header = FALSE, stringsAsFactors = FALSE)$V1

extra_to_rm <- c("proteinGym_rsa_mean_custom", extra_to_rm)

combined <- combined[! combined$Metric %in% extra_to_rm, ]
combined$Metric <- sub('-', '.', combined$Metric)
combined$Metric <- sub('_prob', '', combined$Metric)
combined$Metric <- sub('.prob', '', combined$Metric)
combined$Metric <- sub('_maxscaled', '', combined$Metric)

metrics_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

combined$Metric <- metrics_map[combined$Metric, 'Clean']
combined$Metric[which(combined$Metric == "Grantham (orig.)")] <- 'Grantham'
combined$Metric[which(combined$Metric == "Miyata (orig.)")] <- 'Miyata'
combined$Metric[which(combined$Metric == "This study (median)")] <- 'This study'

combined$Spearman_BH <- p.adjust(combined$Spearman_p, method = 'BH')

combined[combined$Metric == 'Combined (This study + EX)', ]

human_only <- combined[combined$Species == 'Human', ]
tail(human_only[order(human_only$Spearman_rho), ])


ecoli_only <- combined[combined$Species == 'E. coli', ]
tail(ecoli_only[order(ecoli_only$Spearman_rho), ])
