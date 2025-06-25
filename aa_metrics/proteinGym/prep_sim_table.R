rm(list = ls(all.names = TRUE))

proteingym_scores <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_processed.tsv',
                                header=TRUE, sep = '\t', stringsAsFactors = FALSE)

rownames(proteingym_scores) <- paste0(proteingym_scores$aa1, '_', proteingym_scores$aa2)

# First run some sanity checks to make sure it is correlated with known similarity measures.
ex <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity/ex.maxscaled.tsv.gz',
                 header=TRUE, sep = '\t', stringsAsFactors = FALSE)

miyata <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity/miyata_orig.maxscaled.invert.tsv.gz',
                     header=TRUE, sep = '\t', stringsAsFactors = FALSE)

proteingym_scores$ex <- NA
proteingym_scores$miyata <- NA

for (i in 1:nrow(proteingym_scores)) {
 aa1 <- proteingym_scores$aa1[i]
 aa2 <- proteingym_scores$aa2[i]
 proteingym_scores[i, 'ex'] <- ex[ex$aa1 == aa1 & ex$aa2 == aa2, 'similarity']
 proteingym_scores[i, 'miyata'] <- miyata[miyata$aa1 == aa1 & miyata$aa2 == aa2, 'similarity']
}

proteingym_scores$mean_standard_score_prepped <- proteingym_scores$mean_standard_score + abs(min(proteingym_scores$mean_standard_score)) + 0.01
proteingym_scores$mean_standard_score_prepped <- proteingym_scores$mean_standard_score_prepped / (max(proteingym_scores$mean_standard_score_prepped) + 0.01)

proteingym_scores$median_robust_score_prepped <- proteingym_scores$median_robust_score + abs(min(proteingym_scores$median_robust_score)) + 0.01
proteingym_scores$median_robust_score_prepped <- proteingym_scores$median_robust_score_prepped / (max(proteingym_scores$median_robust_score_prepped) + 0.01)

proteingym_scores$median_orderedNorm <- bestNormalize::orderNorm(proteingym_scores$median_robust_score)$x.t
proteingym_scores$median_orderedNorm_prepped <- proteingym_scores$median_orderedNorm + abs(min(proteingym_scores$median_orderedNorm)) + 0.01
proteingym_scores$median_orderedNorm_prepped <- proteingym_scores$median_orderedNorm_prepped / (max(proteingym_scores$median_orderedNorm_prepped) + 0.01)
  
mean_tab <- data.frame(aa1=c(proteingym_scores$aa1, proteingym_scores$aa2),
                           aa2=c(proteingym_scores$aa2, proteingym_scores$aa1),
                           similarity=c(proteingym_scores$mean_standard_score_prepped, proteingym_scores$mean_standard_score_prepped),
                           stringsAsFactors = FALSE)

median_tab <- data.frame(aa1=c(proteingym_scores$aa1, proteingym_scores$aa2),
                       aa2=c(proteingym_scores$aa2, proteingym_scores$aa1),
                       similarity=c(proteingym_scores$median_robust_score_prepped, proteingym_scores$median_robust_score_prepped),
                       stringsAsFactors = FALSE)

orderedNorm_tab <- data.frame(aa1=c(proteingym_scores$aa1, proteingym_scores$aa2),
                         aa2=c(proteingym_scores$aa2, proteingym_scores$aa1),
                         similarity=c(proteingym_scores$median_orderedNorm_prepped, proteingym_scores$median_orderedNorm_prepped),
                         stringsAsFactors = FALSE)

write.table(mean_tab,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity/proteingym_mean_standard.tsv.gz'),
            sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(median_tab,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity/proteingym_median_robust.tsv.gz'),
            sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(orderedNorm_tab,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity/proteingym_median_robust_orderedNorm.tsv.gz'),
            sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

