rm(list = ls(all.names = TRUE))

# Get rank of each distance matrix by a measure (LnL, AIC, or BIC), per lineage and subsample.
compute_model_ranks <- function(df, measure_col) {
  df$Model_rank <- NA
  raw <- list()
  for (lineage in unique(df$lineage)) {
    lineage_df <- df[df$lineage == lineage, ]
    for (subsample in unique(lineage_df$subsample)) {
      subsample_df <- lineage_df[lineage_df$subsample == subsample, ]
      
      if (measure_col == 'lnL') {
        subsample_df <- subsample_df[order(subsample_df[, measure_col], decreasing = TRUE), ]
      } else if (measure_col %in% c('AIC', 'BIC')) {
        subsample_df <- subsample_df[order(subsample_df[, measure_col], decreasing = FALSE), ] 
      } else {
        stop('measure_col must be either lnL, AIC, or BIC!') 
      }
      subsample_df$Model_rank <- 1:nrow(subsample_df)
      raw[[paste(lineage, subsample)]] <- subsample_df
    }
  }
  out_df <- do.call(rbind, raw)
  rownames(out_df) <- NULL
  return(out_df)
}

intab <- read.table(file="/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/combined_output.tsv.gz",
                    header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Ignore coarse DMS custom results, as I decided to just focus on the better version.
# Also remove EX2005-style metrics, which clearly did not work at all for these data, and mean "this study" to avoid confusion.
intab <- intab[! intab$distance_type %in% c("proteinGym_mean_standard_score_rsa_weighted", "proteinGym_median_robust_score_rsa_weighted",
                                            "proteinGym_rsa_median_ex2005_style", "proteinGym_rsa_mean_ex2005_style", "proteinGym_rsa_mean_custom"), ]

metric_id_map <- read.table('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                            header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names=1)

intab$distance_type <- sub('-', '.', intab$distance_type)
intab$distance_type <- sub('_prob', '', intab$distance_type)
intab$distance_type <- sub('_maxscaled', '', intab$distance_type)

intab$distance_type <- metric_id_map[intab$distance_type, 'Clean']

if (length(which(is.na(intab$distance_type))) > 0) {
  stop('Some distance types not found in metric_id_map!')
}

intab[which(intab$lineage == 'Strep'), 'lineage'] <- 'Streptococcus'

# Convert subsample to integer.
intab$subsample_num <- as.integer(gsub("subsample", "", intab$subsample))



# Add in number of sequences and codons per alignment. Number of codons will be used as n for calculating BIC.
alignment_info <- read.table(file="/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/subsampled_alignment_info.tsv.gz",
                             header=TRUE, sep="\t", stringsAsFactors = FALSE)
alignment_info[which(alignment_info$Group == 'Strep'), 'Group'] <- 'Streptococcus'

colnames(intab)[colnames(intab) == "ls"] <- "num_codons"
colnames(intab)[colnames(intab) == "ns"] <- "num_seqs"

miyata_grantham_i <- grep('Miyata|Grantham', intab$distance_type)
# But ignore those that contain "Combined".
miyata_grantham_i <- setdiff(miyata_grantham_i, grep("\\+", intab$distance_type))

miyata_grantham_nonorig_i <- setdiff(miyata_grantham_i, grep("orig.", intab$distance_type))

subset_miyata_grantham <- intab[miyata_grantham_i, ]
intab <- intab[-miyata_grantham_nonorig_i, ]

intab$distance_type <- gsub(' \\(orig.\\)', '', intab$distance_type)

# Split into subsets of 1-10 and 11-20
# Ignore 11-20 for Miyata/Grantham.
subset_1_10 <- intab[intab$subsample_num <= 10, ]
subset_11_20 <- intab[intab$subsample_num > 10, ]
subset_miyata_grantham_1_10 <- subset_miyata_grantham[subset_miyata_grantham$subsample_num <= 10, ]

# Keep only the "This study + EX" (DEX) combined metric and put the others in separate table for plotting.
combined_metrics_to_rm <- unique(grep("\\+", subset_11_20$distance_type, value=TRUE))
combined_only_11_20 <- subset_11_20[c(grep("\\+", subset_11_20$distance_type), which(subset_11_20$distance_type == "DEX")), ]

combined_only_11_20[combined_only_11_20$distance_type == "DEX", "distance_type"] <- "DEX (DMS-EX + EX)"

subset_11_20 <- subset_11_20[-which(subset_11_20$distance_type %in% combined_metrics_to_rm), ]

# Add in results for standard M0 model too, which is the standard model without AAdist.
M0_results <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/M0_combined_output.tsv.gz",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE)

M0_results[M0_results$lineage == 'Strep', 'lineage'] <- 'Streptococcus'

# Convert subsample to integer.
M0_results$subsample_num <- as.integer(gsub("subsample", "", M0_results$subsample))

M0_results$distance_type <- 'Standard M0 (no AAdist)'
M0_results$a <- NA
M0_results$b <- NA

M0_results$num_codons <- NA
M0_results$num_seqs <- NA
for (i in 1:nrow(alignment_info)) {
  lineage <- alignment_info$Group[i]
  subsample_num <- alignment_info$Subsample[i]
  M0_results$num_codons[M0_results$lineage == lineage & M0_results$subsample == subsample_num] <- alignment_info$NumCodons[i]
  M0_results$num_seqs[M0_results$lineage == lineage & M0_results$subsample == subsample_num] <- alignment_info$NumSeqs[i]
}

M0_results <- M0_results[, colnames(subset_11_20)]

M0_1_10 <- M0_results[M0_results$subsample_num <= 10, ]
M0_11_20 <- M0_results[M0_results$subsample_num > 10, ]

subset_1_10 <- rbind(subset_1_10, M0_1_10)
subset_11_20 <- rbind(subset_11_20, M0_11_20)

# Compute BIC:
subset_1_10$BIC <- log(subset_1_10$num_codons) * subset_1_10$num_param - 2 * subset_1_10$lnL
subset_11_20$BIC <- log(subset_11_20$num_codons) * subset_11_20$num_param - 2 * subset_11_20$lnL
combined_only_11_20$BIC <- log(combined_only_11_20$num_codons) * combined_only_11_20$num_param - 2 * combined_only_11_20$lnL
subset_miyata_grantham_1_10$BIC <- log(subset_miyata_grantham_1_10$num_codons) * subset_miyata_grantham_1_10$num_param - 2 * subset_miyata_grantham_1_10$lnL

subset_1_10 <- compute_model_ranks(df = subset_1_10, measure_col = 'BIC')
subset_11_20 <- compute_model_ranks(df = subset_11_20, measure_col = 'BIC')
combined_only_11_20 <- compute_model_ranks(df = combined_only_11_20, measure_col = 'BIC')
subset_miyata_grantham_1_10 <- compute_model_ranks(df = subset_miyata_grantham_1_10, measure_col = 'BIC')
subset_miyata_grantham_1_10$distance_type <- gsub(' \\(orig.\\)', '\\ (original\\)', subset_miyata_grantham_1_10$distance_type)
subset_miyata_grantham_1_10$distance_type <- gsub(' \\(recalc.\\)', '\\ (recalculated\\)', subset_miyata_grantham_1_10$distance_type)

write.table(x=subset_1_10,
            gzfile("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/focal_out_subsamples1_10.tsv.gz"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(x=subset_11_20,
            gzfile("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/focal_out_subsamples11_20.tsv.gz"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(x=combined_only_11_20,
            gzfile("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/other_out_combined_metrics_subsamples11_20.tsv.gz"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(x=subset_miyata_grantham_1_10,
            gzfile("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/other_out_miyata_grantham_subsamples1_10.tsv.gz"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
