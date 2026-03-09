rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ggbeeswarm)
library(cowplot)

# This package is needed to get italicized genera names as facet labels:
library(ggtext) 

subset_1_10 <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/focal_out_subsamples1_10.tsv.gz",
                          header=TRUE, sep="\t", stringsAsFactors = FALSE)

mean(subset_1_10$Model_rank[subset_1_10$distance_type == "DMS-EX"])
mean(subset_1_10$Model_rank[subset_1_10$distance_type == "EX"])

# Get mean difference in AIC (and SD).
aic_diff <- numeric()
for (lineage in unique(subset_1_10$lineage)) {
  lineage_subset <- subset_1_10[subset_1_10$lineage == lineage, ]
   for (subsample_num in unique(lineage_subset$subsample_num)) {
     subsample_subset <- lineage_subset[lineage_subset$subsample_num == subsample_num, ]
     dms_ex_aic <- subsample_subset$AIC[subsample_subset$distance_type == "DMS-EX"]
     ex_aic <- subsample_subset$AIC[subsample_subset$distance_type == "EX"]
     aic_diff <- c(aic_diff, dms_ex_aic - ex_aic) 
   }
}

hist(aic_diff)
wilcox.test(aic_diff)

mean(subset_1_10$AIC[subset_1_10$distance_type == "DMS-EX" & subset_1_10$lineage == "Drosophila"])
mean(subset_1_10$AIC[subset_1_10$distance_type == "EX" & subset_1_10$lineage == "Drosophila"])

mean(subset_1_10$AIC[subset_1_10$distance_type == "DMS-EX" & subset_1_10$lineage == "Mammal"])
mean(subset_1_10$AIC[subset_1_10$distance_type == "EX" & subset_1_10$lineage == "Mammal"])

mean(subset_1_10$AIC[subset_1_10$distance_type == "DMS-EX" & subset_1_10$lineage == "Streptococcus"])
mean(subset_1_10$AIC[subset_1_10$distance_type == "EX" & subset_1_10$lineage == "Streptococcus"])


mean(subset_1_10$AIC[subset_1_10$distance_type == "DMS-EX"])
mean(subset_1_10$AIC[subset_1_10$distance_type == "EX"])


combined_only_11_20 <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/other_out_combined_metrics_subsamples11_20.tsv.gz",
                                  header=TRUE, sep="\t", stringsAsFactors = FALSE)
mean(combined_only_11_20$Model_rank[combined_only_11_20$distance_type == "DEX (DMS-EX + EX)"])

mean(combined_only_11_20$Model_rank[combined_only_11_20$distance_type == "EX + DeMaSk"])
mean(combined_only_11_20$Model_rank[combined_only_11_20$lineage == "Streptococcus" & combined_only_11_20$distance_type == "EX + DeMaSk"])
mean(combined_only_11_20$Model_rank[combined_only_11_20$lineage != "Streptococcus" & combined_only_11_20$distance_type == "EX + DeMaSk"])


subset_11_20 <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/focal_out_subsamples11_20.tsv.gz",
                           header=TRUE, sep="\t", stringsAsFactors = FALSE)
mean(subset_11_20$Model_rank[subset_11_20$distance_type == "DEX"])
table(subset_11_20$Model_rank[subset_11_20$distance_type == "DEX"])


# Get mean difference in AIC (and SD) between DEX and DMS-EX.
aic_diff <- numeric()
for (lineage in unique(subset_11_20$lineage)) {
  lineage_subset <- subset_11_20[subset_11_20$lineage == lineage, ]
  for (subsample_num in unique(lineage_subset$subsample_num)) {
    subsample_subset <- lineage_subset[lineage_subset$subsample_num == subsample_num, ]
    dms_ex_aic <- subsample_subset$AIC[subsample_subset$distance_type == "DMS-EX"]
    dex_aic <- subsample_subset$AIC[subsample_subset$distance_type == "DEX"]
    aic_diff <- c(aic_diff, dms_ex_aic - dex_aic) 
  }
}
median(aic_diff)
mean(aic_diff)
sd(aic_diff)

hist(aic_diff)
wilcox.test(aic_diff)
