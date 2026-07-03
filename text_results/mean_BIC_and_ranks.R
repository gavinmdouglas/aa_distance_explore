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

# Get mean difference in BIC (and SD).
bic_diff <- numeric()
for (lineage in unique(subset_1_10$lineage)) {
  lineage_subset <- subset_1_10[subset_1_10$lineage == lineage, ]
   for (subsample_num in unique(lineage_subset$subsample_num)) {
     subsample_subset <- lineage_subset[lineage_subset$subsample_num == subsample_num, ]
     dms_ex_bic <- subsample_subset$BIC[subsample_subset$distance_type == "DMS-EX"]
     ex_bic <- subsample_subset$BIC[subsample_subset$distance_type == "EX"]
     bic_diff <- c(bic_diff, dms_ex_bic - ex_bic) 
   }
}

hist(bic_diff)
wilcox.test(bic_diff)

mean(subset_1_10$BIC[subset_1_10$distance_type == "DMS-EX" & subset_1_10$lineage == "Drosophila"])
mean(subset_1_10$BIC[subset_1_10$distance_type == "EX" & subset_1_10$lineage == "Drosophila"])

mean(subset_1_10$BIC[subset_1_10$distance_type == "DMS-EX" & subset_1_10$lineage == "Mammal"])
mean(subset_1_10$BIC[subset_1_10$distance_type == "EX" & subset_1_10$lineage == "Mammal"])

mean(subset_1_10$BIC[subset_1_10$distance_type == "DMS-EX" & subset_1_10$lineage == "Streptococcus"])
mean(subset_1_10$BIC[subset_1_10$distance_type == "EX" & subset_1_10$lineage == "Streptococcus"])


mean(subset_1_10$BIC[subset_1_10$distance_type == "DMS-EX"])
mean(subset_1_10$BIC[subset_1_10$distance_type == "EX"])


combined_only_11_20 <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/other_out_combined_metrics_subsamples11_20.tsv.gz",
                                  header=TRUE, sep="\t", stringsAsFactors = FALSE)
mean(combined_only_11_20$Model_rank[combined_only_11_20$distance_type == "DEX (DMS-EX + EX)"])



subset_11_20 <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/focal_out_subsamples11_20.tsv.gz",
                           header=TRUE, sep="\t", stringsAsFactors = FALSE)
mean(subset_11_20$Model_rank[subset_11_20$distance_type == "DEX"])
table(subset_11_20$Model_rank[subset_11_20$distance_type == "DEX"])


# Get mean difference in BIC (and SD) between DEX and DMS-EX.
bic_diff <- numeric()
for (lineage in unique(subset_11_20$lineage)) {
  lineage_subset <- subset_11_20[subset_11_20$lineage == lineage, ]
  for (subsample_num in unique(lineage_subset$subsample_num)) {
    subsample_subset <- lineage_subset[lineage_subset$subsample_num == subsample_num, ]
    dms_ex_bic <- subsample_subset$BIC[subsample_subset$distance_type == "DMS-EX"]
    dex_bic <- subsample_subset$BIC[subsample_subset$distance_type == "DEX"]
    bic_diff <- c(bic_diff, dms_ex_bic - dex_bic) 
  }
}
median(bic_diff)
mean(bic_diff)
sd(bic_diff)

hist(bic_diff)
wilcox.test(bic_diff)
