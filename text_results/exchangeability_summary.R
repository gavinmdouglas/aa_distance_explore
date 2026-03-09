rm(list = ls(all.names = TRUE))

library(reshape2)

by_aa_scaled <- read.table('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/codon_exchangeabilities_prepped_scaled.tsv.gz',
                           sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

sort(colMeans(by_aa_scaled))


# Perform pairwise correlation between all variables.

cor_matrix <- cor(t(by_aa_scaled), method = "spearman")

cor_long <- melt(cor_matrix, varnames = c("metric1", "metric2"), value.name = "rho")
cor_long <- cor_long[-which(cor_long$metric1 == cor_long$metric2), ]

for (i in 1:nrow(cor_long)) {
  m1 <- cor_long$metric1[i]
  m2 <- cor_long$metric2[i]
  sorted_m <- sort(c(m1, m2))
  cor_long$metric1[i] <- sorted_m[1]
  cor_long$metric2[i] <- sorted_m[2]
}

cor_long <- cor_long[! duplicated(cor_long), ]

median(cor_long$rho)
sd(cor_long$rho)
