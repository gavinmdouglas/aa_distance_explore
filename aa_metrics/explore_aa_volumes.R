rm(list = ls(all.names = TRUE))

# Compare mean volumes, and compare standard deviations of these volumes by approach.
# Quick analysis to confirm that the Tsai volumes have the lowest SD on average
# (although not significantly with KW test).

harpaz <- read.table('~/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/literature_data/Harpaz1994_Table1_volumes.tsv',
                            header = TRUE, sep = '\t', stringsAsFactors = FALSE)

tsai <- read.table('~/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/literature_data/Tsai2002_Table3_volumes.tsv',
                   header = TRUE, sep = '\t', stringsAsFactors = FALSE)

harpaz <- harpaz[-which(is.na(harpaz$AA)), ]
rownames(harpaz) <- harpaz$AA

tsai <- tsai[-which(is.na(tsai$AA)), ]
rownames(tsai) <- tsai$AA
tsai <- tsai[rownames(harpaz), ]

# Confirmed that they are all very similar.
# Note that unnamed "volume" column is from Harpaz 1994.
combined_means <- cbind(harpaz, tsai)
combined_means <- combined_means[, c('mean_volume', 'ProtOr_volume', 'Chothia_volume', 'Richards_volume')]

pairs(combined_means)
cor(combined_means)
colMeans(cor(combined_means))

combined_sds <- cbind(harpaz, tsai)
combined_sds <- combined_sds[, c('sd_volume', 'ProtOR_SD', 'Chothia_SD', 'Richards_SD')]

# Note that these are not significantly different:
kruskal.test(combined_sds)

# Nonetheless, seems better to choose the one with lowest SD rather than to take the average volume across approaches.
colMeans(combined_sds)

# Lowest is ProtOR-based approach.
