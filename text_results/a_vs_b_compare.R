rm(list = ls(all.names = TRUE))

subset_11_20 <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/focal_out_subsamples11_20.tsv.gz",
                           header=TRUE, sep="\t", stringsAsFactors = FALSE)

combined_only <- subset_11_20[subset_11_20$distance_type == "Combined (This study + EX)", ]

cor.test(combined_only$a, combined_only$b, method="spearman")

# Decided to just mention mean values per genus.
kruskal.test(combined_only$a ~ combined_only$lineage)
kruskal.test(combined_only$b ~ combined_only$lineage)

mean(combined_only$a[combined_only$lineage == "Streptococcus"])
mean(combined_only$b[combined_only$lineage == "Streptococcus"])

mean(combined_only$a[combined_only$lineage == "Drosophila"])
mean(combined_only$b[combined_only$lineage == "Drosophila"])

mean(combined_only$a[combined_only$lineage == "Mammal"])
mean(combined_only$b[combined_only$lineage == "Mammal"])
