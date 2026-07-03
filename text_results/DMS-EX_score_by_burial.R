rm(list = ls(all.names = TRUE))

combined_buried_exposed <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/proteinGym_asymmetric_concat_rsa_group_similarity.tsv.gz',
                                      header=TRUE, sep="\t", stringsAsFactors = FALSE)

buried_tab <- combined_buried_exposed[combined_buried_exposed$rsa_group == "buried", ]
exposed_tab <- combined_buried_exposed[combined_buried_exposed$rsa_group == "exposed", ]

rownames(buried_tab) <- paste(buried_tab$ref_aa, buried_tab$mut_aa, sep = "_")
rownames(exposed_tab) <- paste(exposed_tab$ref_aa, exposed_tab$mut_aa, sep = "_")

combined_tab <- buried_tab[, "combined_buried_exposed_similiarity", drop = FALSE]
colnames(combined_tab) <- "dms.ex_buried"
combined_tab$dms.ex_exposed <- exposed_tab[rownames(combined_tab), "combined_buried_exposed_similiarity"]

mean(combined_tab$dms.ex_exposed / combined_tab$dms.ex_buried)
median(combined_tab$dms.ex_exposed / combined_tab$dms.ex_buried)
sd(combined_tab$dms.ex_exposed / combined_tab$dms.ex_buried)
wilcox.test(combined_tab$dms.ex_exposed, combined_tab$dms.ex_buried, paired = TRUE)
