rm(list = ls(all.names = TRUE))

# Polarize mutations in human table and output final table.
human_tab <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_snvs_w_prefs_dssp_anc_aa.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# Polarize mutations to get unfolded AF spectrum.
human_tab_forward <- human_tab[which(human_tab$ref_aa == human_tab$unambig_ancestral_AA), ]
human_tab_forward$anc_aa <- human_tab_forward$ref_aa
human_tab_forward$derived_aa <- human_tab_forward$alt_aa
human_tab_forward$AF <- human_tab_forward$variant_count / human_tab_forward$sample_count

human_tab_backward <- human_tab[which(human_tab$alt_aa == human_tab$unambig_ancestral_AA), ]
human_tab_backward$anc_aa <- human_tab_backward$alt_aa
human_tab_backward$variant_count <- human_tab_backward$sample_count - human_tab_backward$variant_count
human_tab_backward$derived_aa <- human_tab_backward$ref_aa
human_tab_backward$AF <- human_tab_backward$variant_count / human_tab_backward$sample_count


human_tab <- rbind(human_tab_forward, human_tab_backward)

human_tab <- human_tab[human_tab$variant_count > 0, ]

# Need to flip RaSP direction:
human_tab$rasp <- human_tab$rasp * -1

# Ignore mutations in start codon position.
human_tab <- human_tab[which(human_tab$pos != 0), ]

human_tab <- human_tab[human_tab$sample_count >= 1400000, ]

human_tab$aa_combo <- apply(human_tab[, c('anc_aa', 'derived_aa')], 1, function(x) {
  paste(x, collapse = '_')
})

# Basic summary statistics:
dim(human_tab)
length(unique(human_tab$ensembl_protein))
min(human_tab$sample_count)

write.table(human_tab,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_snvs_w_prefs_dssp_anc_aa_prep.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
