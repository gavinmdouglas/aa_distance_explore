rm(list = ls(all.names = TRUE))

# Polarize mutations in E. coli table and output final table.
ecoli_tab <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/parsed_allelome_variants_w_prefs_dssp_w_anc.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# Polarize mutations to get unfolded AF spectrum.
ecoli_tab_forward <- ecoli_tab[which(ecoli_tab$consensus_aa == ecoli_tab$Likely_Anc_AA), ]
ecoli_tab_forward$anc_aa <- ecoli_tab_forward$consensus_aa
ecoli_tab_forward$anc_count <- ecoli_tab_forward$consensus_count
ecoli_tab_forward$derived_aa <- ecoli_tab_forward$alt_aa
ecoli_tab_forward$derived_count <- ecoli_tab_forward$alt_count

ecoli_tab_backward <- ecoli_tab[which(ecoli_tab$alt_aa == ecoli_tab$Likely_Anc_AA), ]
ecoli_tab_backward$anc_aa <- ecoli_tab_backward$alt_aa
ecoli_tab_backward$anc_count <- ecoli_tab_backward$alt_count
ecoli_tab_backward$derived_aa <- ecoli_tab_backward$consensus_aa
ecoli_tab_backward$derived_count <- ecoli_tab_backward$consensus_count

ecoli_tab_backward$AF <- 1 - ecoli_tab_backward$AF

ecoli_tab <- rbind(ecoli_tab_forward, ecoli_tab_backward)

ecoli_tab$count <- ecoli_tab$anc_count + ecoli_tab$derived_count

# Only consider sites with at least 2000 unambiguous sites.
sub_tab_filt <- ecoli_tab[which(ecoli_tab$count >= 2000), ]

# Ignore mutations in start codon position.
sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$pos != 0), ]

# Ignore singletons
sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$derived_count > 1), ]
sub_tab_filt <- sub_tab_filt[which(sub_tab_filt$anc_count > 1), ]

sub_tab_filt <- sub_tab_filt[rowSums(is.na(sub_tab_filt)) == 0, ]

# Create AA combinations
sub_tab_filt$aa_combo <- apply(sub_tab_filt[, c('anc_aa', 'derived_aa')], 1, function(x) {
  paste(x, collapse = '_')
})

# Basic summary statistics:
dim(sub_tab_filt)
length(unique(sub_tab_filt$protein))
min(sub_tab_filt$count)

write.table(sub_tab_filt,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/parsed_allelome_variants_w_prefs_dssp_w_anc_prep.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
