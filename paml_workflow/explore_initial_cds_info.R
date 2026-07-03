rm(list = ls(all.names = TRUE))

# Explore initial CDS alignments that potentially could be used for PAML workflow, to see if any should be excluded.
mammal <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/OrthoMaM_seq_info.tsv",
                     header=TRUE, sep = '\t')

fly <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/Drosophila_all_cds_seq_info.tsv",
                     header=TRUE, sep = '\t')

Strep <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/Streptococcus_seq_info.tsv",
                     header=TRUE, sep = '\t')

# Restrict to proteins between 200-1000 codons.
# Have at least 20 sites with syn and non-syn variation (separately).

mammal <- mammal[mammal$alignment_length_codons >= 200 & mammal$alignment_length_codons <= 1000, ]
mammal <- mammal[which(mammal$codonsites_w_syn >= 20 & mammal$codonsites_w_nonsyn >= 20), ]

fly <- fly[fly$alignment_length_codons >= 200 & fly$alignment_length_codons <= 1000, ]
fly <- fly[which(fly$codonsites_w_syn >= 20 & fly$codonsites_w_nonsyn >= 20), ]

Strep <- Strep[Strep$alignment_length_codons >= 200 & Strep$alignment_length_codons <= 1000, ]
Strep <- Strep[which(Strep$codonsites_w_syn >= 20 & Strep$codonsites_w_nonsyn >= 20), ]

write.table(x = gsub('_NT_AL$', '', mammal$name),
            file = '/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/filt_ids/orthomam_filt_ids.txt',
            quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

write.table(x = gsub('.DNA.afa$', '', fly$name),
            file = '/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/filt_ids/drosophila_filt_ids.txt',
            quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

write.table(x = gsub('_final_mask_align_NT$', '', Strep$name),
            file = '/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/filt_ids/streptococcus_filt_ids.txt',
            quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
