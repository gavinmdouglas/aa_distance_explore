rm(list = ls(all.names = TRUE))

# Assess how many different species of Streptococcus are found in dataset
# (which happened to be generated for a different project, but I think would
# be good to include for this work as well).

# Note that this was based on a larger project, and skani was re-run just on the subset of genomes of interest, as a sanity check.
# Panta was also re-run specifically for this subset, to ensure including a wider set of genomes did not influence the result.

strep_genomes <- read.table('~/Strep_genomes.txt',
                            header=FALSE, sep='\t', stringsAsFactors = FALSE)$V1

nonfrequent_species_info <- read.table('/Users/gavin/Drive/research/intergenic_evolution/data/ncbi_download/ncbi_bacteria_archaea_nonmag_dataset_genus_w_10_nonfrequent_species.tsv.gz',
                                       header=TRUE, sep='\t', stringsAsFactors = FALSE, comment.char = '', quote = '')

frequent_species_info <- read.table('/Users/gavin/Drive/research/intergenic_evolution/data/ncbi_download/ncbi_bacteria_archaea_nonmag_dataset_genus_w_10_frequent_species_ALL.tsv.gz',
                                    header=TRUE, sep='\t', stringsAsFactors = FALSE, comment.char = '', quote = '')

species_info <- rbind(nonfrequent_species_info, frequent_species_info)

strep_info <- species_info[which(species_info$accession %in% strep_genomes), ]

unique_species <- unique(strep_info$species)

raw <- list()

for (sp in unique_species) {
 
  sp_subset <- strep_info[which(strep_info$species == sp), ]
  
  chrom_i <- which(sp_subset$assembly_level == 'Chromosome')
  if (length(chrom_i) > 0) {
    raw[[sp]] <- sp_subset[chrom_i[1], , drop = FALSE]
    next
  }
  
  complete_i <- which(sp_subset$assembly_level == 'Complete Genome')
  if (length(complete_i) > 0) {
    raw[[sp]] <- sp_subset[complete_i[1], , drop = FALSE]
    next
  }
  
  scaffold_i <- which(sp_subset$assembly_level == 'Scaffold')
  if (length(scaffold_i) > 0) {
    raw[[sp]] <- sp_subset[scaffold_i[1], , drop = FALSE]
    next
  }
  
  contig_i <- which(sp_subset$assembly_level == 'Contig')
  if (length(contig_i) > 0) {
    raw[[sp]] <- sp_subset[contig_i[1], , drop = FALSE]
    next
  }
  
  print(sp_subset)
  stop('No assembly level found for strep species: ', sp, ' in dataset. Please check the genome type of this species in the dataset.')
}

strep_species_rep <- do.call(rbind, raw)

# De-rep based on ANI as well. Require all genomes to have < 95% pairwise ANI.
strep_skani <- read.table('/Users/gavin/Drive/research/intergenic_evolution/data/skani_output/merged_freq_nonfreq/Streptococcus.txt.gz',
                                header=TRUE, sep='\t', stringsAsFactors = FALSE, quote='', comment.char='')

strep_skani$genome1 <- sapply(strep_skani$Ref_file,
                                    function(x) {
                                      tmp <- strsplit(basename(x), '_')[[1]]
                                      paste(tmp[1], tmp[2], sep = '_')
                                      })

strep_skani$genome2 <- sapply(strep_skani$Query_file,
                                    function(x) {
                                      tmp <- strsplit(basename(x), '_')[[1]]
                                      paste(tmp[1], tmp[2], sep = '_')
                                    })

unique_genomes <- strep_species_rep$accession

strep_skani <- strep_skani[which(strep_skani$genome1 %in% unique_genomes & 
                                            strep_skani$genome2 %in% unique_genomes), ]

strep_skani_high_ANI <- strep_skani[which(strep_skani$ANI >= 95), ]

genomes_to_exclude <- character()
while (nrow(strep_skani_high_ANI) > 0) {
  high_ANI_genome_tallies <- table(c(strep_skani_high_ANI$genome1, strep_skani_high_ANI$genome2))
  max_tally <- max(high_ANI_genome_tallies)
  genome_to_rm <- names(high_ANI_genome_tallies)[which(high_ANI_genome_tallies == max_tally)][1]
  genomes_to_exclude <- c(genomes_to_exclude, genome_to_rm)
  strep_skani_high_ANI <- strep_skani_high_ANI[which(strep_skani_high_ANI$genome1 != genome_to_rm & 
                                                     strep_skani_high_ANI$genome2 != genome_to_rm), ]
}

strep_species_derep <- strep_species_rep[which(!strep_species_rep$accession %in% genomes_to_exclude), ]

# Realized that 71 genomes would be plenty for this analysis, so just kept complete genomes.

strep_species_derep <- strep_species_derep[which(strep_species_derep$assembly_level == 'Complete Genome' | 
                                                 strep_species_derep$assembly_level == 'Chromosome'), ]

# Sanity check afterwards:
unique_genomes <- strep_species_derep$accession

strep_skani <- strep_skani[which(strep_skani$genome1 %in% unique_genomes & 
                                   strep_skani$genome2 %in% unique_genomes), ]

strep_skani_high_ANI <- strep_skani[which(strep_skani$ANI >= 95), ]

# Write out genome set to analyze.
Strep_differing_species_accesions <- unique_genomes
write.table(x = Strep_differing_species_accesions,
            file = '/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_running/Streptococcus/species_accessions.txt',
            sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
