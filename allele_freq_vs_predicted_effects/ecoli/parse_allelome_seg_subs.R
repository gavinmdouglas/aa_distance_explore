rm(list = ls(all.names = TRUE))

library(stringr)

# Parse table of segregating E. coli amino acid variants from the Allelome
# supplementary data.

# As a reminder, the goal is to get a table similar to this:
# protein pos consensus_codon consensus_aa seg_codon seg_aa AC     N           AF
# 1 PF00004-GA4805AA_00441   0             CTC            L       ATC      I  1 60899 1.642063e-05
# 2 PF00004-GA4805AA_00441   5             CCG            P       TCG      S  1 60899 1.642063e-05
# 3 PF00004-GA4805AA_00441  16             TCC            S       CCC      P  1 60899 1.642063e-05

codon_table <- read.table("~/Drive/code_maintenance/codon_tables/genetic_codes/code11.tsv",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names=1)
rownames(codon_table) <- gsub("T", "U", rownames(codon_table))

focal_id_map <- read.table("/Users/gavin/Drive/research/aa_distance/data/Ecoli_focal_seq_id_map.txt",
                           sep = " ", stringsAsFactors = FALSE, header=FALSE)
colnames(focal_id_map) <- c("Blattner_num", "uniprot_id")
rownames(focal_id_map) <- focal_id_map$Blattner_num

allelome_files <- list.files("~/Drive/research/aa_distance/data/allelome_data", pattern = ".csv.gz$", full.names = TRUE)

parse_codon_details <- function(detail_str) {
  # Parse item like {'AAG': 2614.0, 'AAA': 1.0} into named vector
  matches <- stringr::str_match_all(detail_str, "'?([A-Z*]+)'?:\\s*([0-9.]+)")[[1]]
  counts <- as.numeric(matches[, 3])
  names(counts) <- matches[, 2]
  counts
}

raw <- list()

for (alleleome_file in allelome_files) {
 
  intab <- read.table(alleleome_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  Blattner_num <- gsub("_dfz.csv.gz$", "", basename(alleleome_file))
  
  if (! Blattner_num %in% focal_id_map$Blattner_num) { next }
  
  uniprot_id <- focal_id_map[Blattner_num, "uniprot_id"]

  gene_pos <- c()
  gene_consensus_codon <- c()
  gene_consensus_aa <- c()
  gene_alt_codon <- c()
  gene_alt_aa <- c()
  gene_alt_count <- c()
  gene_consensus_count <- c()
  gene_af <- c()
    
  for (i in 1:nrow(intab)) {
    
    row <- intab[i, ]
    consensus_codon <- row$consensus_codon_seq
    consensus_aa <- row$consensus_AAseq
    detail_str <- row$consensus_codon_seq_details
    
    codon_counts <- parse_codon_details(detail_str)
    
    # Only keep valid codons.
    codon_counts <- codon_counts[names(codon_counts) %in% rownames(codon_table)]
    
    # No variants, so just return.
    if (length(codon_counts) < 2) { next }

    non_consensus <- codon_counts[names(codon_counts) != consensus_codon]
    
    # Translate and keep only those differing in AA from consensus
    non_consensus_aa <- codon_table[names(non_consensus), 'amino_acid']
    
    nonsyn <- non_consensus[non_consensus_aa != consensus_aa]
    nonsyn_aa <- non_consensus_aa[non_consensus_aa != consensus_aa]
    
    seg_syn_count <- sum(non_consensus[which(non_consensus_aa == consensus_aa)])

    if (length(nonsyn) == 0) { next }
    
    max_codon_i <- which.max(nonsyn)
    top_codon <- names(nonsyn)[max_codon_i]
    
    # Only consider non-synonymous codons that differ by one mutation from the consensus.
    num_mutations <- sum(strsplit(top_codon, "")[[1]] != strsplit(consensus_codon, "")[[1]])
    if (num_mutations > 1) { next }

    top_aa <- nonsyn_aa[max_codon_i]
    top_count <- as.integer(nonsyn[max_codon_i])
    
    # AF excludes all other non-consensus codons
    consensus_count <- codon_counts[consensus_codon]
    total <- consensus_count + top_count
    af <- top_count / total
    
    gene_pos <- c(gene_pos, row$X)
    gene_consensus_codon <- c(gene_consensus_codon, consensus_codon)
    gene_consensus_aa <- c(gene_consensus_aa, consensus_aa)
    gene_alt_codon <- c(gene_alt_codon, top_codon)
    gene_alt_aa <- c(gene_alt_aa, top_aa)
    gene_alt_count <- c(gene_alt_count, top_count)
    gene_consensus_count <- c(gene_consensus_count, consensus_count)
    gene_af <- c(gene_af, af)

  }
  
  if (length(gene_pos) == 0) { next }
  
  raw[[Blattner_num]] <- data.frame(Blattner_num = Blattner_num,
                                    uniprot_id = uniprot_id,
                                    pos = gene_pos,
                                    consensus_codon = gene_consensus_codon,
                                    consensus_aa = gene_consensus_aa,
                                    alt_codon = gene_alt_codon,
                                    alt_aa = gene_alt_aa,
                                    alt_count = gene_alt_count,
                                    consensus_count = gene_consensus_count,
                                    AF = gene_af)
  
}

results_df <- do.call(rbind, raw)
rownames(results_df) <- NULL

# Remove first position in proteins.
results_df <- results_df[results_df$pos > 0, ]

write.table(results_df,
            gzfile("~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/parsed_allelome_variants.tsv.gz"),
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
