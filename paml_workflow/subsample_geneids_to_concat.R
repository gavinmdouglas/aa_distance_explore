rm(list = ls(all.names = TRUE))

# Get 20 random sets of 12 genes for each lineage, which will be concatenated into single alignments.

parse_samples <- function(geneids, outdir) {
  raw_sampling <- sample(geneids, size = 240, replace=FALSE)
  start_i <- 1
  end_i <- 12
  for (i in 1:20) {
    id_samples <- raw_sampling[start_i:end_i]
    write.table(x=id_samples, file=paste0(outdir, '/subsample', as.character(i), '.txt'),
                col.names = FALSE, row.names = FALSE, quote = FALSE)
    start_i <- start_i + 12
    end_i <- end_i + 12
  }
}

set.seed(18)

mammal <- as.character(read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/filt_ids/orthomam_filt_ids.txt",
                                  stringsAsFactors = FALSE, header=FALSE)$V1)
parse_samples(mammal, outdir = "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/genes_to_concat/Mammal")

fly <- as.character(read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/filt_ids/drosophila_filt_ids.txt",
                                  stringsAsFactors = FALSE, header=FALSE)$V1)
parse_samples(fly, outdir = "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/genes_to_concat/Drosophila")

Strep <- as.character(read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/filt_ids/streptococcus_filt_ids.txt",
                               stringsAsFactors = FALSE, header=FALSE)$V1)
parse_samples(Strep, outdir = "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/seq_info/genes_to_concat/Strep/")
