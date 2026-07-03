rm(list = ls(all.names = TRUE))

RvC_files <- list.files("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_RvC",
                        pattern = "RvC.*\\.tsv\\.gz", full.names = TRUE)

for (RvC_file in RvC_files) {
  metric <- sub('.tsv.gz', '', basename(RvC_file))
  intab <- read.table(RvC_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  intab$similarity <- NA
  intab[which(intab$sub_type == 'Conservative'), 'similarity'] <- 0.99
  intab[which(intab$sub_type == 'Radical'), 'similarity'] <- 0.01
  intab <- intab[, c('aa1', 'aa2', 'similarity')]
  
  # Check for NA values.
  if (any(is.na(intab$similarity))) {
    stop(paste0("NA values found in metric: ", metric, ". Please check the input data."))
  }
  
  if (nrow(intab) != 380) {
    stop(paste0("Unexpected number of rows in metric: ", metric, ". Expected 380, found ", nrow(intab), "."))
  }
  
  
  outfile <- paste0('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_RvC/prepped_similarity_consistent/',
                    metric, '.tsv')
  write.table(x = intab, file = outfile, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  
}
