rm(list = ls(all.names = TRUE))

# For checking how well two different AA similarity tables compare.

library(ggplot2)

file_paths <- c("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/DISTATIS_working/prepped_similarity_consistent/DISTATIS_Custom_EX_VHSE_Miyata_OtherAAOnt.tsv",
                "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/DISTATIS_working/prepped_similarity_consistent/DISTATIS_OLD/DISTATIS_Custom_EX_VHSE_Miyata_OtherAAOnt.tsv.gz",
                "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/ex.tsv.gz")

read_in_sim_tables <- function(file_paths) {
   
  sim_tables <- list()
  for (file_path in file_paths) {
    tab <- read.table(file_path, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    metric_name <- gsub('.tsv.gz', '', basename(file_path))
    metric_name <- gsub('.tsv', '', metric_name)
    
    if (metric_name %in% names(sim_tables)) {
      metric_name <- file_path
    }
    sim_tables[[metric_name]] <- tab
    rownames(sim_tables[[metric_name]]) <- paste0(sim_tables[[metric_name]]$aa1, '_', sim_tables[[metric_name]]$aa2)
  }
  
  # Get combined table with rows matching for values of aa1 and aa2
  combined_table <- sim_tables[[1]][, c('aa1', 'aa2')]
  for (metric_name in names(sim_tables)) {
    combined_table[, metric_name] <- sim_tables[[metric_name]][rownames(combined_table), 'similarity']
  }
  
  
  return(sim_tables)
  
}