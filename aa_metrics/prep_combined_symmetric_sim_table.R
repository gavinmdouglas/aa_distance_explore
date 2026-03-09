# Prepare table of all similarity (including others converted to similarity) measures
# that have been made symmetric. Doing so to make it easy to loop over all metrics.

tmp_in <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/AAOntology-PCs_ASA_Volume.tsv.gz",
                     header=TRUE, sep = '\t', stringsAsFactors = FALSE)

rownames(tmp_in) <- sapply(1:nrow(tmp_in), function(i) {
  paste0(c(tmp_in$aa1[i], tmp_in$aa2[i]), collapse = '_')
})

aa_sim_metrics <- tmp_in[, -c(1, 2, 3), drop = FALSE]

sim_files <- list.files("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent", pattern=".tsv.gz")

simfile_dir <- rep("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/", length(sim_files))

for (i in 1:length(sim_files)) {
  sim_file <- sim_files[i]
  dirpath <- simfile_dir[i]
  
  sim_mat <- read.table(paste0(dirpath, sim_file),
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  metric_name <- sub(".tsv.gz", "", sim_file)
  rownames(sim_mat) <- sapply(1:nrow(sim_mat), function(i) {
    paste0(c(sim_mat$aa1[i], sim_mat$aa2[i]), collapse = '_')
  })
  
  if (! all(rownames(aa_sim_metrics) %in% rownames(sim_mat))) {
    stop(paste0("Not all rownames of aa_sim_metrics are in ", sim))
  }
  
  # Two metrics are expected to be non-symmetrical, and need to be fixed.
  if ((metric_name == 'demask') || (metric_name == 'epstein')) {
    tmp_sim_mat = sim_mat
    for (i in 1:nrow(sim_mat)) {
      aa1 <- sim_mat$aa1[i]
      aa2 <- sim_mat$aa2[i]
      sim_val1 <- sim_mat[sim_mat$aa1 == aa1 & sim_mat$aa2 == aa2, 3]
      sim_val2 <- sim_mat[sim_mat$aa1 == aa2 & sim_mat$aa2 == aa1, 3]
      tmp_sim_mat[i, 3] <- mean(c(sim_val1, sim_val2))
      
    }
    sim_mat <- tmp_sim_mat
  }
  
  # Sanity check that it's symmetrical.
  for (i in 1:nrow(sim_mat)) {
    aa1 <- sim_mat$aa1[i]
    aa2 <- sim_mat$aa2[i]
    sim_val_1 <- sim_mat[i, 3]
    sim_val_2 <- sim_mat[which(sim_mat$aa1 == aa2 & sim_mat$aa2 == aa1), 3]
    if (length(sim_val_2) != 1) {
      stop(paste0("Missing symmetrical value for ", aa1, "_", aa2, " in ", sim_file))
    }
    if (sim_val_1 != sim_val_2) {
      stop(paste0("Asymmetrical values for ", aa1, "_", aa2, " in ", sim_file))
    }
  }
  
  aa_sim_metrics[, metric_name] <- sim_mat[rownames(aa_sim_metrics), 3]
}

write.table(x = aa_sim_metrics,
            file = gzfile("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/combined_symmetric_prepped_similarity_metrics.tsv.gz"),
            sep = "\t",
            col.names = NA,
            row.names = TRUE,
            quote=FALSE)
