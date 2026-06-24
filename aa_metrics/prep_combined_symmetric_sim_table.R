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

# Remove mean-based DMS-EX mesaure to avoid confusion.
aa_sim_metrics <- aa_sim_metrics[, -which(colnames(aa_sim_metrics) == "proteinGym_rsa_mean_custom")]

write.table(x = aa_sim_metrics,
            file = gzfile("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/combined_symmetric_prepped_similarity_metrics.tsv.gz"),
            sep = "\t",
            col.names = NA,
            row.names = TRUE,
            quote=FALSE)

# Also create a table with re-ordered columns and with clean names.
# Output as DISTANCE, as this is how the measures are generally discussed.
metrics_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

colnames(aa_sim_metrics) <- sub('\\.prob$', '', colnames(aa_sim_metrics))
colnames(aa_sim_metrics) <- sub('AAOntology-', 'AAOntology.', colnames(aa_sim_metrics))

colnames(aa_sim_metrics) <- metrics_map[colnames(aa_sim_metrics), 'Clean']

all_nonfocal_DISTATIS <- grep("\\+", colnames(aa_sim_metrics), value = TRUE)
AA_ontology <- grep("AA-Ont.", colnames(aa_sim_metrics), value = TRUE)
AA_ontology <- setdiff(AA_ontology, all_nonfocal_DISTATIS)

experimental_measure <- c("DEX", "EX", "DMS-EX", "DeMaSk")
sub_measure <- c("BLOSUM62", "VTML200", "Xia", "JTT", "WAG", "LG")
RvC_measure <- grep("RvC", colnames(aa_sim_metrics), value=TRUE)

grantham_other_measures <- grep("Grantham", colnames(aa_sim_metrics), value=TRUE)
grantham_other_measures <- grantham_other_measures[-which(grantham_other_measures == "Grantham (orig.)")]
grantham_other_measures <- setdiff(grantham_other_measures, all_nonfocal_DISTATIS)

miyata_other_measures <- grep("Miyata", colnames(aa_sim_metrics), value=TRUE)
miyata_other_measures <- miyata_other_measures[-which(miyata_other_measures == "Miyata (orig.)")]
miyata_other_measures <- setdiff(miyata_other_measures, all_nonfocal_DISTATIS)

physiochemical_measures <- c("Grantham (orig.)", "Miyata (orig.)", "Epstein", "Sneath", "EMPAR", "CSW", "Atchley",
                             "Cruciani", "FASAGI", "Kidera", "Venkatarajan", "VHSE", "zScales")

measure_order <- c(experimental_measure,
                   sub_measure,
                   RvC_measure,
                   physiochemical_measures,
                   AA_ontology,
                   grantham_other_measures,
                   miyata_other_measures,
                   all_nonfocal_DISTATIS)

measure_order[which(duplicated(measure_order))]

setdiff(measure_order, colnames(aa_sim_metrics))
setdiff(colnames(aa_sim_metrics), measure_order)

aa_sim_metrics <- aa_sim_metrics[, measure_order]

aa_dist_metrics <- 1 - aa_sim_metrics

write.table(x = aa_dist_metrics,
            file = gzfile("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/all_distance_measures_symmetric.tsv.gz"),
            sep = "\t",
            col.names = NA,
            row.names = TRUE,
            quote=FALSE)
