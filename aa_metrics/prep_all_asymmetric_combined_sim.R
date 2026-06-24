rm(list = ls(all.names = TRUE))

# Prepare table of all asymmetric metrics (as similarities).

proteinGym_metrics <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/proteinGym_asymmetric_similarity.tsv.gz",
                                 header=TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

ex <-  read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/ex_asymmetric.tsv.gz",
                  header=TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Rows are WT
# Columns are variant
demask <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/demask.txt",
                     header=TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
rownames(demask) <- colnames(demask)

epstein <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/distances/epstein.tsv.gz",
                      header=TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

# Swap in reverse sub score for missing EX values.
for (i in which(is.na(ex$asymmetric_exchange))) {
    na_row <- ex[i, ]
    aa1 <- na_row$aa1
    aa2 <- na_row$aa2
    ex[i, 'asymmetric_exchange'] <- ex[ex$aa1 == aa2 & ex$aa2 == aa1, 'asymmetric_exchange']
}

asymmetric_vals <- proteinGym_metrics
asymmetric_vals$ex <- NA
asymmetric_vals$demask <- NA
asymmetric_vals$epstein <- NA

for (i in 1:nrow(asymmetric_vals)) {
  aa1 <- asymmetric_vals$ref_aa[i]
  aa2 <- asymmetric_vals$mut_aa[i]
  
  ex_val <- ex$asymmetric_exchange[ex$aa1 == aa1 & ex$aa2 == aa2]
  if (length(ex_val) == 1) {
    asymmetric_vals[i, 'ex'] <- ex_val
  } else {
    stop(paste("Missing EX value for", aa1, aa2))
  }
  
  demask_val <- demask[aa1, aa2]
  if (length(demask_val) == 1) {
    asymmetric_vals[i, 'demask'] <- demask_val
  } else {
    stop(paste("Missing DeMaSk value for", aa1, aa2))
  }
  
  epstein_val <- epstein[aa1, aa2]
  if (length(epstein_val) == 1) {
    asymmetric_vals[i, 'epstein'] <- epstein_val
  } else {
    stop(paste("Missing Epstein value for", aa1, aa2))
  }
  
}

asymmetric_vals$epstein <- 1 - asymmetric_vals$epstein

# Then transform all columns, as before:

transform_to_range <- function(numbers, min_val = 0.01, max_val = 0.99) {
  input_min <- min(numbers)
  input_max <- max(numbers)
  return(((numbers - input_min) / (input_max - input_min)) * (max_val - min_val) + min_val)
}

asymmetric_vals$ex <- transform_to_range(asymmetric_vals$ex)
asymmetric_vals$demask <- transform_to_range(asymmetric_vals$demask)
asymmetric_vals$epstein <- transform_to_range(asymmetric_vals$epstein)

colnames(asymmetric_vals)[colnames(asymmetric_vals) == "weighted"] <- "dms.ex_weighted"
colnames(asymmetric_vals)[colnames(asymmetric_vals) == "exposed"] <- "dms.ex_exposed"
colnames(asymmetric_vals)[colnames(asymmetric_vals) == "buried"] <- "dms.ex_buried"

# Also get average values across differing measures.
# DISTATIS not used in this case as these are not symmetric matrices.
# Note these were transformed to same range as others afterwards.
asymmetric_vals$dms.ex_ex <- transform_to_range(rowMeans(asymmetric_vals[, c("dms.ex_weighted", "ex")], na.rm = TRUE))
asymmetric_vals$dms.ex_demask <- transform_to_range(rowMeans(asymmetric_vals[, c("dms.ex_weighted", "demask")], na.rm = TRUE))
asymmetric_vals$ex_demask <- transform_to_range(rowMeans(asymmetric_vals[, c("ex", "demask")], na.rm = TRUE))
asymmetric_vals$dms.ex_ex_demask <- transform_to_range(rowMeans(asymmetric_vals[, c("dms.ex_weighted", "ex", "demask")], na.rm = TRUE))

write.table(x = asymmetric_vals,
            file = gzfile("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/asymmetric_similarity_consistent.tsv.gz"),
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote=FALSE)
