rm(list = ls(all.names = TRUE))

# Sanity checks on prepped tables.
prepped_tabs <- list()

for (tab in list.files('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/prepped',
                       pattern='.tsv.gz$', full.names = TRUE)) {
  category <- gsub('\\..*$', '', basename(tab))
  prepped_tabs[[category]] <- read.table(tab, header = TRUE, sep = '\t')
}

# Look at 10 random samples of RvC of each type.
prepped_tabs$RvC_charge[sample(nrow(prepped_tabs$RvC_charge), 10), ]
prepped_tabs$RvC_polarity_volume[sample(nrow(prepped_tabs$RvC_polarity_volume), 10), ]


# Then look at scatterplots of original metrics vs. the corresponding prepped similarities.
raw_files <- c(list.files('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/distances',
                            pattern='.tsv.gz$', full.names = TRUE),
               list.files('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/similarities',
                            pattern='.tsv.gz$', full.names = TRUE),
               list.files('~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/sub_matrices',
                            pattern='.txt.gz$', full.names = TRUE))

for (raw_file in raw_files) {
  
  base_name <- gsub('\\..*$', '', basename(raw_file))
  
  if (base_name == "vtml200") {
    raw_matrix <- read.table(raw_file, header = TRUE, sep = ' ', comment.char="#", row.names = 1, check.names = FALSE)
  } else if (base_name == "blosum62") {
    raw_matrix <- read.table(raw_file, header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)
  } else {
    raw_matrix <- read.table(raw_file, header = TRUE, sep = '\t', row.names = 1)
  }

  comparison <- prepped_tabs[[base_name]]
  
  colnames(comparison) <- c('aa1', 'aa2', 'similarity')
  
  comparison$raw <- NA
  
  for (i in 1:nrow(comparison)) {
    comparison$raw[i] <- raw_matrix[which(rownames(raw_matrix) == comparison$aa1[i]), which(colnames(raw_matrix) == comparison$aa2[i])]
  }
  
  plot(comparison$raw, comparison$similarity, main = base_name)
  
}
