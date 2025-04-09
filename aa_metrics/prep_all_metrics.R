rm(list = ls(all.names = TRUE))

source('~/Drive/research/aa_distance/aa_distance_explore/aa_metrics/aa_radical_vs_cons_define.R')

# Get easy to parse mappings of the effect of each AA sub., based on:
# radical vs conservative, Grantham's distance, BLOSUM62, and VTML200.

# First do radical vs conservative classifications.
# By charge.
charge_sub_raw <- list()
for (aa1 in names(charge_sub_change)) {
  for (aa2 in names(charge_sub_change[[aa1]])) {
    sub_type <- as.character(charge_sub_change[[aa1]][aa2])
    charge_sub_raw[[paste(aa1, aa2, sub_type, collapse = ',')]] <- data.frame(aa1=aa1, aa2=aa2, sub_type=sub_type)
  }
}

charge_sub_tab <- do.call(rbind, charge_sub_raw)
write.table(x = charge_sub_tab,
            file = gzfile("~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_RvC/RvC_charge.tsv.gz"),
            quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

# By polarity/volume.
polarity_volume_sub_raw <- list()
for (aa1 in names(polarity_volume_sub_change)) {
  for (aa2 in names(polarity_volume_sub_change[[aa1]])) {
    sub_type <- as.character(polarity_volume_sub_change[[aa1]][aa2])
    polarity_volume_sub_raw[[paste(aa1, aa2, sub_type, collapse = ',')]] <- data.frame(aa1=aa1, aa2=aa2, sub_type=sub_type)
  }
}

polarity_volume_sub_tab <- do.call(rbind, polarity_volume_sub_raw)
write.table(x = polarity_volume_sub_tab,
            file = gzfile("~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_RvC/RvC_polarity_volume_Zhang2000.tsv.gz"),
            quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

all_aa <- sort(unique(charge_sub_tab$aa1))

# Then convert distant matrices to long-format max-scaled (with small value added) normalized *similarity* matrices.
dist_file_to_minmax_inverted_simfile <- function(dist_file) {
  
  filebase <- gsub('.tsv.gz', '', basename(dist_file))
  
  intab <- read.table(dist_file, header = TRUE, row.names = 1, sep = '\t')
  
  if (nrow(intab) != 20 || ncol(intab) != 20) { stop('Input table does not have 20 rows and 20 columns!') }
  
  if (! identical(all_aa, sort(colnames(intab))) || ! identical(all_aa, sort(rownames(intab)))) {
    stop('AAs not all found in input table!')
  }
  
  if (any(is.na(intab))) { stop('NAs found in input table!') }
  
  intab <- intab[all_aa, all_aa]

  out_tab <- intab / (max(intab) * 1.01)
  out_tab <- 1 - out_tab
  
  out_tab$aa1 <- rownames(out_tab)
  out_tab <- reshape2::melt(out_tab, id.vars='aa1')
  colnames(out_tab) <- c('aa1', 'aa2', 'similarity')
  
  out_tab <- out_tab[which(out_tab$aa1 != out_tab$aa2), ]
  
  outfile <- paste0('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity/', filebase, '.maxscaled.invert.tsv.gz')
  
  write.table(x = out_tab, file = gzfile(outfile), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
  
  return('Finished.')
}


dist_files <- list.files('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/distances',
                                  pattern='.tsv.gz$', full.names = TRUE)

unique(sapply(dist_files, dist_file_to_minmax_inverted_simfile))


# Now do a similar procedure for the similarity matrices, but based on just max-scaling.
sim_files <- list.files('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities',
                        pattern='.tsv.gz$', full.names = TRUE)

sim_raw <- list()

for (sim_file in sim_files) {
  sim_base <- gsub('.tsv.gz', '', basename(sim_file))
  sim_raw[[sim_base]] <-  read.table(sim_file, header = TRUE, row.names = 1, sep = '\t')
}

sim_raw$csw[sim_raw$csw == 0] <- 0.01 * max(sim_raw$csw)

unique(sapply(names(sim_raw),
              function(x) { 
                
                intab <- sim_raw[[x]]
                
                if (nrow(intab) != 20 || ncol(intab) != 20) { stop('Input table does not have 20 rows and 20 columns!') }
                
                if (! identical(all_aa, sort(colnames(intab))) || ! identical(all_aa, sort(rownames(intab)))) {
                  stop('AAs not all found in input table!')
                }
                
                intab <- intab[all_aa, all_aa]

                out_tab <- intab / max(intab, na.rm = TRUE)
                
                out_tab$aa2 <- rownames(out_tab)
                out_tab <- reshape2::melt(out_tab, id.vars='aa2')
                colnames(out_tab) <- c('aa2', 'aa1', 'similarity')
                out_tab <- out_tab[, c('aa1', 'aa2', 'similarity')]
                
                out_tab <- out_tab[which(out_tab$aa1 != out_tab$aa2), ]
                
                outfile <- paste0('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity/', x, '.maxscaled.tsv.gz')
                
                write.table(x = out_tab, file = gzfile(outfile), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
                
                return('Finished.')
                
                }))


# Finally, read through substitution matrices (BLOMSUM62 and VTML200), and
# convert to probabilities from log-odds, and then write out.
# BLOSUM62
blosum62 <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/sub_matrices/blosum62.txt.gz',
                       header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)

blosum62 <- blosum62[all_aa, all_aa]

blosum62_odds <- exp(blosum62)
blosum62_prob <- blosum62_odds / (1 + blosum62_odds)

blosum62_prob$aa2 <- rownames(blosum62_prob)
blosum62_prob_long <- reshape2::melt(blosum62_prob, id.vars='aa2')
colnames(blosum62_prob_long) <- c('aa2', 'aa1', 'blosum62_prob')
blosum62_prob_long <- blosum62_prob_long[, c('aa1', 'aa2', 'blosum62_prob')]

blosum62_prob_long <- blosum62_prob_long[which(blosum62_prob_long$aa1 != blosum62_prob_long$aa2), ]

write.table(x = blosum62_prob_long,
                   file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity/blosum62.prob.tsv.gz'),
                   quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

# VTML200
vtml200 <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/sub_matrices/vtml200.txt.gz',
                       header = TRUE, sep = ' ', comment.char="#", row.names = 1, check.names = FALSE)

vtml200 <- vtml200[all_aa, all_aa]

vtml200_odds <- exp(vtml200)
vtml200_prob <- vtml200_odds / (1 + vtml200_odds)

vtml200_prob$aa2 <- rownames(vtml200_prob)
vtml200_prob_long <- reshape2::melt(vtml200_prob, id.vars='aa2')
colnames(vtml200_prob_long) <- c('aa2', 'aa1', 'vtml200_prob')
vtml200_prob_long <- vtml200_prob_long[, c('aa1', 'aa2', 'vtml200_prob')]

vtml200_prob_long <- vtml200_prob_long[which(vtml200_prob_long$aa1 != vtml200_prob_long$aa2), ]

write.table(x = vtml200_prob_long,
                   file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity/vtml200.prob.tsv.gz'),
                   quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
