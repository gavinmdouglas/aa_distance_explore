rm(list = ls(all.names = TRUE))

library('DistatisR')

check_symmetrical <- function(mat) {
  if (! all(mat == t(mat))) {
    FALSE
  } else{
    TRUE
  }
}

amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# Grantham
grantham_dist <- read.table("~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/distances/grantham_orig.tsv.gz",
                            header = TRUE, row.names = 1, sep = "\t")
grantham_dist <- grantham_dist / max(grantham_dist)
grantham_dist <- grantham_dist[amino_acids, amino_acids]
check_symmetrical(grantham_dist)

# Sneath
sneath_dist <- read.table("~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/distances/sneath.tsv.gz",
                            header = TRUE, row.names = 1, sep = "\t")
sneath_dist <- sneath_dist / max(sneath_dist)
sneath_dist <- sneath_dist[amino_acids, amino_acids]
check_symmetrical(sneath_dist)

# Epstein
epstein_dist <- read.table("~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/distances/epstein.tsv.gz",
                            header = TRUE, row.names = 1, sep = "\t")
epstein_dist <- epstein_dist / max(epstein_dist)
epstein_dist <- epstein_dist[amino_acids, amino_acids]
check_symmetrical(epstein_dist)

# Take average of off-diagonal values, to make the matrix symmetrical:
for (i in 1:(ncol(epstein_dist) - 1)) {
  for(j in i:ncol(epstein_dist)) {
    if (i != j) {
      mean_val <- (epstein_dist[i, j] + epstein_dist[j, i]) / 2
      epstein_dist[i, j] <- mean_val
      epstein_dist[j, i] <- mean_val 
    }
  }
}
check_symmetrical(epstein_dist)

# Atchley
atchley_dist <- read.table("~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/distances/atchley_minmax.tsv.gz",
                          header = TRUE, row.names = 1, sep = "\t")
atchley_dist <- atchley_dist / max(atchley_dist)
atchley_dist <- atchley_dist[amino_acids, amino_acids]
check_symmetrical(atchley_dist)

# EMPAR
empar_sim <- read.table("~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/similarities/empar.tsv.gz",
                          header = TRUE, row.names = 1, sep = "\t")
empar_sim <- empar_sim / max(empar_sim)
empar_dist <- 1 - empar_sim
empar_dist <- empar_dist[amino_acids, amino_acids]
check_symmetrical(empar_dist)

# Radical vs. conservative
RvC_polvol_breakdown <- read.table("~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/prepped_RvC/RvC_polarity_volume.tsv.gz",
                                   header = TRUE, sep = "\t")
RvC_polvol_dummy <- matrix(0, nrow=20, ncol=20)
for (i in 1:nrow(RvC_polvol_breakdown)) {
  if (RvC_polvol_breakdown[i, 3] == 'Radical') {
    aa1_index <- which(amino_acids == RvC_polvol_breakdown[i, 1])
    aa2_index <- which(amino_acids == RvC_polvol_breakdown[i, 2])
    RvC_polvol_dummy[aa1_index, aa2_index] <- 1
  }
}
check_symmetrical(RvC_polvol_dummy)


RvC_charge_breakdown <- read.table("~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/prepped_RvC/RvC_charge.tsv.gz",
                                   header = TRUE, sep = "\t")
RvC_charge_dummy <- matrix(0, nrow=20, ncol=20)
for (i in 1:nrow(RvC_charge_breakdown)) {
  if (RvC_charge_breakdown[i, 3] == 'Radical') {
    aa1_index <- which(amino_acids == RvC_charge_breakdown[i, 1])
    aa2_index <- which(amino_acids == RvC_charge_breakdown[i, 2])
    RvC_charge_dummy[aa1_index, aa2_index] <- 1
  }
}
check_symmetrical(RvC_charge_dummy)

all_dist_matrices <- array(0, dim = c(20, 20, 7))
all_dist_matrices[, , 1] <- as.matrix(grantham_dist)
all_dist_matrices[, , 2] <- as.matrix(sneath_dist)
all_dist_matrices[, , 3] <- as.matrix(epstein_dist)
all_dist_matrices[, , 4] <- as.matrix(atchley_dist)
all_dist_matrices[, , 5] <- as.matrix(empar_dist)
all_dist_matrices[, , 6] <- as.matrix(RvC_polvol_dummy)
all_dist_matrices[, , 7] <- as.matrix(RvC_charge_dummy)

distasis_output <- DistatisR::distatis(all_dist_matrices, nfact2keep = 20)

# Save eigenvalues and factor scores.
write.table(x = distasis_output$res4Splus$eigValues,
            file = "~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/DISTASIS_working/eigenvalues.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

factor_scores <- data.frame(distasis_output$res4Splus$F[, 1:19])
orig_col <- colnames(factor_scores)
factor_scores$amino_acids <- amino_acids
factor_scores <- factor_scores[, c("amino_acids", orig_col)]
write.table(x = factor_scores,
            file = "~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/DISTASIS_working/distatis_factor_scores.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

# Also, write out factor scores with charge categories,
# as this information was not captured by the other dissimilarity matrices,
# nor in the factors. Will be tested as well.
factor_scores_w_charge <- factor_scores
factor_scores_w_charge$charge <- 0
factor_scores_w_charge$charge[which(amino_acids %in% c('R', 'H', 'K'))] <- 0.1
factor_scores_w_charge$charge[which(amino_acids %in% c('D', 'E'))] <- -0.1
factor_scores_w_charge <- factor_scores_w_charge[, c("amino_acids", "charge", orig_col)]
write.table(x = factor_scores_w_charge,
            file = "~/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/DISTASIS_working/distatis_factor_scores_w_charge.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

# Also save distance matrices for use in related analyses.
write_tab_w_row_name <- function(df, file, new_col_name, sep_char = "\t") {
  orig_col <- colnames(df)
  df[, new_col_name] <- rownames(df)
  df <- df[, c(new_col_name, orig_col)]
  write.table(x = df,
              file = file,
              quote = FALSE, col.names = TRUE, row.names = FALSE, sep = sep_char)
  
}

write_tab_w_row_name(df=grantham_dist,
                     file="/Users/gavin/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/DISTASIS_working/prepped_dissimilarity_matrices/grantham.tsv",
                     new_col_name="amino_acid")
write_tab_w_row_name(df=sneath_dist,
                     file="/Users/gavin/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/DISTASIS_working/prepped_dissimilarity_matrices/sneath.tsv",
                     new_col_name="amino_acid")
write_tab_w_row_name(df=epstein_dist,
                     file="/Users/gavin/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/DISTASIS_working/prepped_dissimilarity_matrices/epstein.tsv",
                     new_col_name="amino_acid")
write_tab_w_row_name(df=atchley_dist,
                     file="/Users/gavin/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/DISTASIS_working/prepped_dissimilarity_matrices/atchley.tsv",
                     new_col_name="amino_acid")
write_tab_w_row_name(df=empar_dist,
                     file="/Users/gavin/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/DISTASIS_working/prepped_dissimilarity_matrices/empar.tsv",
                     new_col_name="amino_acid")

RvC_polvol_dummy <- data.frame(RvC_polvol_dummy)
colnames(RvC_polvol_dummy) <- amino_acids
rownames(RvC_polvol_dummy) <- amino_acids
write_tab_w_row_name(df=RvC_polvol_dummy,
                     file="/Users/gavin/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/DISTASIS_working/prepped_dissimilarity_matrices/RvC_polvol.tsv",
                     new_col_name="amino_acid")

RvC_charge_dummy <- data.frame(RvC_charge_dummy)
colnames(RvC_charge_dummy) <- amino_acids
rownames(RvC_charge_dummy) <- amino_acids
write_tab_w_row_name(df=RvC_charge_dummy,
                     file="/Users/gavin/Drive/ncsu/aa_selection/aa_selection_zenodo/aa_metrics/DISTASIS_working/prepped_dissimilarity_matrices/RvC_charge.tsv",
                     new_col_name="amino_acid")
