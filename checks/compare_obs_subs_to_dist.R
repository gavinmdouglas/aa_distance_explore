rm(list = ls(all.names = TRUE))

d_test <- read.table('~/drosophila_test.tsv', header=TRUE, row.names = 1, sep = '\t')
m_test <- read.table('~/mammal_test.tsv', header=TRUE, row.names = 1, sep = '\t')
s_test <- read.table('~/Strep_test.tsv', header=TRUE, row.names = 1, sep = '\t')

info <- data.frame(matrix(, nrow=0, ncol=5))
colnames(info) <- c('aa1', 'aa2', 'd', 'm', 's')

all_aa <- rownames(d_test)

for (aa1 in all_aa[-length(all_aa)]) {
 
  for (aa2 in all_aa[(which(all_aa == aa1) + 1):length(all_aa)]) {
    d_val <- mean(c(d_test[aa1, aa2], d_test[aa2, aa1]))
    m_val <- mean(c(m_test[aa1, aa2], m_test[aa2, aa1]))
    s_val <- mean(c(s_test[aa1, aa2], s_test[aa2, aa1]))
    
    info <- rbind(info, data.frame(aa1=aa1, aa2=aa2, d=d_val, m=m_val, s=s_val))

  }
  
  rownames(info) <- sapply(1:nrow(info), function(x) {
    paste0(info$aa1[x], '_', info$aa2[x])
  })
}

file_paths <- c("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/DISTATIS_working/prepped_similarity_consistent/DISTATIS_Custom_EX.tsv.gz",
                "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/ex.tsv.gz",
                "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent/grantham_orig.tsv.gz")
                
other_info <- read_in_sim_tables(file_paths)

info$custom_ex <- other_info$DISTATIS_Custom_EX[rownames(info), 'similarity']
info$ex <- other_info$ex[rownames(info), 'similarity']
info$grantham <- other_info$grantham_orig[rownames(info), 'similarity']

# Confirmed that correlation highest between custom ex and observed AA subs between pairs.

# Crude validation that very simplistic approach yields same-ish answer, as a sanity check.
