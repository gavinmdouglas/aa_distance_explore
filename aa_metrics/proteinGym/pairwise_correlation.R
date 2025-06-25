rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(circlize)

gene_vals <- list()

for (infile in list.files(path="/Users/gavin/working_exchange_matrix/",
                          pattern="*.tsv",
                          full.names=TRUE)) {
  intab <- read.table(infile, header=TRUE, sep = '\t', stringsAsFactors = FALSE)
  
  dms_id <- sub("_per_gene_standard_scores_log.tsv$", "", basename(infile))
  
  gene_vals[[dms_id]] <- intab
}

exp_aas <- gene_vals$A4D664_9INFA_Soh_2019$aa_combo
for (i in 1:length(gene_vals)) {
  if (any(gene_vals[[i]]$aa_combo != exp_aas)) {
    stop(paste("Mismatch in aa_combo for", names(gene_vals)[i], "expected", exp_aas, "got", gene_vals[[i]]$aa_combo)) 
  }
}

sort(sapply(gene_vals, function(x) { sum(is.na(x$standardized_score))}))

cor_matrix <- data.frame(matrix(0, nrow=length(gene_vals), ncol=length(gene_vals)))
rownames(cor_matrix) <- names(gene_vals)
colnames(cor_matrix) <- names(gene_vals)
for (i in 1:length(gene_vals)) {
  for (j in 1:length(gene_vals)) {
    if (i != j) {
      cor_matrix[i, j] <- cor(gene_vals[[i]]$standardized_score,
                              gene_vals[[j]]$standardized_score,
                              method='spearman', use='pairwise.complete.obs')
    } else {
      cor_matrix[i, j] <- 1
    }
  }
}

Heatmap(as.matrix(cor_matrix),
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

# Outlier based on correlation (negative in most cases, which is clearly an issue!).
# Also remove one dataset from Weile 2017, so there is only one dataset for this study.

ids_to_exclude <- c("MK01_HUMAN_Brenan_2016", "SUMO1_HUMAN_Weile_2017")

# Keep only one Tsuboyama_2023 dataset, as otherwise it will swamp the signal.
# Kept the dataset with the fewest missing data-points: CBPA2_HUMAN_Tsuboyama_2023_1O6X
Tsuboyama_ids <- grep("Tsuboyama", names(gene_vals), value=TRUE)
Tsuboyama_ids_to_rm <- Tsuboyama_ids[-which(Tsuboyama_ids == "CBPA2_HUMAN_Tsuboyama_2023_1O6X")]
ids_to_exclude <- c(ids_to_exclude, Tsuboyama_ids_to_rm)
Tsuboyama_cor_matrix <- cor_matrix[Tsuboyama_ids, Tsuboyama_ids]

cor_matrix_filt <- cor_matrix[!rownames(cor_matrix) %in% ids_to_exclude, !colnames(cor_matrix) %in% ids_to_exclude]

Heatmap(as.matrix(cor_matrix_filt),
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

write.table(x=rownames(cor_matrix_filt),
            file="/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_filt_dms_ids.txt",
            sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
