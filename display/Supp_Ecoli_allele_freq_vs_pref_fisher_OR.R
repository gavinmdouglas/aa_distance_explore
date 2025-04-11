rm(list = ls(all.names = TRUE))

# Heatmaps to summarize odds ratios results over a wide range of cut-offs,
# to ensure that our results are robust.

library(ComplexHeatmap)
library(circlize)

all_output <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_seg_subs_extreme_prefs_by_freq_fisher_OR.tsv.gz',
                         header=TRUE, sep = '\t', stringsAsFactors = FALSE)

metrics_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

all_output$pref_clean <- metrics_map[all_output$pref, 'Clean']

prefs_to_rm <- c("Grantham (update, recalc.)", "Grantham (update, min-max)", "Grantham (recalc.)", "Grantham (min-max)",
                 "Miyata (recalc.)", "Miyata (min-max)", "Miyata (update, min-max)", "Miyata (update)", "Graur index (inverted)")
all_output <- all_output[which(! all_output$pref_clean %in% prefs_to_rm), ]

all_output$pref_clean[which(all_output$pref_clean == "Grantham (orig.)")] <- 'Grantham'
all_output$pref_clean[which(all_output$pref_clean == "Miyata (orig.)")] <- 'Miyata'

unique_prefs <- unique(all_output$pref_clean)
OR_tab <- data.frame(matrix(NA, nrow = length(unique_prefs), ncol = 12))
rownames(OR_tab) <- unique_prefs

P_tab <- OR_tab

maf_cutoffs = c(0.0001, 0.001, 0.01)
quantile_cutoffs = c(0.001, 0.01, 0.1, 0.25)
maf_grouping_name <- c()
quantile_grouping_name <- c()
col_num <- 1
for (maf_cutoff in maf_cutoffs) {
  for (cutoff in quantile_cutoffs) {
    maf_grouping_name <- c(maf_grouping_name, paste('MAF <', format(maf_cutoff, scientific = FALSE)))
    quantile_grouping_name <- c(quantile_grouping_name, paste('Pref. quantile <', format(cutoff, scientific = FALSE)))
    for (pref in unique_prefs) {
      OR_tab[pref, col_num] <- all_output$fisher_OR[which(all_output$pref_clean == pref & all_output$maf_cutoff == maf_cutoff & all_output$pref_cutoff == cutoff)]
      P_tab[pref, col_num] <- all_output$fisher_p[which(all_output$pref_clean == pref & all_output$maf_cutoff == maf_cutoff & all_output$pref_cutoff == cutoff)]
    }
    col_num <- col_num + 1
  }
}

pref_means <- sort(rowMeans(OR_tab, na.rm = TRUE), decreasing = TRUE)

OR_tab <- OR_tab[names(pref_means), ]
P_tab <- P_tab[names(pref_means), ]

OR_tab_num <- apply(OR_tab, 2, function(x) { sprintf('%.2f', x)})

OR_tab[P_tab >= 0.05] <- NA

ecoli_OR_heatmap <- Heatmap(matrix = as.matrix(OR_tab),
                              heatmap_legend_param = list(title = "Odds ratio"),
                              row_names_side = "left",
                              na_col = 'grey80',
                              column_labels = maf_grouping_name,
                              column_split = quantile_grouping_name,
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              row_title_rot = 0,
                              
                              column_names_rot = 45,
                              column_gap = unit(5, "mm"),
                              col = colorRamp2(c(floor(min(all_output$fisher_OR) * 10) / 10 , 1.0, ceiling(max(all_output$fisher_OR) * 10) / 10 ), c("slateblue1", "white", "red")),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(OR_tab_num[i, j], x, y, gp = gpar(fontsize = 10))
                              })

pdf("~/Drive/research/aa_distance/aa_distance_ms/display/Supp_ecoli_extreme_prefs_by_freq_fisher_OR.pdf",
    width = 10, height = 8)

draw(ecoli_OR_heatmap)

dev.off()

