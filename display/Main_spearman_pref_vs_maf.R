rm(list = ls(all.names = TRUE))

library(ggplot2)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)

ecoli_spearman <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_spearman_pref_vs_maf.tsv.gz',
                             header = TRUE, sep = '\t', stringsAsFactors = FALSE)

human_spearman <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_spearman_pref_vs_maf.tsv.gz',
                             header = TRUE, sep = '\t', stringsAsFactors = FALSE)

ecoli_spearman$Species <- 'E. coli'
human_spearman$Species <- 'Human'

combined <- rbind(ecoli_spearman, human_spearman)

extra_to_rm <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/nonfocal_metrics.txt',
                           header = FALSE, stringsAsFactors = FALSE)$V1

extra_to_rm <- c("proteinGym_rsa_mean_custom", extra_to_rm)

combined <- combined[! combined$Metric %in% extra_to_rm, ]
combined$Metric <- sub('-', '.', combined$Metric)
combined$Metric <- sub('_prob', '', combined$Metric)
combined$Metric <- sub('.prob', '', combined$Metric)
combined$Metric <- sub('_maxscaled', '', combined$Metric)

metrics_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

combined$Metric <- metrics_map[combined$Metric, 'Clean']
combined$Metric[which(combined$Metric == "Grantham (orig.)")] <- 'Grantham'
combined$Metric[which(combined$Metric == "Miyata (orig.)")] <- 'Miyata'
combined$Metric[which(combined$Metric == "This study (median)")] <- 'This study'

combined$Spearman_BH <- p.adjust(combined$Spearman_p, method = 'BH')
combined$significant <- combined$Spearman_BH < 0.05

human_tab <- combined[which(combined$Species == 'Human'), c('Metric', 'Spearman_rho', 'significant')]
ecoli_tab <- combined[which(combined$Species == 'E. coli'), c('Metric', 'Spearman_rho', 'significant')]

rownames(human_tab) <- human_tab$Metric
rownames(ecoli_tab) <- ecoli_tab$Metric

all_metrics <- unique(combined$Metric)

spearman_wide <- data.frame(ecoli_spearman = ecoli_tab[all_metrics, 'Spearman_rho'],
                            human_spearman = human_tab[all_metrics, 'Spearman_rho'])
rownames(spearman_wide) <- all_metrics

sig_wide <- data.frame(ecoli_sig = ecoli_tab[all_metrics, 'significant'],
                       human_sig = human_tab[all_metrics, 'significant'])
rownames(sig_wide) <- all_metrics


scale1 <- scale(spearman_wide$ecoli_spearman)
scale2 <- scale(spearman_wide$human_spearman)
mean_spearman_scale <- numeric()

for (i in 1:nrow(spearman_wide)) {
 mean_spearman_scale <- c(mean_spearman_scale,
                          mean(c(scale1[i], scale2[i]), na.rm = TRUE)) 
}

ordered_i <- order(mean_spearman_scale, decreasing = TRUE)

spearman_wide <- spearman_wide[ordered_i, ]
sig_wide <- sig_wide[ordered_i, ]

# Set any non-sig spearman to NA:
spearman_wide[! sig_wide] <- NA

spearman_wide_num <- apply(spearman_wide, 2, function(x) { sprintf('%.2f', x)})
spearman_wide_num[is.na(spearman_wide)] <- ''

spearman_wide_num[which(rownames(spearman_wide) == "RaSP"), 1] <- 'No data'

spearman_wide <- as.matrix(spearman_wide)


clean_to_group <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                             sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 2)

row_group <- clean_to_group[rownames(spearman_wide), 'Grouping']

group_colours <- c("Physicochemical/structural" = "black",
                   "Experimental"              = "#0072B2",
                   "Substitution-based"        = "#009E73",
                   "Radical vs. conservative"  = "#D55E00",
                   "Deep learning"             = "purple")
label_colours <- group_colours[row_group]

spearman_heatmap <- Heatmap(matrix = spearman_wide,
                            row_names_side = "left",
                            na_col = 'grey80',
                            cluster_rows = FALSE,
                            cluster_columns = FALSE,
                            row_title_rot = 0,
                            name = "Spearman's \u03c1",
                            column_labels = c(expression(italic("E. coli")), "Human"),
                            row_names_gp = gpar(col = label_colours),
                            column_names_rot = 45,
                            column_gap = unit(5, "mm"),
                            col = colorRamp2(c(0, 1), c("white", "dark red")),
                            cell_fun = function(j, i, x, y, width, height, fill) {
                              grid.text(spearman_wide_num[i, j], x, y, gp = gpar(fontsize = 10))
                            })

group_legend <- Legend(
  labels = names(group_colours),
  title = "Measure grouping",
  type = "points",
  pch = NA,
  background = "white",
  border = NA,
  grid_width = unit(0, "mm"),
  labels_gp = gpar(col = group_colours)
)

spearman_heatmap <- grid.grabExpr(draw(spearman_heatmap,
                                   padding = unit(c(2, 5, 2, 5), "mm"),
                                   merge_legend = TRUE,
                                   annotation_legend_list = list(group_legend),
                                   heatmap_legend_side = "right",
                                   annotation_legend_side = "right"))

ggsave(filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Main_allele_freq_spearman.pdf',
       plot = spearman_heatmap,
       device = cairo_pdf,
       width = 5,
       height = 7,
       units = 'in',
       dpi = 600)
