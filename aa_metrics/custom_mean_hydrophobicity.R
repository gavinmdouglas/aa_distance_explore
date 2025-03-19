rm(list = ls(all.names = TRUE))

# Get single hydrophobicity measure per AA based on
# 1) Mathew2008 polar requirement
# 2) Metric #14 from Simm2016.

Simm2016 <- read.table('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/literature_data/Simm2016_TableS1_hydrophobicity.tsv',
                       header = TRUE, sep = '\t', stringsAsFactors = FALSE)

Mathew2008 <- read.table('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/literature_data/Mathew2008_Table1_polar_requirement.tsv',
                         header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 2)

grantham_polarity <- read.table('/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/literature_data/Grantham1974_polarity.tsv',
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 3)

Simm_AA_tri_order <- colnames(Simm2016[6:25])

polarity <- data.frame(aa=Mathew2008[Simm_AA_tri_order, 'AA'],
                       aa_tri=Simm_AA_tri_order,
                       original_grantham_polarity=grantham_polarity[Simm_AA_tri_order, 'Grantham1974_polarity'],
                       hydro=as.numeric(Simm2016[which(Simm2016$Index == 14), Simm_AA_tri_order]),
                       pr=Mathew2008[Simm_AA_tri_order, 'PR_calc_Mathew2008'])

# Multiplied by -1 first, because the the values are strongly negatively correlated.
polarity$hydroneg <- polarity$hydro * -1

# Min-max scale both columns.
polarity$hydroneg_minmax <- (polarity$hydroneg - min(polarity$hydroneg)) / (max(polarity$hydroneg) - min(polarity$hydroneg))
polarity$pr_minmax <- (polarity$pr - min(polarity$pr)) / (max(polarity$pr) - min(polarity$pr))

polarity$mean_minmax <- rowMeans(polarity[, c('hydroneg_minmax', 'pr_minmax')])

write.table(x = polarity,
            file = '/Users/gavin/Drive/ncsu/aa_selection/data/aa_physiochem_metrics/literature_data/polarity_mean_minmax.tsv',
            sep = '\t',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)
