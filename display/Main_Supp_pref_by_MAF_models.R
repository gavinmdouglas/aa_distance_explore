rm(list = ls(all.names = TRUE))

# Plots of R-squared of how well mean allele frequency is explained by DMS-based preference scores, for both E. coli and human variants.
# Key question is whether buried vs. exposed exchangeabilities provide info beyond simply encoding a burial intercept,
# and whether more generally asymmetric exchangeabilities are informative.

library(ggplot2)

ecoli <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_pref_vs_AF_model_results.tsv.gz',
                    header=TRUE, sep="\t", stringsAsFactors = FALSE)

ecoli$species <- "Escherichia coli"

human <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_pref_vs_AF_model_results.tsv.gz',
                    header=TRUE, sep="\t", stringsAsFactors = FALSE)
human$species <- "Human"

intersecting_col <- intersect(colnames(ecoli), colnames(human))

tab <- rbind(ecoli[, intersecting_col], human[, intersecting_col])

# For relabelling metrics:
metric_labels <- list(
  RSA_only           = "RSA group only",
  mut_only           = "Mutation only",
  no_metric          = "RSA group and mutation only",
  dex_symmetric      = "DEX (symmetric)",
  dms_ex_symmetric   = "DMS-EX (symmetric)",
  ex_symmetric       = "EX (symmetric)",
  demask_symmetric   = "DeMaSk (symmetric)",
  grantham_symmetric = "Grantham (symmetric)",
  miyata_symmetric   = "Miyata (symmetric)",
  blosum62_symmetric = "BLOSUM62 (symmetric)",
  vtml200_symmetric  = "VTML200 (symmetric)",
  ex                 = "EX (asymmetric)",
  demask             = "DeMaSk (asymmetric)",
  epstein            = "Epstein (asymmetric)",
  vespag             = "VespaG (asymmetric)",
  thermoMPNN         = "ThermoMPNN (asymmetric)",
  rasp               = "RaSP (asymmetric)",
  dms.ex_weighted    = "DMS-EX weighted (asymmetric)",
  dms_ex_both        = "DMS-EX by RSA group (asymmetric)",
  dms.ex_ex          = "DMS-EX + EX (asymmetric)",
  dms.ex_demask      = "DMS-EX + DeMaSk (asymmetric)",
  ex_demask          = "EX + DeMaSk (asymmetric)",
  dms.ex_ex_demask   = "DMS-EX + EX + DeMaSk (asymmetric)"
)

tab$metric_clean <- sapply(tab$metric, function(x) metric_labels[[x]])

# First plot main figure, focused on subset of metrics needed to investigate specific question
# about importance of symmetry and burial status.

# This should be focused on mean allele frequency models.

focal_metrics <- c("RSA_only", "mut_only", "no_metric", "dex_symmetric", "dms_ex_symmetric", 
                   "dms.ex_weighted", "dms_ex_both")

tab_focal <- tab[tab$metric %in% focal_metrics, ]
tab_focal <- tab_focal[tab_focal$AF_var == "AF", ]

human_tab_focal <- tab_focal[tab_focal$species == "Human", ]
human_order <- human_tab_focal[order(human_tab_focal$adjusted_rsquared), "metric_clean"]

tab_focal$metric_clean <- factor(tab_focal$metric_clean, levels = human_order)

focal_af_models <- ggplot(data = tab_focal, aes(x = adjusted_rsquared, y=metric_clean)) +
                          geom_point(size = 3, colour="grey30") +
                          facet_wrap(~ species, scales = "free_x",
                                     labeller = as_labeller(c("Escherichia coli" = "italic('Escherichia coli')",
                                                              "Human" = "'Human'"),
                                                            label_parsed)) +
                          theme_bw() +
                          labs(x = expression(paste("Adjusted ", R^2)), y = "Model") +
                          theme(strip.background = element_blank()) +
                          scale_x_continuous(limits = c(-0.05, NA),
                                             expand = expansion(mult = c(0, 0.1)))

ggsave(plot = focal_af_models,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Main_focal_metrics_mean_AF_models.pdf',
       width = 7, height = 4, device=cairo_pdf)


# Make the analogous plot, but for proportion rare variants instead, which was an alternative AF-based metric I explored.

supp_focal <- tab[tab$metric %in% focal_metrics, ]
supp_focal <- supp_focal[supp_focal$AF_var == "prop_rare", ]

human_supp_focal <- supp_focal[supp_focal$species == "Human", ]
human_supp_order <- human_supp_focal[order(human_supp_focal$adjusted_rsquared), "metric_clean"]

supp_focal$metric_clean <- factor(supp_focal$metric_clean, levels = human_supp_order)

supp_focal_prop_rare_models <- ggplot(data = supp_focal, aes(x = adjusted_rsquared, y=metric_clean)) +
  geom_point(size = 3, colour="grey30") +
  facet_wrap(~ species, scales = "free_x",
             labeller = as_labeller(c("Escherichia coli" = "italic('Escherichia coli')",
                                      "Human" = "'Human'"),
                                    label_parsed)) +
  theme_bw() +
  labs(x = expression(paste("Adjusted ", R^2)), y = "Model") +
  theme(strip.background = element_blank()) +
  scale_x_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1)))

ggsave(plot = supp_focal_prop_rare_models,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Supp_focal_metrics_prop_rare_models.pdf',
       width = 7, height = 4, device=cairo_pdf)


# Now plot full plots of all metrics, just for broader context which can be provided in supplement.
supp_af <- tab[tab$AF_var == "AF", ]

supp_af_human <- supp_af[supp_af$species == "Human", ]
supp_af_human_order <- supp_af_human[order(supp_af_human$adjusted_rsquared), "metric_clean"]

supp_af$metric_clean <- factor(supp_af$metric_clean, levels = supp_af_human_order)

supp_af_models <- ggplot(data = supp_af, aes(x = adjusted_rsquared, y=metric_clean)) +
  geom_point(size = 3, colour="grey30") +
  facet_wrap(~ species, scales = "free_x",
             labeller = as_labeller(c("Escherichia coli" = "italic('Escherichia coli')",
                                      "Human" = "'Human'"),
                                    label_parsed)) +
  theme_bw() +
  labs(x = expression(paste("Adjusted ", R^2)), y = "Model") +
  theme(strip.background = element_blank()) +
  scale_x_continuous(limits = c(-0.05, NA),
                     expand = expansion(mult = c(0, 0.1)))

ggsave(plot = supp_af_models,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Supp_broader_mean_AF_models.pdf',
       width = 7, height = 6, device=cairo_pdf)


supp_prop_rare <- tab[tab$AF_var == "prop_rare", ]

supp_prop_rare_human <- supp_prop_rare[supp_prop_rare$species == "Human", ]
supp_prop_rare_human_order <- supp_prop_rare_human[order(supp_prop_rare_human$adjusted_rsquared), "metric_clean"]

supp_prop_rare$metric_clean <- factor(supp_prop_rare$metric_clean, levels = supp_prop_rare_human_order)

supp_prop_rare_models <- ggplot(data = supp_prop_rare, aes(x = adjusted_rsquared, y=metric_clean)) +
  geom_point(size = 3, colour="grey30") +
  facet_wrap(~ species, scales = "free_x",
             labeller = as_labeller(c("Escherichia coli" = "italic('Escherichia coli')",
                                      "Human" = "'Human'"),
                                    label_parsed)) +
  theme_bw() +
  labs(x = expression(paste("Adjusted ", R^2)), y = "Model") +
  theme(strip.background = element_blank()) +
  scale_x_continuous(limits = c(-0.05, NA),
                     expand = expansion(mult = c(0, 0.1)))

ggsave(plot = supp_prop_rare_models,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Supp_broader_mean_prop_rare_models.pdf',
       width = 7, height = 6, device=cairo_pdf)
