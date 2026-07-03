rm(list = ls(all.names = TRUE))

# Compute correlations mean preferences vs. minor allele frequencies for each amino acid pair.
metrics_map <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/metrics_raw_to_clean.tsv.gz',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

aa_sim_metrics <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/combined_symmetric_prepped_similarity_metrics.tsv.gz",
                             header=TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

# Compute correlations mean preferences vs. minor allele frequencies for each amino acid pair.
asymmetric_metrics <- read.table("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/asymmetric_similarity_consistent.tsv.gz",
                                 header=TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

rownames(asymmetric_metrics) <- apply(asymmetric_metrics[, c('ref_aa', 'mut_aa')], 1, function(x) {
  paste(x, collapse = '_')
})
asymmetric_metrics <- asymmetric_metrics[, -which(colnames(asymmetric_metrics) %in% c('ref_aa', 'mut_aa'))]

# Add in some unfolded symmetric measures.
# DEX, demask, EX, Grantham, Miyata, BLOSUM62, and VTML
asymmetric_metrics$dex_symmetric <- aa_sim_metrics[rownames(asymmetric_metrics), 'DISTATIS_Custom_EX']
asymmetric_metrics$dms_ex_symmetric <- aa_sim_metrics[rownames(asymmetric_metrics), 'proteinGym_rsa_median_custom']
asymmetric_metrics$ex_symmetric <- aa_sim_metrics[rownames(asymmetric_metrics), 'ex']
asymmetric_metrics$demask_symmetric <- aa_sim_metrics[rownames(asymmetric_metrics), 'demask']
asymmetric_metrics$grantham_symmetric <- aa_sim_metrics[rownames(asymmetric_metrics), 'grantham_orig']
asymmetric_metrics$miyata_symmetric <- aa_sim_metrics[rownames(asymmetric_metrics), 'miyata_orig']
asymmetric_metrics$blosum62_symmetric <- aa_sim_metrics[rownames(asymmetric_metrics), 'blosum62.prob']
asymmetric_metrics$vtml200_symmetric <- aa_sim_metrics[rownames(asymmetric_metrics), 'vtml200.prob']

sub_tab_filt <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/parsed_allelome_variants_w_prefs_dssp_w_anc_prep.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)

sub_tab_filt$log10_AF <- log10(sub_tab_filt$AF)

mean_and_prop_rare_by_aa_combo <- function(indf,
                             focal_col= c("AF", "log10_AF", "vespag", "thermoMPNN"),
                             min_freq=10) {
  
  aa_combo_tallies <- table(indf$aa_combo)
  
  combo_to_keep <- names(aa_combo_tallies)[which(aa_combo_tallies >= min_freq)]
  
  indf <- indf[which(indf$aa_combo %in% combo_to_keep), ]
  
  mean_df <- aggregate(indf[, focal_col], by = list(aa_combo = indf$aa_combo), FUN = mean, na.rm=TRUE)
  
  rownames(mean_df) <- mean_df$aa_combo

  # Get proportion below frequency cut-off for this AA combo.
  # Also get prop transition.
  mean_df$prop_rare <- NA
  mean_df$prop_transition <- NA
  for (aa_comb in rownames(mean_df)) {
    tab_subset <- indf[which(indf$aa_combo == aa_comb), ]
    mean_df[aa_comb, "prop_rare"] <- sum(tab_subset$AF < 0.001) / nrow(tab_subset)
    mean_df[aa_comb, "prop_transition"] <- sum(tab_subset$Ts_or_Tv == 'Ts') / nrow(tab_subset)
  }

  for (aa_metric in colnames(asymmetric_metrics)) {
    mean_df[, aa_metric] <- asymmetric_metrics[rownames(mean_df), aa_metric]
  }
  
  mean_df
}

subtab_filt_asymmetric_buried <- mean_and_prop_rare_by_aa_combo(indf = sub_tab_filt[sub_tab_filt$RSA_group == "buried", ], min_freq = 10)
subtab_filt_asymmetric_exposed <- mean_and_prop_rare_by_aa_combo(indf = sub_tab_filt[sub_tab_filt$RSA_group == "exposed", ], min_freq = 10)

subtab_filt_asymmetric <- rbind(subtab_filt_asymmetric_buried, subtab_filt_asymmetric_exposed)
subtab_filt_asymmetric$RSA_group <- c(rep("buried", nrow(subtab_filt_asymmetric_buried)), rep("exposed", nrow(subtab_filt_asymmetric_exposed)))

subtab_filt_asymmetric$dms_ex_both_naive <- c(subtab_filt_asymmetric_buried$dms.ex_buried, subtab_filt_asymmetric_exposed$dms.ex_exposed)

rownames(subtab_filt_asymmetric) <- paste(subtab_filt_asymmetric$aa_combo, subtab_filt_asymmetric$RSA_group, sep = "_")

combined_buried_exposed <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/proteinGym_asymmetric_concat_rsa_group_similarity.tsv.gz',
                                      header=TRUE, sep="\t", stringsAsFactors = FALSE)

rownames(combined_buried_exposed) <- paste(combined_buried_exposed$ref_aa, combined_buried_exposed$mut_aa, combined_buried_exposed$rsa_group, sep = "_")

subtab_filt_asymmetric$dms_ex_both <- combined_buried_exposed[rownames(subtab_filt_asymmetric), 'combined_buried_exposed_similiarity']




focal_cols <- c("prop_rare", "AF")
raw_spearman_out <- list()

indf <- subtab_filt_asymmetric

for (focal_col in focal_cols) {
  
  tested_names <- character()
  spearman_rho <- numeric()
  spearman_p <- numeric()
  for (coln in colnames(indf)) {
    if (coln == "AF") { next }
    if (coln == "aa_combo") { next }
    if (coln == "prop_rare") { next }
    if (coln == "RSA_group") { next }
    if (coln == "log10_AF") { next }
    if (coln == "prop_transition") { next }
    test_res <- cor.test(indf[, coln], indf[, focal_col], method = 'spearman')
    tested_names <- c(tested_names, coln)
    spearman_rho <- c(spearman_rho, test_res$estimate)
    spearman_p <- c(spearman_p, test_res$p.value)
  }
  
  spearman_tab <- data.frame(focal_col = focal_col,
                             Metric = tested_names,
                             Spearman_rho = spearman_rho,
                             Spearman_p = spearman_p)
  
  raw_spearman_out[[focal_col]] <- spearman_tab[order(spearman_tab$Spearman_rho), ]
  
}

spearman_out <- do.call(rbind, raw_spearman_out)






# Code transition freq as binary, as it's mainly all or none for a given AA combo anyway.
subtab_filt_asymmetric$mainly_transitions <- "No"
subtab_filt_asymmetric$mainly_transitions[subtab_filt_asymmetric$prop_transition > 0.5] <- "Yes"

subtab_filt_asymmetric$mainly_transitions <- factor(subtab_filt_asymmetric$mainly_transitions, c('No', 'Yes'))

subtab_filt_asymmetric$RSA_group <- factor(subtab_filt_asymmetric$RSA_group, c('buried', 'exposed'))

subtab_filt_asymmetric$dms_ex_both_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$dms_ex_both)$x.t
subtab_filt_asymmetric$thermoMPNN_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$thermoMPNN)$x.t
subtab_filt_asymmetric$vespag_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$vespag)$x.t
subtab_filt_asymmetric$dms_ex_symmetric_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$dms_ex_symmetric)$x.t
subtab_filt_asymmetric$dex_symmetric_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$dex_symmetric)$x.t

subtab_filt_asymmetric$prop_rare_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$prop_rare)$x.t
subtab_filt_asymmetric$log10_AF_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$log10_AF)$x.t

metrics <- c(
  "vespag", "thermoMPNN", "dms.ex_weighted",
  "dex_symmetric", "dms_ex_symmetric",
  "dms_ex_both")

results <- lapply(metrics, function(metric) {
  
  subtab_filt_asymmetric[, paste0(metric, "_orderedNorm")] <- bestNormalize::orderNorm(subtab_filt_asymmetric[, metric])$x.t
  metric_name <- paste0(metric, "_orderedNorm")
  
  form <- as.formula(paste0("log10_AF_orderedNorm ~ RSA_group + mainly_transitions + ", metric_name))
  model <- lm(form, data = subtab_filt_asymmetric)
  cf <- coef(model)
  
  data.frame(
    metric = metric,
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = as.numeric(logLik(model)),
    intercept = cf["(Intercept)"],
    RSA_group_exposed = cf["RSA_groupexposed"],
    metric_coef = cf[metric_name]
  )
})

results <- do.call(rbind, results)

tmp1 <- results[order(results$AIC), ]
tmp1$focal <- 'log10_AF'

tmp2 <- results[order(results$AIC), ]
tmp2$focal <- 'prop_rare'

results_combined <- rbind(tmp1, tmp2)






library(ggplot2)
library(dplyr)

# compute delta BIC within each outcome (focal)
plot_df <- results_combined %>%
  group_by(focal) %>%
  mutate(dBIC = BIC - min(BIC)) %>%
  ungroup()

# assign measure type
plot_df <- plot_df %>%
  mutate(type = case_when(
    metric %in% c("vespag", "thermoMPNN", "rasp") ~ "reference (site-level)",
    metric == "dms_ex_both" ~ "burial-conditioned",
    metric %in% c("dms.ex_weighted") ~ "asymmetric",
    TRUE ~ "symmetric"
  ),
  type = factor(type, levels = c("symmetric", "asymmetric",
                                 "burial-conditioned", "reference (site-level)")))

# order measures by dBIC in the informative outcome (log10_AF), apply to both facets
metric_order <- plot_df %>%
  filter(focal == "log10_AF") %>%
  arrange(dBIC) %>%
  pull(metric)
plot_df$metric <- factor(plot_df$metric, levels = rev(metric_order))

# nicer facet labels
plot_df$focal <- factor(plot_df$focal, levels = c("log10_AF", "prop_rare"),
                        labels = c("log10 allele frequency", "proportion rare"))

ggplot(plot_df, aes(x = dBIC, y = metric, fill = type)) +
  geom_col(width = 0.7) +
  facet_wrap(~ focal) +
  scale_fill_manual(values = c(
    "symmetric"               = "#2a78d6",
    "asymmetric"              = "#eb6834",
    "burial-conditioned"      = "#4a3aa7",
    "reference (site-level)"  = "#898781")) +
  labs(x = expression(Delta*"BIC from best model"), y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    strip.text = element_text(face = "plain", hjust = 0)
  )







# buried stratum
buried_tab <- subtab_filt_asymmetric[subtab_filt_asymmetric$RSA_group == "buried", ]
exposed_tab <- subtab_filt_asymmetric[subtab_filt_asymmetric$RSA_group == "exposed", ]

cor_in_stratum <- function(df, metric, outcome = "log10_AF") {
  ok <- complete.cases(df[[metric]], df[[outcome]])
  ct <- cor.test(df[[metric]][ok], df[[outcome]][ok], method = "spearman")
  data.frame(metric = metric, n = sum(ok), rho = ct$estimate, p = ct$p.value)
}

# pick the environment-appropriate dms.ex variant per stratum
buried_metrics  <- c("dex_symmetric", "dms_ex_symmetric", "ex_symmetric", "demask_symmetric",
                     "grantham_symmetric", "miyata_symmetric", "blosum62_symmetric", "vtml200_symmetric",
                     "dms.ex_buried", "ex", "demask", "vespag", "thermoMPNN", "epstein", "dms.ex_weighted")
exposed_metrics <- sub("dms.ex_buried", "dms.ex_exposed", buried_metrics)

buried_res  <- do.call(rbind, lapply(buried_metrics,  function(m) cor_in_stratum(buried_tab, m)))
exposed_res <- do.call(rbind, lapply(exposed_metrics, function(m) cor_in_stratum(exposed_tab, m)))

buried_res[order(-abs(buried_res$rho)), ]
exposed_res[order(-abs(exposed_res$rho)), ]









results <- lapply(metrics, function(metric) {
  
  subtab_filt_asymmetric[, paste0(metric, "_orderedNorm")] <- bestNormalize::orderNorm(subtab_filt_asymmetric[, metric])$x.t
  
  metric_name <- paste0(metric, "_orderedNorm")
  
  form <- as.formula(
    paste0(
      "prop_rare_orderedNorm ~ RSA_group + ", metric_name,
      " + mainly_transitions"
    )
  )
  
  model <- lm(form, data = subtab_filt_asymmetric)
  
  model_summary <- summary(model)
  
  data.frame(
    metric = metric,
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = as.numeric(logLik(model)),
    Rsquared = model_summary$r.squared,
    intercept = model_summary$coefficients["(Intercept)", "Estimate"],
    RSA_group_exposed = model_summary$coefficients["RSA_groupexposed", "Estimate"],
    metric = model_summary$coefficients[metric_name, "Estimate"],
    mainly_transitions = model_summary$coefficients["mainly_transitionsYes", "Estimate"])
})

results <- do.call(rbind, results)

results[order(results$AIC), ]



# write.table(x=spearman_out,
#             file=gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_spearman_pref_vs_maf_asymmetric_w_singletons.tsv.gz'),
#             sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(x=spearman_out,
            file=gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_spearman_pref_vs_maf_asymmetric_no_singletons.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)