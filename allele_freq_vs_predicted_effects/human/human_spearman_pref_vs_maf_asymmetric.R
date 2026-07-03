rm(list = ls(all.names = TRUE))

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

sub_tab_filt <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_snvs_w_prefs_dssp_anc_aa_prep.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)

sub_tab_filt$log10_AF <- log10(sub_tab_filt$AF)

# Note the sum-scaled variable are simply internal checks for how much the signal disappears if per-site each "preference" is scaled to sum to 1.
# Meaning differences across sites/regions are ignored.
prefs <- c("vespa", "vespag_sum_scaled", "thermoMPNN", "rasp")

mean_and_prop_rare_by_aa_combo <- function(indf,
                                           focal_col= c("AF", "log10_AF", "vespag", "vespag_sum_scaled", "thermoMPNN", "rasp"),
                                           min_freq=50) {

  aa_combo_tallies <- table(indf$aa_combo)
  
  combo_to_keep <- names(aa_combo_tallies)[which(aa_combo_tallies >= min_freq)]
  
  indf <- indf[which(indf$aa_combo %in% combo_to_keep), ]
  
  mean_df <- aggregate(indf[, focal_col], by = list(aa_combo = indf$aa_combo), FUN = mean, na.rm=TRUE)
  
  rownames(mean_df) <- mean_df$aa_combo
  mean_df$prop_rare <- NA
  for (aa_comb in rownames(mean_df)) {
   
    tab_subset <- indf[which(indf$aa_combo == aa_comb), ]
    mean_df[aa_comb, "prop_rare"] <- sum(tab_subset$AF < 1e-06) / nrow(tab_subset)
    mean_df[aa_comb, "prop_transition"] <- sum(tab_subset$Ts_or_Tv == 'transition') / nrow(tab_subset)
    mean_df[aa_comb, "prop_CpG"] <- sum(tab_subset$CpG_or_not == 'CpG') / nrow(tab_subset)
  }
  
  for (aa_metric in colnames(asymmetric_metrics)) {
    mean_df[, aa_metric] <- asymmetric_metrics[rownames(mean_df), aa_metric]
  }
  
  mean_df
}

subtab_filt_asymmetric <- mean_and_prop_rare_by_aa_combo(indf = sub_tab_filt)
subtab_filt_asymmetric_buried <- mean_and_prop_rare_by_aa_combo(indf = sub_tab_filt[sub_tab_filt$RSA_group == "buried", ])
subtab_filt_asymmetric_exposed <- mean_and_prop_rare_by_aa_combo(indf = sub_tab_filt[sub_tab_filt$RSA_group == "exposed", ])

all_asymmetric_aa_combos <- table(c(subtab_filt_asymmetric$aa_combo, subtab_filt_asymmetric_buried$aa_combo, subtab_filt_asymmetric_exposed$aa_combo))
ubiq_asymmetric_aa_combos <- names(all_asymmetric_aa_combos)[which(all_asymmetric_aa_combos == 3)]
subtab_filt_asymmetric <- subtab_filt_asymmetric[which(subtab_filt_asymmetric$aa_combo %in% ubiq_asymmetric_aa_combos), ]
subtab_filt_asymmetric_buried <- subtab_filt_asymmetric_buried[which(subtab_filt_asymmetric_buried$aa_combo %in% ubiq_asymmetric_aa_combos), ]
subtab_filt_asymmetric_exposed <- subtab_filt_asymmetric_exposed[which(subtab_filt_asymmetric_exposed$aa_combo %in% ubiq_asymmetric_aa_combos), ]

subtab_filt_asymmetric <- rbind(subtab_filt_asymmetric_buried, subtab_filt_asymmetric_exposed)

subtab_filt_asymmetric$RSA_group <- c(rep("buried", nrow(subtab_filt_asymmetric_buried)), rep("exposed", nrow(subtab_filt_asymmetric_exposed)))

rownames(subtab_filt_asymmetric) <- paste(subtab_filt_asymmetric$aa_combo, subtab_filt_asymmetric$RSA_group, sep = "_")

combined_buried_exposed <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/proteinGym_asymmetric_concat_rsa_group_similarity.tsv.gz',
                                      header=TRUE, sep="\t", stringsAsFactors = FALSE)

rownames(combined_buried_exposed) <- paste(combined_buried_exposed$ref_aa, combined_buried_exposed$mut_aa, combined_buried_exposed$rsa_group, sep = "_")

subtab_filt_asymmetric$dms_ex_both <- combined_buried_exposed[rownames(subtab_filt_asymmetric), 'combined_buried_exposed_similiarity']

# Code transition freq as binary, as it's mainly all or none for a given AA combo anyway.
subtab_filt_asymmetric$mainly_mut_type <- "transversion"
subtab_filt_asymmetric$mainly_mut_type[subtab_filt_asymmetric$prop_transition > 0.5 & subtab_filt_asymmetric$prop_CpG <= 0.5] <- "Non-CpG transition"
subtab_filt_asymmetric$mainly_mut_type[subtab_filt_asymmetric$prop_transition > 0.5 & subtab_filt_asymmetric$prop_CpG > 0.5] <- "CpG transition"

subtab_filt_asymmetric$mainly_mut_type <- factor(subtab_filt_asymmetric$mainly_mut_type, c('transversion', 'Non-CpG transition', 'CpG transition'))

subtab_filt_asymmetric$RSA_group <- factor(subtab_filt_asymmetric$RSA_group, c('buried', 'exposed'))

subtab_filt_asymmetric$prop_rare_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$prop_rare)$x.t
subtab_filt_asymmetric$AF_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$AF)$x.t

metrics <- c("dms.ex_weighted",
  "dex_symmetric", "dms_ex_symmetric",
  "dms_ex_both")

results <- lapply(metrics, function(metric) {
  subtab_filt_asymmetric[, paste0(metric, "_orderedNorm")] <- bestNormalize::orderNorm(subtab_filt_asymmetric[, metric])$x.t
  metric_name <- paste0(metric, "_orderedNorm")
  
  form <- as.formula(paste0("AF_orderedNorm ~ RSA_group + mainly_mut_type + ", metric_name))
  model <- lm(form, data = subtab_filt_asymmetric)
  cf <- coef(model)
  
  data.frame(metric = metric,
             AIC = AIC(model),
             BIC = BIC(model),
             logLik = as.numeric(logLik(model)),
             intercept = cf["(Intercept)"],
             RSA_group_exposed = cf["RSA_groupexposed"],
             metric_coef = cf[metric_name])
})

results <- do.call(rbind, results)

results[order(results$AIC), ]

tmp1 <- results[order(results$AIC), ]
tmp1$focal <- 'log10_AF'

tmp2 <- results[order(results$AIC), ]
tmp2$focal <- 'prop_rare'

results_combined <- rbind(tmp1, tmp2)

 
library(ggplot2)
library(dplyr)

plot_df <- results_combined %>%
  group_by(focal) %>%
  mutate(dBIC = BIC - min(BIC)) %>%
  ungroup() %>%
  mutate(type = case_when(
    metric %in% c("vespag", "thermoMPNN") ~ "reference (site-level)",
    metric == "dms_ex_both" ~ "burial-conditioned",
    metric == "dms.ex_weighted" ~ "asymmetric",
    TRUE ~ "symmetric"
  ),
  type = factor(type, levels = c("symmetric", "asymmetric",
                                 "burial-conditioned", "reference (site-level)")))

metric_order <- plot_df %>%
  filter(focal == "log10_AF") %>%
  arrange(dBIC) %>%
  pull(metric)
plot_df$metric <- factor(plot_df$metric, levels = rev(metric_order))

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











groupings <- c("all", "buried", "exposed")
focal_cols <- c("prop_rare", "AF")
raw_spearman_out <- list()

for (grouping in groupings) {
  if (grouping == "all") {
    indf <- subtab_filt_asymmetric
  } else if (grouping == "buried") {
    indf <- subtab_filt_asymmetric_buried
  } else if (grouping == "exposed") {
    indf <- subtab_filt_asymmetric_exposed
  }

  for (focal_col in focal_cols) {
   
    tested_names <- character()
    spearman_rho <- numeric()
    spearman_p <- numeric()
    for (coln in colnames(indf)) {
      if (coln == "AF") { next }
      if (coln == "aa_combo") { next }
      if (coln == "prop_rare") { next }
      if (coln == "mainly_mut_type") { next }
      if (coln == "RSA_group") { next }
      test_res <- cor.test(indf[, coln], indf[, focal_col], method = 'spearman')
      tested_names <- c(tested_names, coln)
      spearman_rho <- c(spearman_rho, test_res$estimate)
      spearman_p <- c(spearman_p, test_res$p.value)
    }
    
    spearman_tab <- data.frame(focal_col = focal_col,
                               grouping = grouping,
                               Metric = tested_names,
                               Spearman_rho = spearman_rho,
                               Spearman_p = spearman_p)
    
    raw_spearman_out[[paste(grouping, focal_col)]] <- spearman_tab[order(spearman_tab$Spearman_rho), ]

  }
  
}

spearman_out <- do.call(rbind, raw_spearman_out)

write.table(x=spearman_out,
            file=gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/human_variants/human_spearman_pref_vs_maf_asymmetric.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
