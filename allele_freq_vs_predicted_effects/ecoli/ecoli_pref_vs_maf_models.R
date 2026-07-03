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

# Code transition freq as binary, as it's mainly all or none for a given AA combo anyway.
subtab_filt_asymmetric$mainly_transitions <- "No"
subtab_filt_asymmetric$mainly_transitions[subtab_filt_asymmetric$prop_transition > 0.5] <- "Yes"

subtab_filt_asymmetric$mainly_transitions <- factor(subtab_filt_asymmetric$mainly_transitions, c('No', 'Yes'))

subtab_filt_asymmetric$RSA_group <- factor(subtab_filt_asymmetric$RSA_group, c('buried', 'exposed'))

subtab_filt_asymmetric$prop_rare_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$prop_rare)$x.t
subtab_filt_asymmetric$AF_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$AF)$x.t

subtab_filt_asymmetric$prop_rare_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$prop_rare)$x.t
subtab_filt_asymmetric$AF_orderedNorm <- bestNormalize::orderNorm(subtab_filt_asymmetric$AF)$x.t

metrics <- c("RSA_only", "mut_only", "no_metric", "dms.ex_weighted",
             "dex_symmetric", "dms_ex_symmetric",
             "dms_ex_both", "vespag", "thermoMPNN", "ex",
             "demask", "epstein", "dms.ex_ex",
             "dms.ex_demask", "ex_demask", "dms.ex_ex_demask",
             "ex_symmetric", "demask_symmetric",
             "grantham_symmetric", "miyata_symmetric",
             "blosum62_symmetric", "vtml200_symmetric")

AF_vars <- c("AF", "prop_rare")

results_raw <- list()
for (AF_var in AF_vars) {
  raw <- lapply(metrics, function(metric) {
    
    if (metric == "RSA_only") {
      form <- as.formula(paste0(AF_var, "_orderedNorm ~ RSA_group"))
    } else if (metric == "mut_only") {
      form <- as.formula(paste0(AF_var, "_orderedNorm ~ mainly_transitions"))
    } else if (metric == "no_metric") {
      form <- as.formula(paste0(AF_var, "_orderedNorm ~ RSA_group + mainly_transitions"))
    } else {
      subtab_filt_asymmetric[, paste0(metric, "_orderedNorm")] <- bestNormalize::orderNorm(subtab_filt_asymmetric[, metric])$x.t
      metric_name <- paste0(metric, "_orderedNorm")
      form <- as.formula(paste0(AF_var, "_orderedNorm ~ RSA_group + mainly_transitions + ", metric_name))
    }
    model <- lm(form, data = subtab_filt_asymmetric)
    cf <- coef(model)
    
    if (metric == "RSA_only") {
      metric_coef <- "NA"
      RSA_group_coef <- cf["RSA_groupexposed"]
      mainly_transitions_coef <- NA
    } else if (metric == "mut_only") {
      metric_coef <- "NA"
      RSA_group_coef <- "NA"
      mainly_transitions_coef <- cf["mainly_transitions"]
    } else if (metric == "no_metric") {
      metric_coef <- "NA"
      RSA_group_coef <- cf["RSA_groupexposed"]
      mainly_transitions_coef <- cf["mainly_transitions"]
    } else {
      metric_coef <- cf[metric_name]
      RSA_group_coef <- cf["RSA_groupexposed"]
      mainly_transitions_coef <- cf["mainly_transitions-CpG transition"]
    }
    
    data.frame(AF_var = AF_var,
               metric = metric,
               adjusted_rsquared = summary(model)$adj.r.squared,
               AIC = AIC(model),
               BIC = BIC(model),
               logLik = as.numeric(logLik(model)),
               intercept = cf["(Intercept)"],
               RSA_group_exposed = RSA_group_coef,
               mainly_transitions_coef = mainly_transitions_coef,
               metric_coef = metric_coef)
  })
  results_raw[[AF_var]] <- do.call(rbind, raw)
}

results <- do.call(rbind, results_raw)

write.table(x=results,
            file=gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/ecoli_pref_vs_AF_model_results.tsv.gz'),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
