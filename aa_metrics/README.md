### Parsing amino acid metrics

- aa_pairwise_diff_func.R - General function for computing amino acid (AA) pairwise distances based on a set of input variables and a specified normalization (can be used to compute different versions of Grantham distances for instance).

- aa_radical_vs_cons_define.R - R objects defined that contain definitions of radical and conservative substitutions of different types.

- codon_mean_exchangeabilitity.py - Script to get metric similar to Graur's stability (but the inverse) for each codon based on the 'similarity' metrics. It is the average similarity of all amino acids that arise through all possible single mutations that cause substitutions for a given codon (written to: aa_metrics/codon_mean_exchangeabilities.tsv.gz).

- compute_custom_atchley_based_distance.R - Commands to compute AA distances based on Atchley factors.

- compute_grantham_miyata_variants.R - Commands to compute different types of Grantham and Miyata distances.

- custom_mean_hydrophobicity.R - Compute updated polarity measure to use in Grantham and Miyata calculations.

- explore_aa_volumes.R - Quick comparison of AA volumes reported in different papers.

- parse_PARD_metrics.py - Python script to parse several distance/similarity metrics from the PARD package.

- prep_all_metrics.R - Commands to process all the AA metrics of interest, to output them to a consistent (long) format, which is a simple similarity measure ranging from > 0 to <= 1.

- sanity_check_prepped_metrics.R - Some quick checks on the prepped versions of the metric files.