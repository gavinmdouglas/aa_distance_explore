Code for manuscript entitled "DEX: an amino acid exchangeability measure for codon substitution modelling and selection inference" by Gavin Douglas and Louis-Marie Bobay.

Repository structure:
- `aa_metrics` - Code for processing various AA measures. Includes code for parsing custom measure from proteinGym database.
- `allele_freq_vs_predicted_effects` - Workflow for prepping segregating non-synonymous substitution allele frequencies across _E. coli_ and humans. See README in each respective folder for details.
- `compute_prefs` - Code for parsing a sequence to measure preferences (now just used for parsing VespaG preferences)
- `display` - Code for making figures
- `fasta_processing` - General code for processing FASTAs
- `functions.py` - Python script with a few functions used across multiple scripts
- `paml_workflow` - Code for running PAML workflow
- `text_results` - Quick scripts to reproducibly produce specific details/results to report in the main text.
