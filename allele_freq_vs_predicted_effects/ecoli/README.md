### _E. coli_ analysis

The original E. coli analysis was based on the Vigue 2022 dataset, but this was made up of only protein domains. It was also hard to link these to the alphafold-based structure proteins based on the E. coli K12 strain proteome. Accordingly, I re-ran the E. coli allele frequency workflow based on the "allelome" dataset from Catoiu and colleagues (2023): https://doi.org/10.1073/pnas.2218835120.

This dataset is actually much more convenient to process as it provides the segregating codons and allele frequencies directly.

1. Parse allelome files to get single table

Did so locally with this R script:
```
parse_allelome_seg_subs.R
```

1. Infer amino acid preferences

Infer amino acid preferences per site based on consensus sequences for method that takes sequence context into account.

*VespaG*
(Internal note: VespaG was run on SP5000)

Ran commit version 660b17c6964eb6db8e3f3bee2b8bbd3f9f574d23.

Note that option that normalized predictions to be between 0-1 is set.

Also, VespaG needs to be run from the tool directory (currently anyway).

```
cd /mfs/gdouglas/local/prg/VespaG

python -m vespag predict \
	--input /mfs/gdouglas/projects/aa_selection/vespag/ecoli/translated_aligned_consensus.faa \
	--output /mfs/gdouglas/projects/aa_selection/vespag/ecoli/vespag_out \
	--normalize

gzip /mfs/gdouglas/projects/aa_selection/vespag/ecoli/vespag_out/*csv
```

Then (in the same folder as in earlier commands) to get the VespaG predictions in preference-format:
```
python ~/scripts/aa_distance_explore/compute_prefs/vespag_to_pref.py \
	-i Ecoli_focal_seqs_vespag_output \
	-o prefs/vespag \
	--ref_fill NA
```

1. Mean codon exchangeabilities per site

```
mkdir codon_e_tmp
for CODON in strain_data/codon_aligned/*fna; do
    FASTABASE=$( basename $CODON .fna )

    echo "python ~/scripts/aa_distance_explore/fasta_processing/protein_site_exchangeability.py \
            --id $FASTABASE \
            -c strain_data/codon_aligned_consensus.fna.gz \
            -a strain_data/translated_aligned_consensus.faa.gz \
            -v ~/projects/aa_distance/ecoli_variants/prefs/vespag/ \
            --codon_alignments strain_data/codon_aligned \
            -o codon_e_tmp/$FASTABASE.tsv" >> protein_site_exchangeability_cmds.sh
done

cat protein_site_exchangeability_cmds.sh | parallel -j 32 --joblog protein_site_exchangeability_cmds.log

python ~/local/parallel_joblog_summary/joblog_summary.py \
            --cmds protein_site_exchangeability_cmds.sh \
            --log protein_site_exchangeability_cmds.log
```

Get header (by re-running command with the `--header` option too) and combine individual files.
```
python ~/scripts/aa_distance_explore/fasta_processing/protein_site_exchangeability.py \
            --id $FASTABASE \
            -c strain_data/codon_aligned_consensus.fna.gz \
            -a strain_data/translated_aligned_consensus.faa.gz \
            -v ~/projects/aa_distance/ecoli_variants/prefs/vespag/ \
            --codon_alignments strain_data/codon_aligned \
            --header \
            -o HEADER_tmp.tsv
head -n 1 HEADER_tmp.tsv > strain_data/per_codon_vespag_exchangeability.tsv
cat codon_e_tmp/*tsv >> strain_data/per_codon_vespag_exchangeability.tsv
rm -r codon_e_tmp HEADER_tmp.tsv
gzip strain_data/per_codon_vespag_exchangeability.tsv
```

1. Create combined tables

Get combined table of substitutions and VespaG preferences.
```
python ~/scripts/aa_distance_explore/fasta_processing/combine_subs_and_prefs.py \
    -s strain_data/ecoli_seg_subs.tsv.gz \
    -f prefs/ \
    -o strain_data/ecoli_seg_subs_w_vespag.tsv

gzip strain_data/ecoli_seg_subs_w_vespag.tsv
```

Also get a table of mean exchangeability per site per preference type,
split by whether the site is always invariant or is segregating (at any frequency).
This is a hard-coded script.
```
python /home6/gmdougla/scripts/aa_distance_explore/allele_freq_vs_predicted_effects/ecoli/ecoli_subset_invariant_and_subs.py > \
    strain_data/ecoli_per_codon_vespag_exchange_invariant_vs_freq.tsv

gzip strain_data/ecoli_per_codon_vespag_exchange_invariant_vs_freq.tsv
```
