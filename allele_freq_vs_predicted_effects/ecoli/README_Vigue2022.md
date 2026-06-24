### _E. coli_ analysis

This analyis was on previously [prepped protein-coding alignments](https://zenodo.org/records/5774192) across _E. coli_ strains. The downloaded folder was renamed to `Vigue2022_dataset`.

Note that most of the scripts used below are general purpose and are in the `fasta_processing/` directory.

`explore_maf_vs_prefs.Rmd` is a simple exploratory R notebook for some initial analyses.


1. Translate codons

```
mkdir -p strain_data/translated

for FASTA in Vigue2022_dataset/homologous_sequences/local_strains/PF*fasta; do

    BASE=$( basename $FASTA .fasta )

    echo "python ~/scripts/aa_distance_explore/fasta_processing/translate.py \
            -f $FASTA \
            -o strain_data/translated/$BASE.faa" \
            --silent \
            --skip_check_start \
            >> translate_cmds.sh
done

cat translate_cmds.sh | parallel -j 32 --joblog translate_cmds.log

# Sanity check that commands finished successfully.
python ~/local/parallel_joblog_summary/joblog_summary.py \
            --cmds translate_cmds.sh \
            --log translate_cmds.log
```

1. Align translated sequences.

Align translated sequence with MUSCLE.

Most sequences were of the same length, but there were of slightly different length. Ran with muscle (super5 mode for speed).

```
mkdir strain_data/translated_aligned

for FASTA in strain_data/translated/*faa*; do

    FASTAFILE=$( basename $FASTA )

    echo "muscle \
            -super5 $FASTA \
            -output strain_data/translated_aligned/$FASTAFILE" \
            >> align_cmds.sh
done

cat align_cmds.sh | parallel -j 32 --joblog align_cmds.log

python ~/local/parallel_joblog_summary/joblog_summary.py \
            --cmds align_cmds.sh \
            --log align_cmds.log
```

1. Determine consensus amino acid sequences across all genes.

Do this for each individual protein to make it easy to parallelize.

```
mkdir consensus_aa_tmp

for AA in strain_data/translated_aligned/*faa; do
    FASTABASE=$( basename $AA .faa )

    echo "python ~/scripts/aa_distance_explore/fasta_processing/parse_aa_consensus.py \
            -i $AA \
            -o consensus_aa_tmp/$FASTABASE.faa" \
        >> consensus_aa_cmds.sh
done

cat consensus_aa_cmds.sh | parallel -j 32 --joblog consensus_aa_cmds.log

python ~/local/parallel_joblog_summary/joblog_summary.py \
            --cmds consensus_aa_cmds.sh \
            --log consensus_aa_cmds.log

cat consensus_aa_tmp/*faa > strain_data/translated_aligned_consensus.faa
rm -r consensus_aa_tmp
```

1. Codon-aligned nucleotide sequences and get consensus

```
mkdir strain_data/codon_aligned

for AA in strain_data/translated_aligned/*faa; do
    FASTABASE=$( basename $AA .faa )
    NUCL="Vigue2022_dataset/homologous_sequences/local_strains/$FASTABASE.fasta"

    echo "python ~/scripts/aa_distance_explore/fasta_processing/align_codons_from_aligned_aa.py \
	        -a $AA \
	        -n $NUCL \
	        -o strain_data/codon_aligned/$FASTABASE.fna" \
        >> align_codons_from_aligned_aa_cmds.sh
done

cat align_codons_from_aligned_aa_cmds.sh | parallel -j 32 --joblog align_codons_from_aligned_aa_cmds.log

python ~/local/parallel_joblog_summary/joblog_summary.py \
            --cmds align_codons_from_aligned_aa_cmds.sh \
            --log align_codons_from_aligned_aa_cmds.log
```

To get consensus codon sequence (ran per individual protein to speed up):

```
mkdir consensus_codon_tmp
AA="strain_data/translated_aligned_consensus.faa.gz"
for CODON in strain_data/codon_aligned/*fna; do
    FASTABASE=$( basename $CODON .fna )

    echo "python ~/scripts/aa_distance_explore/fasta_processing/parse_codon_consensus.py \
            -c $CODON \
            -a $AA \
            -o consensus_codon_tmp/$FASTABASE.fna" \
        >> parse_codon_consensus_cmds.sh
done

cat parse_codon_consensus_cmds.sh | parallel -j 32 --joblog parse_codon_consensus_cmds.log

python ~/local/parallel_joblog_summary/joblog_summary.py \
            --cmds parse_codon_consensus_cmds.sh \
            --log parse_codon_consensus_cmds.log

cat consensus_codon_tmp/*fna > strain_data/codon_aligned_consensus.fna
rm -r consensus_codon_tmp
```

1. Call segregating amino acid polymorphisms

Call them relative to consensus sequence. Again, run individually to parallelize easily.

Do so at codon level, and only keep most frequency segregating codon per site, restricted to those 

```
mkdir seg_subs_tmp

for FASTA in strain_data/codon_aligned/*fna; do
    FASTABASE=$( basename $FASTA .fna )

    echo "python ~/scripts/aa_distance_explore/fasta_processing/segregating_subs.py \
                -c strain_data/codon_aligned_consensus.fna.gz \
                -f $FASTA \
                -o seg_subs_tmp/$FASTABASE.tsv" \
        >> seg_subs_cmds.sh
done

cat seg_subs_cmds.sh | parallel -j 32 --joblog seg_subs_cmds.log

python ~/local/parallel_joblog_summary/joblog_summary.py \
            --cmds consensus_aa_cmds.sh \
            --log consensus_aa_cmds.log
```

Get header (by re-running command with the `--header` option too) and combine individual files.
```
python ~/scripts/aa_distance_explore/fasta_processing/segregating_subs.py \
    -c strain_data/translated_aligned_consensus.faa.gz \
    -f strain_data/translated_aligned/PF00004-GA4805AA_00441.faa \
    -o seg_subs_header.tsv \
    --header
head -n 1 seg_subs_header.tsv > strain_data/ecoli_seg_subs.tsv
cat seg_subs_tmp/*tsv >> strain_data/ecoli_seg_subs.tsv
rm -r seg_subs_tmp seg_subs_header.tsv
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
	-i vespag_out \
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
