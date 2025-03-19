Commands run to explore allele frequencies (and invariant sites) for non-synonymous mutations in human proteins.

First, downloaded all gnomAD VCFs and gencode protein FASTA.

Then parsed out Gencode protein sequences that match RaSP sequences. Keep only one matching protein per gene (which resulted in six matching proteins being excluded).
```
# First, to get RaSP protein FASTA.
python /home6/gmdougla/scripts/aa_selection_ms/allele_freq_vs_predicted_effects/human/rasp_human_seqs_to_fasta.py
gzip /home6/gmdougla/projects/aa_selection/prerun_struc_prefs/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.fasta

# Then to run main command.
python /home6/gmdougla/scripts/aa_selection_ms/allele_freq_vs_predicted_effects/human/rasp_seqs_vs_gencode.py
```

That command produced 'gencode.v45.pc_translations_rasp.fa' below, which contains the Gencode proteins of interest (intersecting with RaSP predictions).

Then parsed out non-synonymous mutations per chromosome:
```
cd ~/projects/aa_selection/human_variants/nonsyn_snvs/
mkdir by_chrom

for VCF in ~/projects/aa_selection/human_variants/gnomAD_v4.1/vcfs/*.vcf.bgz; do
	BASENAME=$( basename $VCF .vcf.bgz )
	CHR=${BASENAME/gnomad.exomes.v4.1.sites./}
	echo "python ~/scripts/aa_selection_ms/allele_freq_vs_predicted_effects/human/parse_gnomAD_vcf.py \
	        -v $VCF \
	        -f ~/projects/aa_selection/human_variants/gencode/gencode.v45.pc_translations_rasp.fa.gz \
	        > ~/projects/aa_selection/human_variants/nonsyn_snvs/by_chrom/$CHR.tsv" \
        >> parse_gnomAD_cmds.sh
done

cat parse_gnomAD_cmds.sh | parallel -j 32 --joblog parse_gnomAD_cmds.log '{}'

python ~/local/parallel_joblog_summary/joblog_summary.py --cmds parse_gnomAD_cmds.sh --log parse_gnomAD_cmds.log

parallel -j 24 --eta 'gzip {}' ::: by_chrom/*tsv
```

Parse the non-synonymous mutations identified, output a new FASTA of the protein sequences of all proteins with at least one non-synonymous mutation called (written to 'gencode_w_snv.faa').
Also, output a FASTA with the corresponding gencode CDS (nucleotide) sequences for these proteins, which will be downloaded (written to 'gencode_w_snv.fna').

These files were produced with this command:
```
python /home6/gmdougla/scripts/aa_selection_ms/allele_freq_vs_predicted_effects/human/prep_gencode_fastas_w_snvs.py
```


Infer amino acid preferences per site based on consensus sequences.

*All similarity metrics*

Get amino acid preferences based on all 'similarity' metrics.

```
# Split into 32 chunks to make it easily parallelizable:
python /home6/gmdougla/scripts/aa_selection_ms/fasta_processing/split_fasta.py -f gencode_w_snv.faa -n 32 -o gencode_w_snv_split

# Quick sanity check:
/home6/gmdougla/scripts/aa_selection_ms/allele_freq_vs_predicted_effects/human/split_seqs_sanity_checks.py

mkdir prefs_by_chunk
for FASTA in gencode_w_snv_split/*faa; do
	CHUNK=$( basename $FASTA .faa )
	echo "python ~/scripts/aa_selection_ms/compute_prefs/fasta_to_prefs.py \
				-f $FASTA \
				-o prefs_by_chunk/$CHUNK \
				-r /home6/gmdougla/projects/aa_selection/aa_metrics/prepped_RvC \
				-s /home6/gmdougla/projects/aa_selection/aa_metrics/prepped_similarity \
				-p 0.05 \
				-i NA" >> fasta_to_prefs_cmds.sh
done

cat fasta_to_prefs_cmds.sh | parallel -j 32 --joblog fasta_to_prefs_cmds.log

python ~/local/parallel_joblog_summary/joblog_summary.py \
            --cmds fasta_to_prefs_cmds.sh \
            --log fasta_to_prefs_cmds.log

# Then move them all into one folder.
mkdir prefs
for METRIC in prefs_by_chunk/split_0/*; do
	METRIC=$( basename $METRIC )
	mkdir prefs/$METRIC
	mv prefs_by_chunk/*/$METRIC/* prefs/$METRIC
done
```

VespaG threw error due to X residues. Ran this command to convert all X residues to M (which were in the beginning of sequence). Note that predictions for these sites should be skipped, just added to get it working. Also replaced all 'U' residues with 'C' (which again was just to get VespaG running, and predictions for these sites were ignored in downstream analyses). Written to `gencode_w_snv_noambig.faa`.

```
python ~/allele_freq_vs_predicted_effects/human/check_ambig_and_replace.py
```

Then ran VespaG.

```
# Internal note: On SP5000
cd /mfs/gdouglas/local/prg/VespaG

python -m vespag predict \
	--input /mfs/gdouglas/projects/aa_selection/vespag/human/gencode_w_snv_noambig.faa \
	--output /mfs/gdouglas/projects/aa_selection/vespag/human/vespag_out \
	--normalize

parallel -j 60 --eta 'gzip {}' ::: /mfs/gdouglas/projects/aa_selection/vespag/human/vespag_out/*csv
```

Then (in the same folder as in earlier commands) to get the VespaG predictions in preference-format:
```
python ~/scripts/aa_selection_ms/compute_prefs/vespag_to_pref.py \
	-i vespag_out \
	-o prefs/vespag \
	--ref_fill NA
```


Mean codon exchangeabilities per site

```
grep '>' gencode_w_snv.fna | sed 's/>//g' > gencode_w_snv_ids.txt

mkdir codon_e_tmp
for PROTID in $( cat gencode_w_snv_ids.txt ); do
    echo "python ~/scripts/aa_selection_ms/allele_freq_vs_predicted_effects/human/protein_site_exchangeability_human.py \
            --id $PROTID \
            -c gencode_w_snv.fna \
			-a gencode_w_snv.faa \
            -e ~/projects/aa_selection/aa_metrics/codon_exchangeabilities.tsv.gz \
            -v prefs/vespag/ \
			-r prefs/rasp/ \
            -o codon_e_tmp/$PROTID.tsv" \
		>> protein_site_exchangeability_cmds.sh
done

cat protein_site_exchangeability_cmds.sh | parallel -j 32 --joblog protein_site_exchangeability_cmds.log

python ~/local/parallel_joblog_summary/joblog_summary.py \
            --cmds protein_site_exchangeability_cmds.sh \
            --log protein_site_exchangeability_cmds.log
```

Get header (by re-running command with the `--header` option too) and combine individual files.
```
python ~/scripts/aa_selection_ms/allele_freq_vs_predicted_effects/human/protein_site_exchangeability_human.py \
	--id ENSP00000005178.5 \
	-c gencode_w_snv.fna \
	-a gencode_w_snv.faa \
	-e ~/projects/aa_selection/aa_metrics/codon_exchangeabilities.tsv.gz \
	-v prefs/vespag/ \
	-r prefs/rasp/ \
	--header \
	-o HEADER_tmp.tsv
head -n 1 HEADER_tmp.tsv > per_codon_exchangeability.tsv
cat codon_e_tmp/*tsv >> per_codon_exchangeability.tsv
rm -r codon_e_tmp HEADER_tmp.tsv
```

Create combined tables

Get combined table of substitutions and preferences.
```
python ~/scripts/aa_selection_ms/allele_freq_vs_predicted_effects/human/combine_subs_and_prefs_human.py \
	> ~/projects/aa_selection/human_variants/nonsyn_snvs/human_snvs_w_prefs.tsv
```

Also get a table of mean exchangeability per site per preference type,
split by whether the site is always invariant or is segregating (at any frequency).
This is a hard-coded script.  
```
python /home6/gmdougla/scripts/aa_selection_ms/allele_freq_vs_predicted_effects/human/human_subset_invariant_and_freq.py \
	> human_per_codon_exchangeability_invariant_vs_freq.tsv
```
