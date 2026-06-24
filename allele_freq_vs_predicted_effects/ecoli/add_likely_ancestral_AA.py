import os
import gzip
import sys
import pandas as pd
import csv

prot_id_map = pd.read_csv('/Users/gavin/Drive/research/aa_distance/data/Ecoli_focal_seq_id_map.txt',
                          sep=' ', header=None, index_col=1)

intab = pd.read_csv('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/parsed_allelome_variants_w_prefs_dssp.tsv.gz',
                    sep='\t', compression='gzip', header=0)

intab['Blattner'] = intab['protein'].map(prot_id_map[0])

intab['Eferg_AA'] = 'NA'
intab['Sent_AA'] = 'NA'
intab['Likely_Anc_AA'] = 'NA'

unique_blattner = intab['Blattner'].unique()

prepped_df = {}

for bnumber in unique_blattner:
    subtab = intab[intab['Blattner'] == bnumber]

    subtab.index = subtab['pos']

    allelome_filepath = "/Users/gavin/Drive/research/aa_distance/polarize_allelome/data/allelome_w_relatives/" + bnumber + ".csv.gz"
    with gzip.open(allelome_filepath, 'rt') as in_fh:
        headerline = in_fh.readline().strip().split(',')
        consensus_aa_idx = headerline.index('consensus_AAseq')
        Eferg_aa_idx = headerline.index('Eferg_AA')
        Sent_aa_idx = headerline.index('Sent_AA')

        pos = -1
        for line in in_fh:
            fields = next(csv.reader([line]))
            consensus_aa = fields[consensus_aa_idx]

            if consensus_aa == '-' or consensus_aa == '*':
                continue
            pos += 1
            Eferg_aa = fields[Eferg_aa_idx]
            Sent_aa = fields[Sent_aa_idx]

            if pos in subtab.index:
                # Ensure only one index matches this pos.
                if isinstance(subtab.loc[pos], pd.DataFrame):
                    sys.exit(f"Error: multiple entries in subtab for {bnumber} at position {pos}")
                alt_aa = subtab.loc[pos, 'alt_aa']

                if consensus_aa != subtab.loc[pos, 'consensus_aa']:
                    sys.exit(f"Error: consensus AA {consensus_aa} does not match expected AA {subtab.loc[pos, 'consensus_aa']} at position {pos} for {bnumber}")
                else:
                    subtab.loc[pos, 'Eferg_AA'] = Eferg_aa
                    subtab.loc[pos, 'Sent_AA'] = Sent_aa

                    if Eferg_aa == Sent_aa:
                        subtab.loc[pos, 'Likely_Anc_AA'] = Eferg_aa
                    else:
                        subtab.loc[pos, 'Likely_Anc_AA'] = "NA"

    prepped_df[bnumber] = subtab

combined_df = pd.concat(prepped_df.values())
combined_df.to_csv('~/Drive/research/aa_distance/aa_distance_zenodo/allele_freq_vs_predicted_effects/ecoli_variants/parsed_allelome_variants_w_prefs_dssp_w_anc.tsv.gz',
                   index=False, compression='gzip', sep='\t')
