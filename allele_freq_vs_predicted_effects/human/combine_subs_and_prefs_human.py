#!/usr/bin/python3

import sys
import os
from collections import defaultdict
import pandas as pd
from pyfaidx import Fasta

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from functions import possible_aa, read_fasta

# Hard-coded script specifically for Human SNVs, as these are in a specific format.
# Keep only most frequent substitution per protein position.

# Read in human genome (and make chromosome names match those in the SNV table).
ref = Fasta('/home6/gmdougla/projects/aa_distance/human_variants/full_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa')

def get_context_base(chrom, pos, offset):
    i = pos - 1 + offset
    if i < 0 or i >= len(ref[chrom]):
        return 'N'
    return str(ref[chrom][i]).upper()

prot_seqs = read_fasta('/home6/gmdougla/projects/aa_distance/human_variants/nonsyn_snvs/gencode_w_snv.faa')
considered_proteins = set(list(prot_seqs.keys()))

snvfiles = os.listdir('/home6/gmdougla/projects/aa_distance/human_variants/nonsyn_snvs/by_chrom/')
snv_tabs = []
for snvfile in snvfiles:

    snv_tab = pd.read_csv('/home6/gmdougla/projects/aa_distance/human_variants/nonsyn_snvs/by_chrom/' + snvfile, sep='\t', header=0)
    snv_tabs.append(snv_tab)

snv_tab = pd.concat(snv_tabs, ignore_index=True)

# Remove rows with non-possible AA substitutions.
snv_tab = snv_tab[(snv_tab['ref_aa'].isin(possible_aa)) & (snv_tab['alt_aa'].isin(possible_aa))]

# Keep rows with protein_id in considered_proteins.
snv_tab = snv_tab[snv_tab['protein_id'].isin(considered_proteins)]

# Make sure no multi-base mutations are present.
snv_tab = snv_tab[(snv_tab['ref_base'].str.len() == 1) & (snv_tab['alt_base'].str.len() == 1)]

# Keep most frequent non-synonymous mutations per protein position.
snv_tab = snv_tab.sort_values(by='variant_count', ascending=False).drop_duplicates(subset=['protein_id', 'protein_pos'], keep='first')

# Make chromosome names match those in the reference genome.
snv_tab['chrom'] = snv_tab['chrom'].str.replace('^chr', '', regex=True)

pref_master_folder = '/home6/gmdougla/projects/aa_distance/human_variants/nonsyn_snvs/prefs'
pref_folders = sorted([f for f in os.listdir(pref_master_folder) if os.path.isdir(os.path.join(pref_master_folder, f)) and 'thermoMPNN' not in f])

if len(pref_folders) == 0:
    sys.exit('Error: No preference folders found.')

header_fields = ['protein_id', 'protein_pos', 'ref_base', 'alt_base', 'ref_aa', 'alt_aa', 'variant_count', 'sample_count', 'Ts_or_Tv', 'CpG_or_not'] + pref_folders
print('\t'.join(header_fields))

for focal_prot in sorted(snv_tab['protein_id'].unique()):

    # Read in preference files from ML tools for this protein.
    prefs_per_protein = defaultdict(dict)
    for pref_folder in pref_folders:
        pref_file = os.path.join(pref_master_folder, pref_folder, focal_prot + '.csv')
        if not os.path.isfile(pref_file):
            sys.exit('Error: No preference file found for ' + focal_prot + ' in ' + pref_folder)
        with open(pref_file, 'r') as pref_fh:
            pref_header = pref_fh.readline().strip().split(',')
            for pref_line in pref_fh:
                pref_line = pref_line.strip().split(',')
                pos = int(pref_line[0]) - 1
                prefs_per_protein[pref_folder][pos] = {}
                for i in range(1, len(pref_line)):
                    prefs_per_protein[pref_folder][pos][pref_header[i]] = pref_line[i]

    prot_tab = snv_tab[snv_tab['protein_id'] == focal_prot]

    for _, row in prot_tab.iterrows():
        ref_base = row['ref_base']
        alt_base = row['alt_base']
        protein_pos = int(row['protein_pos'])
        ref_aa = row['ref_aa']
        alt_aa = row['alt_aa']

        chrom = row['chrom']
        genome_pos = int(row['genome_pos'])

        if ref_base != str(ref[chrom][genome_pos - 1]).upper():
            sys.exit('Error: Reference base does not match the reference genome.')

        # Identify as transition or transversion.
        if {ref_base, alt_base} in ({'A', 'G'}, {'C', 'T'}):
            mut_type = 'transition'
        else:
            mut_type = 'transversion'

        # Identify as at a CpG site or not.
        cpg_or_not = 'Not_CpG'
        if ref_base == 'C':
            next_base = get_context_base(chrom, genome_pos, 1)
            if next_base == 'G':
                cpg_or_not = 'CpG'
        elif ref_base == 'G':
            prev_base = get_context_base(chrom, genome_pos, -1)
            if prev_base == 'C':
                cpg_or_not = 'CpG'

        outline = [focal_prot, str(protein_pos), ref_base, alt_base, ref_aa, alt_aa, str(row['variant_count']), str(row['sample_count']), mut_type, cpg_or_not]

        for pref_folder in pref_folders:
            pos = protein_pos - 1
            if pos in prefs_per_protein[pref_folder]:
                # Make sure ref_aa value is empty.
                if 'RvC' not in pref_folder and pref_folder != 'rasp' and prefs_per_protein[pref_folder][pos][ref_aa] != '':
                    print(pref_folder)
                    print(pos)
                    print(focal_prot)
                    print(prefs_per_protein[pref_folder][pos][ref_aa])
                    print(prefs_per_protein[pref_folder][pos])
                    sys.exit('Error: Expected reference amino acid is not empty, which it should be (if it was filled with NA previously)')
                else:
                    outline.append(prefs_per_protein[pref_folder][pos][alt_aa])
            else:
                sys.exit('Error: No preference found for ' + focal_prot + ' at position ' + str(pos) + ' in ' + pref_folder)
        print('\t'.join(outline))
