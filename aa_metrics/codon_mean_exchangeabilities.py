#!/usr/bin/env python

import sys
import os
from collections import defaultdict
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import (read_sim_matrix,
                       continuous_prefs_dict_sanity,
                       code11_codon_to_aa)

# Compute mean exchangeabilities based on average of amino acid similarities possible for all possible single (nn-synonymous) nucleotide changes.
sim_folder = '/home6/gmdougla/projects/aa_distance/aa_metrics/prepped_similarity_consistent'

sim_files = [f for f in os.listdir(sim_folder) if f.endswith('.tsv.gz')]
if len(sim_files) == 0:
    sys.exit('Error: No similarity effect files found.')

sub_maps = []
tab_ids = []
for sim_file in sorted(sim_files):
    sim_dict = read_sim_matrix(sub_sim_file=os.path.join(sim_folder, sim_file),
                                identity_set='NA')

    dict_check = continuous_prefs_dict_sanity(sim_dict)
    if not dict_check:
        sys.exit('Error: Similarity table is not formatted correctly: ' + os.path.join(sim_folder, sim_file))

    sub_maps.append(sim_dict)
    tab_ids.append(sim_file.split('.')[0])

# Read Grantham distances into sub_map as well (to compute Graur stability index, for comparison).
grantham = pd.read_csv('/home6/gmdougla/projects/aa_distance/aa_metrics/distances/grantham_orig.tsv.gz',
                       sep='\t', header=0, index_col=0)

# Get dictionary of tuple of (row_val, col_val) as key and value as cell.
grantham_dict = {}
for row in grantham.index:
    for col in grantham.columns:
        grantham_dict[(row, col)] = grantham.loc[row, col]
sub_maps.append(grantham_dict)
tab_ids.append('grantham_dist_Graur_stability_index')

bases = ['A', 'C', 'G', 'T']

# Loop through every possible codon.
if len(code11_codon_to_aa) != 64:
    sys.exit('Error: Code11 codon to amino acid mapping is not complete.')

print('codon\tamino_acid\tnum_nonsyn_single\t' + '\t'.join(tab_ids))
for codon in sorted(code11_codon_to_aa.keys()):
    aa1 = code11_codon_to_aa[codon]

    if aa1 == '*':
        continue

    num_nonsyn_single = 0
    similarity_sum = defaultdict(float)
    for i in range(3):
        for base in bases:
            if codon[i] == base:
                continue
            new_codon = codon[:i] + base + codon[i + 1:]
            aa2 = code11_codon_to_aa[new_codon]
            if aa1 == aa2 or aa2 == '*':
                continue
            else:
                num_nonsyn_single += 1
                for sub_map, effect_type in zip(sub_maps, tab_ids):
                    similarity_sum[effect_type] += sub_map[(aa1, aa2)]

    outline = [codon, aa1, str(num_nonsyn_single)]
    for effect_type in tab_ids:
        outline.append(str(similarity_sum[effect_type] / num_nonsyn_single))
    print('\t'.join(outline))
