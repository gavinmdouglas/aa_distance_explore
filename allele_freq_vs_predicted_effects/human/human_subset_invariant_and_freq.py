#!/usr/bin/env python

import sys
from collections import defaultdict
import gzip

# Subset get combined table of subs and invariant sites.

# First, read through sub table.
freq_subs = set()
freq_sub_maf = defaultdict(list)
freq_subs_AC = defaultdict(list)
freq_sub_consensus_aa = dict()

with gzip.open('/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/human_snvs_w_prefs.tsv.gz', 'rt') as subs_fh:
    sub_header = subs_fh.readline().strip().split('\t')
    col_to_i = {}
    for i, col in enumerate(sub_header):
        col_to_i[col] = i
    for subs_fh_line in subs_fh:
        subs_fh_line = subs_fh_line.strip().split('\t')
        position = subs_fh_line[col_to_i['protein_id']] + ';' + subs_fh_line[col_to_i['protein_pos']]
        freq = float(subs_fh_line[col_to_i['variant_count']]) / float(subs_fh_line[col_to_i['sample_count']])
        freq_subs.add(position)
        freq_sub_maf[position].append(freq)
        freq_subs_AC[position].append(int(subs_fh_line[col_to_i['variant_count']]))
        freq_sub_consensus_aa[position] = subs_fh_line[col_to_i['ref_aa']]

# Then move through mean exchangeabilitity table.
with gzip.open('/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/per_codon_exchangeability.tsv.gz', 'rt') as ex_fh:
    header = ex_fh.readline().strip().split('\t')
    col_to_i = {}
    for i, col in enumerate(header):
        col_to_i[col] = i

    new_header = header[:3] + ['Max_AC', 'Max_AF'] + header[3:]
    print('\t'.join(new_header))

    for line in ex_fh:
        line = line.strip().split('\t')
        position = line[col_to_i['Protein']] + ';' + line[col_to_i['Pos']]
        if position in freq_subs:
            if freq_sub_consensus_aa[position] != line[col_to_i['Consensus_AA']]:
                print('Warning: Consensus AA values do not match for position ' + position + '\n', file=sys.stderr)
            max_freq = max(freq_sub_maf[position])
            max_AC = max(freq_subs_AC[position])
            print('\t'.join(line[:3] + [str(max_AC), str(max_freq)] + line[3:]))
        else:
            print('\t'.join(line[:3] + ['0', '0'] + line[3:]))
