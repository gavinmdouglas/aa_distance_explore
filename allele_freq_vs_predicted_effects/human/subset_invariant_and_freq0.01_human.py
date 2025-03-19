#!/usr/bin/env python

import sys
from collections import defaultdict
import gzip

# Subset get combined table of subs > 1% and invariant sites.

# First, read through sub table.
other_subs = set()
freq_subs = set()
freq_sub_maf = defaultdict(list)
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
        if freq >= 0.01:
            freq_subs.add(position)
            freq_sub_maf[position].append(str(freq))
            freq_sub_consensus_aa[position] = subs_fh_line[col_to_i['ref_aa']]
        else:
            other_subs.add(position)

# Then move through mean exchangeabilitity table.
with gzip.open('/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/per_codon_exchangeability.tsv.gz', 'rt') as ex_fh:
    header = ex_fh.readline().strip().split('\t')
    col_to_i = {}
    for i, col in enumerate(header):
        col_to_i[col] = i

    new_header = header[:3] + ['AF', 'Class'] + header[3:]
    print('\t'.join(new_header))

    for line in ex_fh:
        line = line.strip().split('\t')
        position = line[col_to_i['Protein']] + ';' + line[col_to_i['Pos']]
        if position in freq_subs:
            if freq_sub_consensus_aa[position] != line[col_to_i['Consensus_AA']]:
                print('Warning: Consensus AA values do not match for position ' + position + '\n', file=sys.stderr)

            freq_str = ';'.join(freq_sub_maf[position])
            new_line = line[:3] + [freq_str, 'atleast0.01'] + line[3:]
            print('\t'.join(new_line))

        elif position in other_subs:
            continue
        else:
            new_line = line[:3] + ['NA', 'Invariant'] + line[3:]
            print('\t'.join(new_line))
