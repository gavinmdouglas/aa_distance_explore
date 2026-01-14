#!/usr/bin/env python

from collections import defaultdict
import gzip
import sys

# Subset get combined table of all subs (of any frequency) and invariant sites,
# with at least 60,000 strains unambiguously represented at that site.
# Also, as a sanity check, make sure the N values match between the tables.

# First, read through sub table.
freq_subs = set()
freq_sub_maf = defaultdict(list)
freq_sub_AC = defaultdict(list)
freq_sub_N = dict()
freq_sub_consensus_aa = dict()

with gzip.open('/home6/gmdougla/projects/aa_distance/ecoli_variants/strain_data/ecoli_seg_subs.tsv.gz', 'rt') as subs_fh:
    sub_header = subs_fh.readline().strip().split('\t')
    col_to_i = {}
    for i, col in enumerate(sub_header):
        col_to_i[col] = i
    for subs_fh_line in subs_fh:
        subs_fh_line = subs_fh_line.strip().split('\t')
        position = subs_fh_line[col_to_i['protein']] + ';' + subs_fh_line[col_to_i['pos']]
        freq = float(subs_fh_line[col_to_i['AF']])
        AC = int(subs_fh_line[col_to_i['AC']])
        freq_subs.add(position)
        freq_sub_maf[position].append(freq)
        freq_sub_AC[position].append(AC)
        freq_sub_N[position] = int(subs_fh_line[col_to_i['N']])
        freq_sub_consensus_aa[position] = subs_fh_line[col_to_i['consensus_aa']]

# Then move through mean exchangeabilitity table.
with gzip.open('/home6/gmdougla/projects/aa_distance/ecoli_variants/strain_data/per_codon_vespag_exchangeability.tsv.gz', 'rt') as ex_fh:
    header = ex_fh.readline().strip().split('\t')
    col_to_i = {}
    for i, col in enumerate(header):
        col_to_i[col] = i

    new_header = header[:4] + ['Max_AC', 'Max_AF'] + header[4:]
    print('\t'.join(new_header))
    for line in ex_fh:
        line = line.strip().split('\t')
        if int(line[col_to_i['N']]) < 60000:
            continue

        position = line[col_to_i['Protein']] + ';' + line[col_to_i['Pos']]
        if position in freq_subs:
            if freq_sub_N[position] != int(line[col_to_i['N']]):
                print(freq_sub_N[position], line[col_to_i['N']], file=sys.stderr)
                print('Warning: N values do not match for position ' + position + '\n', file=sys.stderr)

            if freq_sub_consensus_aa[position] != line[col_to_i['Consensus_AA']]:
                print('Warning: Consensus AA values do not match for position ' + position + '\n', file=sys.stderr)

            max_AC = max(freq_sub_AC[position])
            max_freq = max(freq_sub_maf[position])
            print('\t'.join(line[:4] + [str(max_AC), str(max_freq)] + line[4:]))
        else:
            print('\t'.join(line[:4] + ['0', '0'] + line[4:]))
