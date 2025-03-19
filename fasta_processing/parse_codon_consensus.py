#!/usr/bin/python3

import argparse
import sys
import os
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import (read_fasta,
                       write_fasta,
                       code11_codon_to_aa)


def main():

    parser = argparse.ArgumentParser(
        description='''
        Parse file with codon aligned FASTAs.
        Also parse a single FASTA with the consensus amino acid sequences for each protein.
        Sequence IDs (codon FASTA prefix and consensus amino acid header) must match.
        Determine consensus codon sequence for each protein.
        Output a single FASTA containing all codon consensus sequences.
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-c', '--codon',
                        metavar='PATH',
                        type=str,
                        help='File with codon-aligned for a single protein.',
                        required=True)

    parser.add_argument('-a', '--aa',
                        metavar='PATH',
                        type=str,
                        help='Path to FASTA of amino acid consensus sequences for the same protein',
                        required=True)

    parser.add_argument('-o', '--output',
                        metavar='OUTPUT',
                        type=str,
                        help='Path to output FASTA for consensus sequence.',
                        required=True)

    args = parser.parse_args()

    codon_seqs = read_fasta(args.codon)
    seq_id = os.path.basename(args.codon).split('.')[0]
    aa_consensus_all = read_fasta(args.aa)
    aa_consensus = aa_consensus_all[seq_id]
    codon_consensus = ''
    if '-' in aa_consensus:
        sys.exit('Gap in AA consensus sequence for: ' + seq_id)

    aa_len = len(aa_consensus)
    codon_tallies_by_pos = defaultdict(lambda: defaultdict(int))
    for codon_seq in codon_seqs.values():
        if len(codon_seq) != aa_len * 3:
            sys.exit('Sequence length does not match. Expected ' + str(aa_len * 3) + ' bp, but found ' + str(len(codon_seq)) + ' bp.')

        for i in range(aa_len):
            codon = codon_seq[i * 3: (i * 3) + 3]
            codon_tallies_by_pos[i][codon] += 1

    for i in range(aa_len):
        codon_tally = codon_tallies_by_pos[i]
        codon_tally_by_aa = defaultdict(int)
        codon_aa_groups = defaultdict(list)
        for codon in codon_tally.keys():
            if codon in code11_codon_to_aa.keys():
                aa = code11_codon_to_aa[codon]
            else:
                aa = 'X'
            codon_tally_by_aa[aa] += codon_tally[codon]
            codon_aa_groups[aa].append(codon)
        max_aa = max(codon_tally_by_aa, key=codon_tally_by_aa.get)
        if max_aa != aa_consensus[i]:
            for f in codon_tally.keys():
                print(f + ' ' + str(codon_tally[f]) + ' ' + code11_codon_to_aa[f])
            print('Determined consensus codon AA: ' + max_aa)
            print('Expected consensus AA: ' + aa_consensus[i])
            sys.exit('Error: Most common amino acid does not match consensus for ' + seq_id + ' at position ' + str(i) + '.')
        top_codon = max(codon_aa_groups[max_aa], key=codon_tally.get)
        if code11_codon_to_aa[top_codon] == aa_consensus[i]:
            codon_consensus += top_codon
        else:
            sys.exit('This should have been caught earlier...')

    # Make sure final sequence is the expected length:
    if len(codon_consensus) != aa_len * 3:
        sys.exit('Error: Expected ' + str(aa_len * 3) + ' bp, but found ' + str(len(codon_consensus)) + ' bp.')

    # Also make sure all translated codons match the consensus AA sequence:
    for i in range(aa_len):
        codon = codon_consensus[i * 3: (i * 3) + 3]
        aa = code11_codon_to_aa[codon]
        if aa != aa_consensus[i]:
            sys.exit('Error: Translated codon does not match consensus AA sequence for ' + seq_id + ' at position ' + str(i) + '.')

    out_dict = {seq_id: codon_consensus}
    write_fasta(out_dict, args.output)


if __name__ == '__main__':
    main()
