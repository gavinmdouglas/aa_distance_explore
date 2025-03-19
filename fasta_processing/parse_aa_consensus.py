#!/usr/bin/python3

import argparse
import sys
import os
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import read_fasta, write_fasta, possible_aa


def main():

    parser = argparse.ArgumentParser(
        description='''
        Parse FASTA with aligned amino acids for a single protein.
        Determine consensus amino acid sequence and output.
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-i', '--input',
                        metavar='PATH',
                        type=str,
                        help='Path to single amino acid alignment.',
                        required=True)

    parser.add_argument('-o', '--output',
                        metavar='OUTPUT',
                        type=str,
                        help='Path to output FASTA with protein consensus sequence.',
                        required=True)

    args = parser.parse_args()

    seqs = read_fasta(args.input)
    file_base = os.path.basename(args.input).split('.')[0]

    consensus_seq = ''

    aa_tallies_by_pos = defaultdict(lambda: defaultdict(int))
    seq_len = len(list(seqs.values())[0])

    for seq_id in seqs:
        if len(seqs[seq_id]) != seq_len:
            print(file_base, file=sys.stderr)
            print(seq_id, file=sys.stderr)
            sys.exit('Sequence length does not match. Expected ' + str(seq_len) + ' bp, but found ' + str(len(seqs[seq_id])) + ' bp.')

        for i, aa in enumerate(seqs[seq_id]):
            aa_tallies_by_pos[i][aa] += 1

    for i in range(seq_len):
        aa_tally = aa_tallies_by_pos[i]
        consensus_aa = max(aa_tally, key=aa_tally.get)
        if consensus_aa not in possible_aa:
            sys.exit('Error: Non-standard amino acid found as consensus sequence: ' + consensus_aa)
        consensus_seq += consensus_aa

    if len(consensus_seq) != seq_len:
        sys.exit('Error: Consensus sequence length does not match input sequence length.')

    output_dict = {file_base: consensus_seq}
    write_fasta(output_dict, args.output)


if __name__ == '__main__':
    main()
