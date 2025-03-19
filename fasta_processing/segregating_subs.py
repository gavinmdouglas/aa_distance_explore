#!/usr/bin/python3

import argparse
import sys
import os
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import read_fasta, possible_aa


def main():

    parser = argparse.ArgumentParser(
        description='''
        Call segregating substitutions in a FASTA file relative to a consensus sequence.
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-c', '--consensus',
                        metavar='PATH',
                        type=str,
                        help='FASTA containing consensus sequence for the AAs in the FASTA (others can be there too). The sequence header must match the FASTA prefix exactly.',
                        required=True)

    parser.add_argument('-f', '--fasta',
                        metavar='PATH',
                        type=str,
                        help='Path to FASTA with aligned amino acids for single protein.',
                        required=True)

    parser.add_argument('-o', '--output',
                        metavar='OUTPUT',
                        type=str,
                        help='Path to output seg. sub. table for this protein.',
                        required=True)

    parser.add_argument('--header',
                        action='store_true',
                        help='Print header.')

    args = parser.parse_args()

    consensus_seqs = read_fasta(args.consensus)

    with open(args.output, 'w') as out:
        if args.header:
            out.write('\t'.join(['protein', 'pos', 'consensus_aa', 'seg_aa', 'AC', 'N', 'AF']) + '\n')

        aa_tallies_by_pos = defaultdict(lambda: defaultdict(int))
        file_base = os.path.basename(args.fasta).split('.')[0]

        if file_base not in consensus_seqs:
            sys.exit('Error: Consensus sequence not found for ' + file_base)

        seqs = read_fasta(args.fasta)

        seq_len = len(list(seqs.values())[0])

        for seq_id in seqs:
            if len(seqs[seq_id]) != seq_len:
                print(file_base, file=sys.stderr)
                print(seq_id, file=sys.stderr)
                sys.exit('Sequence length does not match. Expected ' + str(seq_len) + ' bp, but found ' + str(len(seqs[seq_id])) + ' bp.')

            for i, aa in enumerate(seqs[seq_id]):
                aa_tallies_by_pos[i][aa] += 1

        consensus_seq = consensus_seqs[file_base]
        for i in range(seq_len):
            aa_tally = aa_tallies_by_pos[i]
            consensus_aa = consensus_seq[i]
            if consensus_aa not in possible_aa:
                continue

            if len(aa_tally) == 0:
                sys.exit('No amino acids found at position ' + str(i) + ' in ' + file_base + '?!')
            elif len(aa_tally) == 1:
                if list(aa_tally.keys())[0] != consensus_aa:
                    sys.exit('Only one amino acid found at position ' + str(i) + ' in ' + file_base + '. Expected ' + consensus_aa + ', but found ' + list(aa_tally.keys())[0])
                else:
                    continue
            else:
                seg_aa = set()
                num_strains = 0
                max_aa = max(aa_tally, key=aa_tally.get)
                if max_aa != consensus_aa:
                    sys.exit('Error: Most common amino acid does not match consensus for ' + file_base + ' at position ' + str(i) + '.')
                for aa in aa_tally:
                    if aa not in possible_aa:
                        continue
                    num_strains += aa_tally[aa]
                    if aa != consensus_aa and aa in possible_aa:
                        seg_aa.add(aa)

                for seg in seg_aa:
                    out.write('\t'.join([file_base, str(i), consensus_aa, seg, str(aa_tally[seg]), str(num_strains), str(aa_tally[seg] / num_strains)]) + '\n')


if __name__ == '__main__':
    main()
