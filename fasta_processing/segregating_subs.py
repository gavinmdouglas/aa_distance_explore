#!/usr/bin/python3

import argparse
import sys
import os
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import (read_fasta,
                       possible_aa,
                       read_genetic_code,
                       translate_DNA)


def main():

    parser = argparse.ArgumentParser(
        description='''
        Call segregating AA substitutions in a FASTA file relative to a consensus sequence.
        Call the most frequent segregating substitution per site (so only one per codon site).
        Do so based on DNA codon alignments, so that independent nucleotide mutations are considered as different AA substitutions.
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-c', '--consensus',
                        metavar='PATH',
                        type=str,
                        help='FASTA containing consensus sequence for protein-coding DNA sequences in the FASTA (others can be there too). The sequence header must match the FASTA prefix exactly.',
                        required=True)

    parser.add_argument('-f', '--fasta',
                        metavar='PATH',
                        type=str,
                        help='Path to FASTA with aligned protein-coding DNA sequences for single protein.',
                        required=True)

    parser.add_argument('-o', '--output',
                        metavar='OUTPUT',
                        type=str,
                        help='Path to output seg. sub. table for this protein.',
                        required=True)

    parser.add_argument('-g', '--genetic_code',
                        metavar='PATH',
                        type=str,
                        help='Path to genetic code file.',
                        required=False,
                        default='/home6/gmdougla/db/codon_tables/genetic_codes/code11.tsv')

    parser.add_argument('--header',
                        action='store_true',
                        help='Print header.')

    args = parser.parse_args()

    codon_to_aa, possible_start_codons = read_genetic_code(args.genetic_code)

    consensus_seqs = read_fasta(args.consensus)

    file_base = os.path.basename(args.fasta).split('.')[0]
    if file_base not in consensus_seqs:
        sys.exit('Error: Consensus sequence not found for ' + file_base)

    seqs = read_fasta(args.fasta)
    seq_len = len(list(seqs.values())[0])

    codon_tallies_by_pos = defaultdict(lambda: defaultdict(int))

    consensus_seq = consensus_seqs[file_base]

    consensus_seq_aa = translate_DNA(in_dna=consensus_seq,
                                    codon_to_aa=codon_to_aa,
                                    possible_start_codons=possible_start_codons,
                                    check_start=False,
                                    silent=True)

    for seq_id in seqs:
        if len(seqs[seq_id]) != seq_len:
            print(file_base, file=sys.stderr)
            print(seq_id, file=sys.stderr)
            sys.exit('Sequence length does not match. Expected ' + str(seq_len) + ' bp, but found ' + str(len(seqs[seq_id])) + ' bp.')

        for aa_i, codon_i in enumerate(range(0, seq_len, 3)):
            codon_seq = seqs[seq_id][codon_i:codon_i + 3]

            consensus_aa = consensus_seq_aa[aa_i]
            if consensus_aa not in possible_aa:
                continue

            codon_tallies_by_pos[aa_i][codon_seq] += 1

    with open(args.output, 'w') as out:
        if args.header:
            out.write('\t'.join(['protein', 'pos', 'consensus_codon', 'consensus_aa', 'seg_codon', 'seg_aa', 'AC', 'N', 'AF']) + '\n')

        for aa_i in range(len(consensus_seq_aa)):

            if aa_i not in codon_tallies_by_pos:
                continue

            consensus_aa = consensus_seq_aa[aa_i]
            consensus_codon = consensus_seqs[file_base][aa_i * 3: aa_i * 3 + 3]
        
            codon_tallies = codon_tallies_by_pos[aa_i]

            top_codon_tally = 0
            top_codon = None
            num_strains = 0
            max_codon = max(codon_tallies, key=codon_tallies.get)
            if max_codon != consensus_codon:
                print('Warning: Most common codon does not match consensus for ' + file_base + ' at AA position ' + str(aa_i) + '. Most common codon: ' + max_codon + ', consensus codon: ' + consensus_codon + ')', file=sys.stderr)
                print('Skipping this position.', file=sys.stderr)
                continue

            for codon_seq in codon_tallies:
                num_strains += codon_tallies[codon_seq]

                # Check if this is non-ambiguously different AA.
                if codon_seq not in codon_to_aa:
                    continue

                codon_aa = codon_to_aa[codon_seq]

                if codon_aa == consensus_aa:
                    continue

                if codon_tallies[codon_seq] > top_codon_tally:
                    top_codon_tally = codon_tallies[codon_seq]
                    top_codon = codon_seq

            if top_codon is None:
                continue
            
            # Check if top codon differs by exactly one nucleotide.
            differences = sum(1 for a, b in zip(top_codon, consensus_codon) if a != b)
            if differences != 1:
                continue

            top_codon_aa = codon_to_aa[top_codon]
            out.write('\t'.join([file_base, str(aa_i), consensus_codon, consensus_aa, top_codon, top_codon_aa, str(top_codon_tally), str(num_strains), str(top_codon_tally / num_strains)]) + '\n')


if __name__ == '__main__':
    main()
