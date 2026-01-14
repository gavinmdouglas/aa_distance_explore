#!/usr/bin/python3

import argparse
import sys
import os
import pandas as pd
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import read_fasta, code11_codon_to_aa, possible_aa


def main():

    parser = argparse.ArgumentParser(
        description='''
        Determine amino acid mean exchangeability for protein consensus sequences,
        based on (usually) VespaG preferences per site. Will run one protein ID at a time, to save time.

        NOTE: An earlier version of this script used codon exchangeability tables with pre-calcualted values. But this produced a lot of redundant lines...
        Can just cross-reference these static values rather than repeating them for every position.
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-i', '--id',
                        metavar='PATH',
                        type=str,
                        help='Protein ID.',
                        required=True)

    parser.add_argument('-c', '--codon',
                        metavar='PATH',
                        type=str,
                        help='File with codon consensus sequences.',
                        required=True)

    parser.add_argument('-a', '--aa',
                        metavar='PATH',
                        type=str,
                        help='Path to FASTA of amino acid consensus sequences.',
                        required=True)

    # Path to codon alignments.
    parser.add_argument('--codon_alignments',
                        metavar='PATH',
                        type=str,
                        help='Path to codon alignments. Used to get the number of unambiguously covered sequences per site.',
                        required=True)

    parser.add_argument('-v', '--vespag',
                        metavar='PATH',
                        type=str,
                        help='Path to VespaG folder.',
                        required=True)

    parser.add_argument('-o', '--output',
                        metavar='OUTPUT',
                        type=str,
                        help='Path to output table.',
                        required=True)

    parser.add_argument('--header',
                        action='store_true',
                        help='Print header.')

    args = parser.parse_args()

    seq_id = args.id
    codon_seqs = read_fasta(args.codon)
    gene = codon_seqs[seq_id]
    aa_seqs = read_fasta(args.aa)
    protein = aa_seqs[seq_id]

    vespag_file = os.path.join(args.vespag, seq_id + '.csv')
    vespag_prefs = pd.read_csv(vespag_file, sep=',', index_col=0, header=0)

    # Get number of sequences with non-ambig codons per pos.
    cov_per_pos = defaultdict(int)
    aligned_codons = read_fasta(os.path.join(args.codon_alignments, seq_id + '.fna'))
    gene_len = len(gene)
    for gene_seq in aligned_codons.values():
        if len(gene_seq) != gene_len:
            sys.exit('Error: Alignment length does not match gene length.')
        for i in range(len(protein)):
            tmp_codon = gene_seq[i * 3: (i * 3) + 3]
            if tmp_codon in code11_codon_to_aa.keys() and code11_codon_to_aa[tmp_codon] in possible_aa:
                cov_per_pos[i] += 1

    with open(args.output, 'w') as out_fh:
        if args.header:
            print('\t'.join(['Protein', 'Pos', 'Consensus_codon', 'Consensus_AA', 'N'] + ['VespaG']), file=out_fh)

        bases = ['A', 'C', 'G', 'T']

        for i in range(len(protein)):
            aa1 = protein[i]
            codon = gene[i * 3: (i * 3) + 3]

            outline = [seq_id, str(i), codon, aa1, str(cov_per_pos[i])]

            if code11_codon_to_aa[codon] not in possible_aa:
                continue
            elif aa1 != code11_codon_to_aa[codon]:
                print(i)
                print(codon)
                print('Expected AA: ' + aa1)
                print('Codon AA: ' + code11_codon_to_aa[codon])
                sys.exit('Error: codon and amino acid do not match.')

            # And then get VespaG.
            num_nonsyn_single = 0
            similarity_sum = 0.0
            for j in range(3):
                for base in bases:
                    if codon[j] == base:
                        continue
                    new_codon = codon[:j] + base + codon[j + 1:]
                    aa2 = code11_codon_to_aa[new_codon]
                    if aa1 == aa2 or aa2 == '*':
                        continue
                    else:
                        num_nonsyn_single += 1
                        similarity_sum += vespag_prefs.loc[i + 1, aa2]

            if num_nonsyn_single > 0:
                outline.append(str(similarity_sum / num_nonsyn_single))
            else:
                outline.append('NA')

            print('\t'.join(outline), file=out_fh)


if __name__ == '__main__':
    main()
