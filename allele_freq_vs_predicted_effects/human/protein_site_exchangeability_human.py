#!/usr/bin/python3

import argparse
import sys
import os
import pandas as pd
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from functions import read_fasta, possible_aa, read_genetic_code

code1_codon_to_aa, possible_start_codons = read_genetic_code('/home6/gmdougla/db/codon_tables/genetic_codes/code1.tsv')

def main():
    parser = argparse.ArgumentParser(
        description='''
        Determine amino acid mean exchangeability for protein consensus sequences,
        (1) VespaG and (2) RaSP preferences per site.
        Will run one protein ID at a time, to save time.

        Originally included all exchangeability metrics, but this lead to unnecessary redundancy.

        Modified for human SNVs, as these are in a specific format.
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

    parser.add_argument('-v', '--vespag',
                        metavar='PATH',
                        type=str,
                        help='Path to VespaG folder.',
                        required=True)
    
    parser.add_argument('-r', '--rasp',
                        metavar='PATH',
                        type=str,
                        help='Path to RaSP folder.',
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

    rasp_file = os.path.join(args.rasp, seq_id + '.csv')
    rasp_prefs = pd.read_csv(rasp_file, sep=',', index_col=0, header=0)

    with open(args.output, 'w') as out_fh:
        if args.header:
            print('\t'.join(['protein', 'pos', 'consensus_codon', 'consensus_aa', 'vespag', 'rasp']), file=out_fh)

        bases = ['A', 'C', 'G', 'T']

        for i in range(len(protein)):
            aa1 = protein[i]
            codon = gene[i * 3: (i * 3) + 3]
            outline = [seq_id, str(i), codon, aa1]


            if codon not in code1_codon_to_aa.keys() or code1_codon_to_aa[codon] not in possible_aa:
                continue
            elif aa1 != code1_codon_to_aa[codon]:
                if i == 0 and aa1 == 'M':
                    # Alternative start codon - skip.
                    continue
                print(i)
                print(codon)
                print('Expected AA: ' + aa1)
                print('Codon AA: ' + code1_codon_to_aa[codon])
                sys.exit('Error: codon and amino acid do not match.')

            # For VespaG.
            num_nonsyn_single = 0
            similarity_sum = 0.0
            for j in range(3):
                for base in bases:
                    if codon[j] == base:
                        continue
                    new_codon = codon[:j] + base + codon[j + 1:]
                    aa2 = code1_codon_to_aa[new_codon]
                    if aa1 == aa2 or aa2 == '*':
                        continue
                    else:
                        num_nonsyn_single += 1
                        similarity_sum += vespag_prefs.loc[i + 1, aa2]
            outline.append(str(similarity_sum / num_nonsyn_single))

            # For rasp.
            num_nonsyn_single = 0
            similarity_sum = 0.0
            for j in range(3):
                for base in bases:
                    if codon[j] == base:
                        continue
                    new_codon = codon[:j] + base + codon[j + 1:]
                    aa2 = code1_codon_to_aa[new_codon]
                    if aa1 == aa2 or aa2 == '*':
                        continue
                    else:
                        num_nonsyn_single += 1
                        similarity_sum += np.power(10, rasp_prefs.loc[i + 1, aa2] * -1)
            outline.append(str(similarity_sum / num_nonsyn_single))

            print('\t'.join(outline), file=out_fh)


if __name__ == '__main__':
    main()
