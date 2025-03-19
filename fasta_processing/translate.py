#!/usr/bin/python3

import argparse
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import(read_fasta,
                      write_fasta,
                      read_genetic_code,
                      translate_DNA)

def main():

    parser = argparse.ArgumentParser(

        description='''
        Translate CDS (DNA) FASTA to AA.
        Requires that a genetic code file from this repository is specified: https://github.com/gavinmdouglas/codon_tables
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter

    )

    parser.add_argument('-f', '--fasta',
                        metavar='PATH',
                        type=str,
                        help='FASTA file containing CDS (nucleotide sequences).',
                        required=True)

    parser.add_argument('-o', '--output',
                        metavar='OUTPUT',
                        type=str,
                        help='Path to output FASTA containing amino acids.',
                        required=True)

    parser.add_argument('-g', '--genetic_code',
                        metavar='PATH',
                        type=str,
                        help='Path to genetic code file.',
                        required=False,
                        default='/home6/gmdougla/db/codon_tables/genetic_codes/code11.tsv')

    parser.add_argument('--skip_check_start',
                        action='store_true',
                        help='Do not check for start codons and throw error if not present at beginning (which is the default behaviour).')
    
    parser.add_argument('--silent',
                        action='store_true',
                        help='Silent mode, no warnings regarding translation to stdout.')

    args = parser.parse_args()

    codon_to_aa, possible_start_codons = read_genetic_code(args.genetic_code)

    cds = read_fasta(args.fasta)

    aa = {}
    for gene in cds:
        aa[gene] = translate_DNA(in_dna=cds[gene],
                                 codon_to_aa=codon_to_aa,
                                 possible_start_codons=possible_start_codons,
                                 check_start=not args.skip_check_start,
                                 silent=args.silent)

        # Remove stop signal at end, if present.
        if aa[gene][-1] == '*':
            aa[gene] = aa[gene][:-1]

    write_fasta(aa, args.output)


if __name__ == '__main__':
    main()
