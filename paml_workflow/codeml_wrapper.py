#!/usr/bin/env python

import os
import argparse
from Bio.Phylo.PAML import codeml


def main():
    parser = argparse.ArgumentParser(
        description='''Run CODEML on a FASTA alignment.''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--fasta',
                        help='Path to FASTA file.',
                        required=True,
                        metavar='PATH',
                        type=str)

    parser.add_argument('-t', '--tree',
                        help='Path to tree file.',
                        required=True,
                        metavar='PATH',
                        type=str)

    parser.add_argument('-o', '--output_dir',
                        help='Path to output folder.',
                        required=True,
                        metavar='PATH',
                        type=str)

    parser.add_argument('-m', '--model',
                        default=0,
                        help='Model for omega estimation in PAML vernacular. Default value is 0 (all branches have the same rate).',
                        required=False,
                        metavar='INT',
                        type=int)

    parser.add_argument('-c', '--CodonFreq',
                        default=2,
                        help='Codon frequency model. Default value is 2 (codon frequencies are estimated from nucleotide frequencies at each position separately).',
                        required=False,
                        metavar='INT',
                        type=int)

    parser.add_argument('-a', '--aaDist',
                        default=0,
                        help='How AA distances are included in model (0 by default - equal distances). 1 indicates Grantham distances should be used.',
                        required=False,
                        metavar='INT',
                        type=int)

    parser.add_argument('--nssites',
                        default=0,
                        help='Number of sites to be used in the analysis. Default value is 0 (all sites).',
                        required=False,
                        metavar='INT',
                        type=int)

    parser.add_argument('--noisy',
                        default=0,
                        help='Noisy output. Default value is 0 (minimal info).',
                        required=False,
                        metavar='INT',
                        type=int)

    parser.add_argument('--verbose',
                        default=0,
                        help='Verbose output. Default value is 1.',
                        required=False,
                        metavar='INT',
                        type=int)

    args = parser.parse_args()

    cml = codeml.Codeml()
    cml.alignment = args.fasta
    cml.tree = args.tree
    cml.out_file = os.path.join(args.output_dir, 'codeml_results.txt')
    cml.working_dir = args.output_dir

    # Indicates codons, one alignment, remove sites with ambiguous/missing data.
    # Universal genetic code, use user-specified tree. No molecular clock.
    cml.set_options(seqtype=1)
    cml.set_options(ndata=1)
    cml.set_options(cleandata=1)
    cml.set_options(icode=0)
    cml.set_options(runmode=0)
    cml.set_options(clock=0)

    cml.set_options(noisy=args.noisy)
    cml.set_options(verbose=args.verbose)
    cml.set_options(model=args.model)
    cml.set_options(NSsites=[args.nssites])
    cml.set_options(CodonFreq=args.CodonFreq)
    cml.set_options(aaDist=args.aaDist)

    cml.run(verbose=True)


if __name__ == '__main__':
    status = main()
