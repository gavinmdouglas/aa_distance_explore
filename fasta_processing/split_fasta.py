#!/usr/bin/env python

import argparse
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import read_fasta, write_fasta


def main():
    parser = argparse.ArgumentParser(

        description=
        '''
        Split sequences in FASTA into N separate files.
        ''',

        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--fasta',
                        help='Input FASTA.',
                        required=True,
                        metavar='PATH',
                        type=str)

    parser.add_argument('-o', '--out_folder',
                        help='Output folder for output FASTAs.',
                        required=True,
                        metavar='PATH',
                        type=str)

    parser.add_argument('-n', '--num_splits',
                        help='Number of files to output.',
                        required=True,
                        metavar='INT',
                        type=int)

    args = parser.parse_args()

    seqs = read_fasta(args.fasta, cut_header=True)

    all_ids = sorted(list(seqs.keys()))

    if args.num_splits <= 1:
        raise ValueError('Number of splits must be greater than 1.')

    split_size = len(all_ids) // args.num_splits
    remainder = len(all_ids) % args.num_splits
    splits = []
    start = 0
    
    for i in range(args.num_splits):
        end = start + split_size + (1 if i < remainder else 0)
        splits.append(all_ids[start:end])
        start = end

    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)

    for i, split in enumerate(splits):
        split_seqs = {seqid: seqs[seqid] for seqid in split}
        write_fasta(split_seqs, os.path.join(args.out_folder, 'split_%d.faa' % i))


if __name__ == '__main__':
    main()
