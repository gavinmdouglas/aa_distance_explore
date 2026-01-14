#!/usr/bin/env python

import argparse
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import (read_sim_matrix,
                       read_fasta,
                       seq_sub_prob_prefs,
                       continuous_prefs_dict_sanity)


def main():
    parser = argparse.ArgumentParser(

        description=
        '''
        Parse a FASTA containing amino acid sequences, and generate an AA preference table for each sequence.
        Do so based on all radical vs conservative effects in specified tables, as well as any similarity metrics in separate specified tables.
        ''',

        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--fasta',
                        help='FASTA of protein sequences (*as amino acids*).',
                        required=True,
                        metavar='FILE',
                        type=str)

    parser.add_argument('-o', '--out_folder',
                        help='Output folder for outfiles (a subfolder will be created for each effect type).',
                        required=True,
                        metavar='PATH',
                        type=str)

    parser.add_argument('-s', '--sim_file',
                        help='''
                            File containing gzipped table breaking down similarities of all pairwise amino acids.
                            Usually symmetrical, but not necessarily.
                            "aa1" is ancestral.
                            Also, the category name will be the file prefix until the first ".".
                        ''',
                        required=True,
                        metavar='FOLDER',
                        type=str)

    parser.add_argument('-i', '--identity_fill',
                        help='''
                            Fill in similarity for amino acid pairs where the identity is the same.
                            Options are "max", "median", "NA", or "one".
                            Default is "one".
                        ''',
                        required=False,
                        metavar='STR',
                        type=str,
                        default='one')

    args = parser.parse_args()

    if not os.path.exists(args.fasta):
        sys.exit('Error: FASTA file does not exist.')
    
    seqs = read_fasta(args.fasta, cut_header=True)

    if len(seqs) == 0:
        sys.exit('Error: No sequences found in FASTA.')
    
    sim_dict = read_sim_matrix(sub_sim_file=args.sim_file,
                               identity_set=args.identity_fill)

    dict_check = continuous_prefs_dict_sanity(sim_dict)
    if not dict_check:
        sys.exit('Error: Similarity table is not formatted correctly: ' + args.sim_file)

    pref_id = os.path.basename(args.sim_file).split('.')[0]

    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)

    outfolder = args.out_folder + '/' + pref_id
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    for seq_id, seq in seqs.items():
        outfile = outfolder + '/' + seq_id + '.csv'
        seq_sub_prob_prefs(seq=seq, effect_tab=sim_dict, outfile=outfile)

if __name__ == '__main__':
    main()
