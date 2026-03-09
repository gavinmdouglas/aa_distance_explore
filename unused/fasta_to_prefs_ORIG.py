#!/usr/bin/env python

import argparse
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import (read_rad_vs_cons_effects,
                       read_sim_matrix,
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

    parser.add_argument('-r', '--rc_folder',
                        help='''
                            Folder containing gzipped tables breaking down radical vs. conservative substitutions based on amino acid pairs.
                            Generally symmetrical, but could in principle by assymetrical (if one could justify it!).
                            "aa1" is ancestral.
                            File must end in ".tsv.gz".
                            Also, the category name will be the file prefix until the first ".".
                        ''',
                        required=True,
                        metavar='FOLDER',
                        type=str)

    parser.add_argument('-s', '--sim_folder',
                        help='''
                            Folder containing gzipped tables breaking down similarities of all pairwise amino acids.
                            Usually symmetrical, but not necessarily.
                            "aa1" is ancestral.
                            File must end in ".tsv.gz".
                            Also, the category name will be the file prefix until the first ".".
                        ''',
                        required=True,
                        metavar='FOLDER',
                        type=str)
    
    # Float between >0 and <0.5 to fill in as similarity for radical vs conservative effects.
    parser.add_argument('-p', '--pref_radical_fill',
                        help='''
                            Fill in similarity for amino acid pairs where substitutions are considered radical (must be > 0 and < 0.5).
                            The corresponding conservative value is simply 1 - this value.
                            Default is 0.05.
                        ''',
                        required=False,
                        metavar='FLOAT',
                        type=float,
                        default=0.05)

    parser.add_argument('-i', '--identity_fill',
                        help='''
                            Fill in similarity for amino acid pairs where the identity is the same.
                            Options are "max", "median", "NA", or "one".
                            Default is "max".
                        ''',
                        required=False,
                        metavar='STR',
                        type=str,
                        default='max')

    args = parser.parse_args()

    if not 0 < args.pref_radical_fill < 0.5:
        sys.exit('Error: Fill in value for radical preferences must be > 0 and < 0.5.')

    if not os.path.exists(args.fasta):
        sys.exit('Error: FASTA file does not exist.')

    # Make sure that input folders exist and there are matching files within.
    if not os.path.exists(args.rc_folder):
        sys.exit('Error: Radical vs conservative effect folder does not exist.')
    elif not os.path.exists(args.sim_folder):
        sys.exit('Error: Similarity effect folder does not exist.')
    
    rc_files = [f for f in os.listdir(args.rc_folder) if f.endswith('.tsv.gz')]
    sim_files = [f for f in os.listdir(args.sim_folder) if f.endswith('.tsv.gz')]
    if len(rc_files) == 0:
        sys.exit('Error: No radical vs conservative effect files found.')
    elif len(sim_files) == 0:
        sys.exit('Error: No similarity effect files found.')
    
    seqs = read_fasta(args.fasta, cut_header=True)

    if len(seqs) == 0:
        sys.exit('Error: No sequences found in FASTA.')

    sub_maps = []
    tab_ids = []
    for rc_file in rc_files:
        r_vs_c_dict = read_rad_vs_cons_effects(sub_file=os.path.join(args.rc_folder, rc_file),
                                               swap_rad_cons=True,
                                               rad_pref=args.pref_radical_fill)
        
        dict_check = continuous_prefs_dict_sanity(r_vs_c_dict)
        if not dict_check:
            sys.exit('Error: Radical vs conservative effect table is not formatted correctly: ' + os.path.join(args.rc_folder, rc_file))

        sub_maps.append(r_vs_c_dict)
        tab_ids.append(rc_file.split('.')[0])

    for sim_file in sim_files:
        sim_dict = read_sim_matrix(sub_sim_file=os.path.join(args.sim_folder, sim_file),
                                   identity_set=args.identity_fill)

        dict_check = continuous_prefs_dict_sanity(sim_dict)
        if not dict_check:
            sys.exit('Error: Similarity table is not formatted correctly: ' + os.path.join(args.sim_folder, sim_file))

        sub_maps.append(sim_dict)
        tab_ids.append(sim_file.split('.')[0])

    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)

    for effect_type in tab_ids:
        if not os.path.exists(args.out_folder + '/' + effect_type):
            os.makedirs(args.out_folder + '/' + effect_type)

    for seq_id, seq in seqs.items():
        for sub_map, effect_type in zip(sub_maps, tab_ids):
            outfile = args.out_folder + '/' + effect_type + '/' + seq_id + '.csv'
            prefs = seq_sub_prob_prefs(seq=seq, effect_tab=sub_map, outfile=outfile)


if __name__ == '__main__':
    main()
