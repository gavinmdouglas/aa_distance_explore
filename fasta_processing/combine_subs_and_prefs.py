#!/usr/bin/python3

import argparse
import sys
import os
from collections import defaultdict
import gzip

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import read_fasta, possible_aa


def main():

    parser = argparse.ArgumentParser(
        description='''
        Combine table of segregating AA substitutions with information on preferences of that substitution at that site.
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-s', '--subs',
                        metavar='PATH',
                        type=str,
                        help='Table of subs produced by "segregating_subs.py" (gzipped).',
                        required=True)

    parser.add_argument('-f', '--folder',
                        metavar='PATH',
                        type=str,
                        help='Path to folder containing preference categories. Expecting at least one subfolder (named as the preference). Should be one table per protein in each subfolder.',
                        required=True)

    parser.add_argument('-o', '--output',
                        metavar='OUTPUT',
                        type=str,
                        help='Path to combined output table.',
                        required=True)

    args = parser.parse_args()

    unique_proteins = set()
    pos_per_protein = defaultdict(set)
    subs_info = defaultdict(lambda: defaultdict(list))
    with gzip.open(args.subs, 'rt') as subs_fh:
        subs_header = subs_fh.readline().strip().split('\t')
        for sub_line in subs_fh:
            sub_line = sub_line.strip().split('\t')
            pos = int(sub_line[1])
            subs_info[sub_line[0]][pos].append(sub_line[2:])
            unique_proteins.add(sub_line[0])
            pos_per_protein[sub_line[0]].add(pos)

    pref_folders = sorted([f for f in os.listdir(args.folder) if os.path.isdir(os.path.join(args.folder, f))])

    if len(pref_folders) == 0:
        sys.exit('Error: No preference folders found.')

    with open(args.output, 'w') as out_fh:
        out_fh.write('\t'.join(subs_header) + '\t' + '\t'.join(pref_folders) + '\n')

        for protein_id in sorted(list(unique_proteins)):
            prefs_per_protein = defaultdict(dict)
            for pref_folder in pref_folders:
                pref_file = os.path.join(args.folder, pref_folder, protein_id + '.csv')
                if not os.path.isfile(pref_file):
                    sys.exit('Error: No preference file found for ' + protein_id + ' in ' + pref_folder)
                with open(pref_file, 'r') as pref_fh:
                    pref_header = pref_fh.readline().strip().split(',')
                    for pref_line in pref_fh:
                        pref_line = pref_line.strip().split(',')
                        pos = int(pref_line[0]) - 1
                        prefs_per_protein[pref_folder][pos] = {}
                        for i in range(1, len(pref_line)):
                            prefs_per_protein[pref_folder][pos][pref_header[i]] = pref_line[i]
            
            for pos in sorted(list(pos_per_protein[protein_id])):
                raw_sub_infos = subs_info[protein_id][pos]
                for raw_sub_info in raw_sub_infos:
                    outline = [protein_id, str(pos)] + raw_sub_info
                    derived_AA = raw_sub_info[1]
                    exp_ref_aa = raw_sub_info[0]
                    for pref_folder in pref_folders:
                        if pos in prefs_per_protein[pref_folder]:
                            # Make sure exp_ref_aa value is empty.
                            if not 'RvC' in pref_folder and prefs_per_protein[pref_folder][pos][exp_ref_aa] != '':
                                print(pref_folder)
                                print(pos)
                                print(protein_id)
                                print(prefs_per_protein[pref_folder][pos][exp_ref_aa])
                                sys.exit('Error: Expected reference amino acid is not empty, which it should be (if it was filled with NA previously)')
                            outline.append(prefs_per_protein[pref_folder][pos][derived_AA])
                        else:
                            sys.exit('Error: No preference found for ' + protein_id + ' at position ' + str(pos) + ' in ' + pref_folder)
                    out_fh.write('\t'.join(outline) + '\n')


if __name__ == '__main__':
    main()
