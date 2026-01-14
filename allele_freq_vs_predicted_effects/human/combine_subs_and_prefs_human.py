#!/usr/bin/python3

import sys
import os
from collections import defaultdict
import gzip

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from functions import possible_aa

# Hard-coded script specifically for Human SNVs, as these are in a specific format.
# Keep only most frequent substitution per protein position.

considered_proteins = set()

with open('/home6/gmdougla/projects/aa_distance/human_variants/nonsyn_snvs/gencode_w_snv.faa', 'r') as fasta_fh:
    for line in fasta_fh:
        if line.startswith('>'):
            considered_proteins.add(line.strip()[1:])

pos_per_protein = defaultdict(set)
subs_info = defaultdict(lambda: defaultdict(list))

snvfiles = os.listdir('/home6/gmdougla/projects/aa_distance/human_variants/nonsyn_snvs/by_chrom/')
for snvfile in snvfiles:
    with gzip.open('/home6/gmdougla/projects/aa_distance/human_variants/nonsyn_snvs/by_chrom/' + snvfile, 'rt') as subs_fh:
        subs_header = subs_fh.readline().strip().split('\t')
        for sub_line in subs_fh:
            sub_line = sub_line.strip().split('\t')
            protein_pos = int(sub_line[5]) - 1
            protein_id = sub_line[4]
            anc_aa = sub_line[6]
            derived_aa = sub_line[7]
            if anc_aa not in possible_aa or derived_aa not in possible_aa:
                continue
            if protein_id in considered_proteins:
                subs_info[protein_id][protein_pos].append(sub_line[6:])
                pos_per_protein[protein_id].add(protein_pos)

for protein_id in subs_info.keys():
    for protein_pos in subs_info[protein_id].keys():
        raw_sub_infos = subs_info[protein_id][protein_pos]
        max_count = 0
        i_to_keep = None
        for i in range(len(raw_sub_infos)):
            count = int(raw_sub_infos[i][2])
            if count > max_count:
                max_count = count
                i_to_keep = i
        subs_info[protein_id][protein_pos] = [raw_sub_infos[i_to_keep]]

pref_master_folder = '/home6/gmdougla/projects/aa_distance/human_variants/nonsyn_snvs/prefs'
pref_folders = sorted([f for f in os.listdir(pref_master_folder) if os.path.isdir(os.path.join(pref_master_folder, f))])

if len(pref_folders) == 0:
    sys.exit('Error: No preference folders found.')

print('protein_id\tprotein_pos\t' + '\t'.join(subs_header[6:]) + '\t' + '\t'.join(pref_folders))
for protein_id in sorted(list(considered_proteins)):
    prefs_per_protein = defaultdict(dict)
    for pref_folder in pref_folders:
        pref_file = os.path.join(pref_master_folder, pref_folder, protein_id + '.csv')
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
                    if exp_ref_aa not in possible_aa or derived_AA not in possible_aa:
                        outline.append('NA')
                    elif 'RvC' not in pref_folder and pref_folder != 'rasp' and prefs_per_protein[pref_folder][pos][exp_ref_aa] != '':
                        print(pref_folder)
                        print(pos)
                        print(protein_id)
                        print(prefs_per_protein[pref_folder][pos][exp_ref_aa])
                        print(prefs_per_protein[pref_folder][pos])
                        sys.exit('Error: Expected reference amino acid is not empty, which it should be (if it was filled with NA previously)')
                    else:
                        outline.append(prefs_per_protein[pref_folder][pos][derived_AA])
                else:
                    sys.exit('Error: No preference found for ' + protein_id + ' at position ' + str(pos) + ' in ' + pref_folder)
            print('\t'.join(outline))
