#!/usr/bin/python3

import sys
import os
from collections import defaultdict
import gzip

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from functions import possible_aa

# Combine table of segregating AA substitutions with information on preferences of that substitution at that site.
# Unlike E. coli case, just start from previously generated table and add in thermoMPNN and DSSP classifications into buried and exposed.

# First get id map (IDs used in AlphaFold models vs ENSEMBL genes)
id_mapfile = "/home6/gmdougla/projects/aa_distance/human_variants/gencode/rasp_gencode_map.tsv"
alpha_to_ensembl = {}
ensembl_to_alpha = {}
past_ids = set()
with open(id_mapfile, 'r') as id_map_fh:
    header = id_map_fh.readline()
    for id_line in id_map_fh:
        id_line = id_line.strip().split('\t')
        alpha_id = id_line[0].split(';')[0]
        ensembl_id = id_line[1].split('|')[0]
        if alpha_id in past_ids or ensembl_id in past_ids:
            sys.exit("ID already present")
        else:
            past_ids.add(alpha_id)
            past_ids.add(ensembl_id)
        alpha_to_ensembl[alpha_id] = ensembl_id
        ensembl_to_alpha[ensembl_id] = alpha_id

subs_info = defaultdict(lambda: defaultdict(list))

# protein_id      protein_pos     ref_base        alt_base        ref_aa  alt_aa  variant_count   sample_count    Ts_or_Tv        CpG_or_not      rasp    vespag  vespag_sum_scaled
# ENSP00000000233.5       33      A       G       I       V       136     1461854 transition      Not_CpG 0.8672465405609342      0.2100  0.2944
# ENSP00000000233.5       150     A       G       S       G       69      1461680 transition      Not_CpG 0.4735385367838163      0.2535  0.0898
subs_file = "/home6/gmdougla/projects/aa_distance/human_variants/nonsyn_snvs/human_snvs_w_prefs.tsv.gz"
with gzip.open(subs_file, 'rt') as subs_fh:
    subs_header = subs_fh.readline().strip().split('\t')[2:]
    for sub_line in subs_fh:
        sub_line = sub_line.strip().split('\t')
        ensembl_id = sub_line[0]
        pos = int(sub_line[1])
        if pos in subs_info[ensembl_id]:
            sys.exit('Error: Duplicate position found for ' + ensembl_id + ' at position ' + str(pos) + ' in substitutions file.')
        subs_info[ensembl_id][pos] = sub_line[2:]

# Read through RSA files.
rsa_folder = "/home6/gmdougla/projects/aa_distance/human_variants/UP000005640_9606_HUMAN_v4_dssp_rsa"
rsa_files = [f for f in os.listdir(rsa_folder) if os.path.isfile(os.path.join(rsa_folder, f))]

# First loop through all and figure out all RSA files per protein.
protein_to_rsa_files = defaultdict(dict)
for rsa_file in rsa_files:
    rsa_split = rsa_file.split('-')
    rsa_uniprot = rsa_split[1]
    fragment_num = int(rsa_split[2][1:])
    protein_to_rsa_files[rsa_uniprot][fragment_num] = rsa_file

# Fragments are split across different coordinates, which I grepped from the PDBs into UP000005640_9606_HUMAN_v4_pdb_DBREF_grep.txt
# AF-A0A024R1R8-F1-model_v4.pdb:DBREF  XXXX A    1    64  UNP    A0A024R1R8 A0A024R1R8_HUMAN     1     64
# AF-A0A024RBG1-F1-model_v4.pdb:DBREF  XXXX A    1   181  UNP    A0A024RBG1 NUD4B_HUMAN      1    181
dbref_info_file = "/home6/gmdougla/projects/aa_distance/human_variants/UP000005640_9606_HUMAN_v4_pdb_DBREF_grep.txt"
fragment_start = {}
with open(dbref_info_file, 'r') as dbref_fh:
    for line in dbref_fh:
        line = line.strip().split()
        pdb_file = line[0].split(':')[0]
        pdb_file_split = pdb_file.split('-')
        uniprot_id = pdb_file_split[1]
        if uniprot_id not in protein_to_rsa_files:
            continue
        fragment_num = int(pdb_file_split[2][1:])
        fragment_start[(uniprot_id, fragment_num)] = int(line[8])

# Then loop through proteins.
pos_to_rsa_group = dict()
for rsa_uniprot in protein_to_rsa_files.keys():
    max_fragment = max(protein_to_rsa_files[rsa_uniprot].keys())

    if rsa_uniprot not in alpha_to_ensembl:
        continue
    ensembl_protein_id = alpha_to_ensembl[rsa_uniprot]
    pos_to_rsa_group[ensembl_protein_id] = dict()

    for i in range(1, max_fragment + 1):
        if i not in protein_to_rsa_files[rsa_uniprot].keys():
            sys.exit("Missing RSA file for " + rsa_uniprot + " fragment " + str(i))

        fragment_start_pos = fragment_start[(rsa_uniprot, i)] - 1

        # Expected format:
        # pdb_position,chain,structure,pdb_aa,rsa
        # 1,A,,M,0.9107142857142857
        rsa_filepath = os.path.join(rsa_folder, protein_to_rsa_files[rsa_uniprot][i])
        with gzip.open(rsa_filepath, 'rt') as rsa_fh:
            rsa_header = rsa_fh.readline().strip().split(',')
            if rsa_header != ['pdb_position', 'chain', 'structure', 'pdb_aa', 'rsa']:
                sys.exit('Error: RSA file header does not match expected format. Found: ' + ','.join(rsa_header) + ' Expected: pdb_position,chain,structure,pdb_aa,rsa')
            for rsa_line in rsa_fh:
                rsa_line = rsa_line.strip().split(',')
                
                # Note that last pos added to position as proteins split into multiple fragments always start from 1, even though they can be in middle.
                pos = int(rsa_line[0]) + fragment_start_pos
                rsa_value = float(rsa_line[4])
                pdb_aa = rsa_line[3]

                if pos not in subs_info[ensembl_protein_id]:
                    continue

                if pdb_aa != subs_info[ensembl_protein_id][pos][2]:
                    sys.exit('Mismatch of AA for ' + ensembl_protein_id + ' at position ' + str(pos) + '. Expected: ' + subs_info[ensembl_protein_id][pos][2] + ' Found: ' + pdb_aa)

                if pos in pos_to_rsa_group[ensembl_protein_id]:
                    pos_to_rsa_group[ensembl_protein_id][pos].append(rsa_value)
                else:
                    pos_to_rsa_group[ensembl_protein_id][pos] = [rsa_value]

# Then read in ThermoMPNN.
pref_folder_path = "/home6/gmdougla/projects/aa_distance/human_variants/nonsyn_snvs/prefs/thermoMPNN"

# Just like RSA files, need to figure out which files are fragments of same protein first.
protein_to_pref_files = defaultdict(dict)
for pref_file in os.listdir(pref_folder_path):
    if not os.path.isfile(os.path.join(pref_folder_path, pref_file)):
        continue
    pref_split = pref_file.split('-')
    pref_uniprot = pref_split[0]
    fragment_num = int(pref_split[1][1:].replace('.csv', ''))
    protein_to_pref_files[pref_uniprot][fragment_num] = pref_file


thermompnn_values = defaultdict(dict)
for rsa_uniprot in protein_to_pref_files.keys():
    max_fragment = max(protein_to_pref_files[rsa_uniprot].keys())

    if rsa_uniprot not in alpha_to_ensembl:
        continue
    ensembl_protein_id = alpha_to_ensembl[rsa_uniprot]

    for i in range(1, max_fragment + 1):
        if i not in protein_to_pref_files[rsa_uniprot].keys():
            sys.exit("Missing preference file for " + rsa_uniprot + " fragment " + str(i))

        fragment_start_pos = fragment_start[(rsa_uniprot, i)] - 1

        pref_filepath = os.path.join(pref_folder_path, protein_to_pref_files[rsa_uniprot][i])
        with open(pref_filepath, 'r') as pref_fh:
            pref_header = pref_fh.readline().strip().split(',')

            col_to_i = {pref_header[i]: i for i in range(len(pref_header))}
            for pref_line in pref_fh:
                pref_line = pref_line.strip().split(',')
                pos = int(pref_line[0]) + fragment_start_pos
                
                if pos not in subs_info[ensembl_protein_id]:
                    continue

                ref_aa = subs_info[ensembl_protein_id][pos][2]
                derived_aa = subs_info[ensembl_protein_id][pos][3]

                if ref_aa not in possible_aa or derived_aa not in possible_aa:
                    continue

                ref_aa_i = col_to_i[ref_aa]
                derived_aa_i = col_to_i[derived_aa]

                if pref_line[ref_aa_i] != '':
                    sys.exit('Expected ref AA to have empty pref setting.')

                alt_aa_pref = float(pref_line[derived_aa_i])

                if pos in thermompnn_values[ensembl_protein_id]:
                    thermompnn_values[ensembl_protein_id][pos].append(alt_aa_pref)
                else:
                    thermompnn_values[ensembl_protein_id][pos] = [alt_aa_pref]


output_tab = "/home6/gmdougla/projects/aa_distance/human_variants/nonsyn_snvs/human_snvs_w_prefs_dssp.tsv.gz"
with gzip.open(output_tab, 'wt') as out_fh:
    out_fh.write('\t'.join(['ensembl_protein', 'uniprot_id', 'pos']) + '\t' + '\t'.join(subs_header) + '\t' + 'thermoMPNN' + '\t' + 'RSA_group' + '\n')

    for ensembl_protein_id in sorted(list(subs_info.keys())):
        rsa_uniprot = ensembl_to_alpha[ensembl_protein_id]

        for pos in sorted(list(subs_info[ensembl_protein_id].keys())):

            ref_aa = subs_info[ensembl_protein_id][pos][2]
            derived_aa = subs_info[ensembl_protein_id][pos][3]
            if ref_aa not in possible_aa or derived_aa not in possible_aa:
                continue

            if pos not in thermompnn_values[ensembl_protein_id]:
                sys.exit('Pos ' + str(pos) + ' not found in thermompnn_values for ' + ensembl_protein_id)
            if pos not in pos_to_rsa_group[ensembl_protein_id]:
                sys.exit('Pos ' + str(pos) + ' not found in pos_to_rsa_group for ' + ensembl_protein_id)

            outline = [ensembl_protein_id, rsa_uniprot, str(pos)] + subs_info[ensembl_protein_id][pos]

            mean_thermompnn = sum(thermompnn_values[ensembl_protein_id][pos]) / len(thermompnn_values[ensembl_protein_id][pos])
            mean_rsa = sum(pos_to_rsa_group[ensembl_protein_id][pos]) / len(pos_to_rsa_group[ensembl_protein_id][pos])

            if mean_rsa < 0.2:
                rsa_group = 'buried'
            else:
                rsa_group = 'exposed'
            
            outline += [str(mean_thermompnn), rsa_group]
            out_fh.write('\t'.join(outline) + '\n')
