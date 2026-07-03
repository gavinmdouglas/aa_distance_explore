#!/usr/bin/python3

import sys
import os
from collections import defaultdict
import gzip

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from functions import possible_aa, classify_Ts_or_Tv

# Combine table of segregating AA substitutions with information on preferences of that substitution at that site.
# Also parse DSSP values and determine whether each site is buried or exposed.
unique_uniprot = set()
pos_per_protein = defaultdict(set)
subs_info = defaultdict(lambda: defaultdict(list))
uniprot_to_blattner = {}

# Blattner_num    uniprot_id      pos     consensus_codon consensus_aa    alt_codon    alt_aa  alt_count       consensus_count AF
# b0003   P00547  6       CCG     P       ACG     T       2       2536    0.000788022064617809
# b0003   P00547  7       GCU     A       ACU     T       1       2617    0.000381970970206264
subs_file = "/home6/gmdougla/projects/aa_distance/ecoli_variants/parsed_allelome_variants.tsv.gz"
with gzip.open(subs_file, 'rt') as subs_fh:
    subs_header = subs_fh.readline().strip().split('\t')[2:]
    for sub_line in subs_fh:
        sub_line = sub_line.strip().split('\t')
        blattner_id = sub_line[0]
        uniprot_id = sub_line[1]
        uniprot_to_blattner[uniprot_id] = blattner_id
        pos = int(sub_line[2])
        if pos in subs_info[uniprot_id]:
            sys.exit('Error: Duplicate position found for ' + uniprot_id + ' at position ' + str(pos) + ' in substitutions file.')
        subs_info[uniprot_id][pos] = sub_line[3:]
        unique_uniprot.add(uniprot_id)
        pos_per_protein[uniprot_id].add(pos)

# RSA < 0.2 set to buried
# RSA >= 0.2 set to exposed
rsa_folder = "/home6/gmdougla/projects/aa_distance/ecoli_variants/UP000000625_83333_ECOLI_v4_dssp_rsa"
rsa_files = [f for f in os.listdir(rsa_folder) if os.path.isfile(os.path.join(rsa_folder, f))]

for rsa_file in rsa_files:
    rsa_uniprot = rsa_file.split('-')[1]
    if rsa_uniprot not in unique_uniprot:
        continue

    aa_matches = 0
    total_pos = 0

    # Expected format:
    # pdb_position,chain,structure,pdb_aa,rsa
    # 1,A,,M,0.9107142857142857
    with gzip.open(os.path.join(rsa_folder, rsa_file), 'rt') as rsa_fh:
        rsa_header = rsa_fh.readline().strip().split(',')
        if rsa_header != ['pdb_position', 'chain', 'structure', 'pdb_aa', 'rsa']:
            sys.exit('Error: RSA file header does not match expected format. Found: ' + ','.join(rsa_header) + ' Expected: pdb_position,chain,structure,pdb_aa,rsa')
        for rsa_line in rsa_fh:
            rsa_line = rsa_line.strip().split(',')
            
            pos = int(rsa_line[0]) - 1
            if pos in pos_per_protein[rsa_uniprot]:
                total_pos += 1
                rsa_value = float(rsa_line[4])

                # Make sure amino acid matches.
                if rsa_line[3] == subs_info[rsa_uniprot][pos][1]:
                    aa_matches += 1

                if len(subs_info[rsa_uniprot][pos]) > 7:
                        sys.exit('Error: Duplicate position found for ' + rsa_uniprot + ' at position ' + str(pos) + ' in substitutions file.')

                subs_info[rsa_uniprot][pos].append(str(rsa_value))

                if rsa_value < 0.2:
                    subs_info[rsa_uniprot][pos].append('buried')
                else:
                    subs_info[rsa_uniprot][pos].append('exposed')

    mismatches = total_pos - aa_matches
    if mismatches > 1:
        sys.exit("More than 1 AA mismatch found for " + rsa_uniprot + " between substitutions file and RSA file. Found " + str(mismatches) + " mismatches out of " + str(total_pos) + " positions.")

pref_folder_path = "/home6/gmdougla/projects/aa_distance/ecoli_variants/prefs"
pref_folders = sorted([f for f in os.listdir(pref_folder_path) if os.path.isdir(os.path.join(pref_folder_path, f))])

if len(pref_folders) == 0:
    sys.exit('Error: No preference folders found.')

output_tab = "/home6/gmdougla/projects/aa_distance/ecoli_variants/parsed_allelome_variants_w_prefs_dssp.tsv.gz"
with gzip.open(output_tab, 'wt') as out_fh:
    subs_header.insert(5, 'Ts_or_Tv')
    out_fh.write('\t'.join(['protein']) + '\t' + '\t'.join(subs_header) + '\t' + 'RSA' + '\t' + 'RSA_group' + '\t' + '\t'.join(pref_folders) + '\n')

    for protein_id in sorted(list(unique_uniprot)):
        prefs_per_protein = defaultdict(dict)
        mismatch_counts = defaultdict(int)
        for pref_folder in pref_folders:

            if pref_folder == "vespag" or pref_folder == "vespag_sum_scaled":
                protein_id_match = uniprot_to_blattner[protein_id] + '_' + protein_id
            else:
                protein_id_match = protein_id

            pref_file = os.path.join(pref_folder_path, pref_folder, protein_id_match + '.csv')
            if not os.path.isfile(pref_file):
                sys.exit('Error: No preference file found for ' + protein_id_match + ' in ' + pref_folder)
            with open(pref_file, 'r') as pref_fh:
                pref_header = pref_fh.readline().strip().split(',')
                for pref_line in pref_fh:
                    pref_line = pref_line.strip().split(',')
                    pos = int(pref_line[0]) - 1
                    prefs_per_protein[pref_folder][pos] = {}
                    for i in range(1, len(pref_line)):
                        prefs_per_protein[pref_folder][pos][pref_header[i]] = pref_line[i]
        
        for pos in sorted(list(pos_per_protein[protein_id])):
            raw_sub_info = subs_info[protein_id][pos]

            if len(raw_sub_info) != 9 or raw_sub_info[-1] not in ['buried', 'exposed']:
                print(raw_sub_info)
                sys.exit('Error: Expected 9 columns of raw substitution info followed by RSA group, but found ' + str(len(raw_sub_info)) + ' columns for ' + protein_id + ' at position ' + str(pos))
            
            outline = [protein_id, str(pos)] + raw_sub_info
            anc_codon = outline[2]
            des_codon = outline[4]
            transition_type = classify_Ts_or_Tv(anc_codon, des_codon)
            outline.insert(6, transition_type)

            derived_AA = raw_sub_info[3]
            exp_ref_aa = raw_sub_info[1]
            for pref_folder in pref_folders:
                if pos in prefs_per_protein[pref_folder]:
                    # Make sure exp_ref_aa value is empty.
                    if prefs_per_protein[pref_folder][pos][exp_ref_aa] != '':
                        mismatch_counts[pref_folder] += 1

                    if prefs_per_protein[pref_folder][pos][derived_AA] == '':
                         outline.append('NA')
                    else:
                        outline.append(prefs_per_protein[pref_folder][pos][derived_AA])
                else:
                    sys.exit('Error: No preference found for ' + protein_id + ' at position ' + str(pos) + ' in ' + pref_folder)
                
            out_fh.write('\t'.join(outline) + '\n')

        for pref_folder in pref_folders:
            if mismatch_counts[pref_folder] > 1:
                sys.exit("More than one mismatch found for " + pref_folder + " in " + protein_id + ". Found " + str(mismatch_counts[pref_folder]) + " mismatches.")
