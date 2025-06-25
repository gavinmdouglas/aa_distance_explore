import os
import sys
from collections import defaultdict
import pandas as pd
import numpy as np

possible_aa = set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])

# Suppress future warnings.
#import warnings
#warnings.simplefilter(action='ignore', category=FutureWarning)

# Parse DMS enrichment values to amino acid preferences for input to phydms.
# Read through the DMS CSVs to get the tested positions.
# NOTE: This was the original version of the script, which prepared the DMS data assuming it was all log-fold change data compared to the reference.
# I later realized there is a mixture of different DMS data-types, so a different approach was used.
dms_folder = "/Users/gavin/Drive/research/aa_distance/data/proteinGym/DMS_ProteinGym_substitutions/"
dms_files = os.listdir(dms_folder)

files_to_ignore = set(['POLG_HCVJF_Qi_2014.csv', 'POLG_PESV_Tsuboyama_2023_2MXD.csv', 'POLG_CXB3N_Mattenberger_2021.csv', 'OXDA_RHOTO_Vanella_2023_expression.csv',
                       'P53_HUMAN_Giacomelli_2018_Null_Etoposide.csv', 'P53_HUMAN_Giacomelli_2018_Null_Nutlin.csv', 'P53_HUMAN_Giacomelli_2018_WT_Nutlin.csv'])

genes = set()

empty_sub_scores = {}
sorted_aa = sorted(list(possible_aa))
for i in range(len(sorted_aa) - 1):
    for j in range(i + 1, len(sorted_aa)):
        aa1 = sorted_aa[i]
        aa2 = sorted_aa[j]
        empty_sub_scores[aa1 + aa2] = []

mean_per_gene = empty_sub_scores.copy()

for dms_file in dms_files:
    if dms_file in files_to_ignore:
        continue

    site_wt = dict()
    site_to_effect = dict()

    # Make sure that all sites (except position 1) were tested for all possible derived amino acids.
    with open(os.path.join(dms_folder, dms_file), "r") as dms_fh:
        header = dms_fh.readline()
        if header.strip() != "mutant,mutated_sequence,DMS_score,DMS_score_bin" and header.strip() != "mutant,mutated_sequence,DMS_score,DMS_bin_score,DMS_score_bin":
            sys.exit("Error: Unexpected header in DMS file: " + dms_file + " " + header)

        min_dms_score = 1000000
        max_dms_score = -1000000

        seq_len = None

        raw_scores = []

        for line in dms_fh:
            line = line.strip().split(",")
            mutant = line[0]
            mutated_seq = line[1]
            dms_score = float(line[2])
            dms_score_bin = line[3]

            if seq_len and len(mutated_seq) != seq_len:
                sys.exit("Error: Different sequence lengths in DMS file: " + dms_file)
            else:
                seq_len = len(mutated_seq)

            if ":" in mutant:
                continue

            wt_aa = mutant[0]
            pos = int(mutant[1:-1])
            mut_aa = mutant[-1]

            if wt_aa not in possible_aa or mut_aa not in possible_aa:
                sys.exit("Error: Unexpected amino acid in DMS file: " + dms_file + " for line " + mutant)

            if pos == 1:
                continue

            if dms_score < min_dms_score:
                min_dms_score = dms_score
            if dms_score > max_dms_score:
                max_dms_score = dms_score

            if pos not in site_wt.keys():
                site_wt[pos] = wt_aa
                site_to_effect[pos] = dict()
            else:
                if site_wt[pos] != wt_aa:
                    sys.exit("Error: Different reference amino acid for same position in " + dms_file)
            
            if pos not in site_to_effect.keys():
                site_to_effect[pos] = dict()
            
            if mut_aa not in site_to_effect[pos].keys():
                site_to_effect[pos][mut_aa] = dms_score
                raw_scores.append(dms_score)
            else:
                sys.exit("Error: Duplicate DMS test for " + dms_file + " at position " + str(pos) + " for amino acid " + mut_aa)

        num_sites_tested = 0
        site_num_aa_tested = []
        for i in range(2, seq_len + 1):
            if i not in site_wt.keys():
                break
            num_sites_tested += 1

            site_num_aa_tested.append(len(site_to_effect[i].keys()))

        prop_sites_tested = num_sites_tested / (seq_len - 1)

        if num_sites_tested >= 20 and prop_sites_tested >= 0.95 and np.mean(site_num_aa_tested) >= 15:
            dms_file_split = dms_file.split("_")
            if dms_file_split[0] in genes:
                sys.exit("Error: Duplicate gene in DMS file: " + dms_file_split[0])
            print(dms_file)
