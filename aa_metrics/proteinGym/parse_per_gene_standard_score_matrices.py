import os
import sys
from collections import defaultdict
import pandas as pd
import numpy as np
import csv

possible_aa = set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
unique_aa = sorted(list(possible_aa))

dms_to_transform = dict()
with open('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_score_type_unambig.csv', 'r') as f:
    csv_reader = csv.DictReader(f)
    for line in csv_reader:
        dms_to_transform[line['DMS_id']] = line['dms_score_transform']

dms_folder = "/Users/gavin/Drive/research/aa_distance/data/proteinGym/DMS_ProteinGym_substitutions/"
dms_ids = sorted(list(dms_to_transform.keys()))

mean_per_gene = defaultdict(list)

for dms_id in dms_ids:
    site_wt = dict()
    site_to_effect = dict()
    all_vals = []
    
    # Make sure that all sites (except position 1) were tested for all possible derived amino acids.
    dms_file = os.path.join(dms_folder, dms_id + '.csv')
    with open(dms_file, "r") as dms_fh:
        header = dms_fh.readline()
        if header.strip() != "mutant,mutated_sequence,DMS_score,DMS_score_bin" and header.strip() != "mutant,mutated_sequence,DMS_score,DMS_bin_score,DMS_score_bin":
            sys.exit("Error: Unexpected header in DMS file: " + dms_file + " " + header)

        seq_len = None

        for line in dms_fh:
            line = line.strip().split(",")
            mutant = line[0]
            mutated_seq = line[1]
            raw_dms_score = float(line[2])
            dms_score_bin = line[3]

            if dms_to_transform[dms_id] == "ln":
                dms_score = raw_dms_score
            elif dms_to_transform[dms_id] == "log2":
                dms_score = np.log(2 ** raw_dms_score)
            elif dms_to_transform[dms_id] == "raw ratio":
                if raw_dms_score <= 0:
                    raw_dms_score = 1e-10
                dms_score = np.log(raw_dms_score)
            else:
                sys.exit("Error: Unexpected DMS score transformation type: " + dms_to_transform[dms_id])

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
                all_vals.append(dms_score)
            else:
                sys.exit("Error: Duplicate DMS test for " + dms_file + " at position " + str(pos) + " for amino acid " + mut_aa)

        mean_score = np.mean(all_vals)
        std_score = np.std(all_vals)

        gene_vals = defaultdict(list)
        for site in site_to_effect.keys():
            ref_aa = site_wt[site]
            for aa2 in site_to_effect[site].keys():
                aa_combo = (ref_aa, aa2)
                standardized_score = (site_to_effect[site][aa2] - mean_score) / std_score
                gene_vals[aa_combo].append(standardized_score)

        gene_out = '/Users/gavin/working_exchange_matrix/' + dms_id + '_per_gene_standard_scores_log.tsv'
        with open(gene_out, 'w') as gene_fh:
            gene_fh.write("aa_combo\tstandardized_score\n")
            for aa1 in unique_aa:
                for aa2 in unique_aa:
                    aa_combo = (aa1, aa2)
                    if aa1 == aa2:
                        continue
                    elif aa_combo not in gene_vals:
                        standardized_score = 'NA'
                    else:
                        standardized_score = np.mean(gene_vals[aa_combo])

                    gene_fh.write(f"{','.join(aa_combo)}\t{standardized_score}\n")