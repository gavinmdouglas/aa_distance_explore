import os
import sys
from collections import defaultdict
import numpy as np
import csv
from scipy import stats

possible_aa = set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])

dms_ids = []
with open('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_score_type.csv', 'r') as f:
    csv_reader = csv.DictReader(f)
    for line in csv_reader:
        if line['DMS_id'] == 'NCAP_I34A1_Doud_2015':
            continue
        else:
            dms_ids.append(line['DMS_id'])


dms_folder = "/Users/gavin/Drive/research/aa_distance/data/proteinGym/DMS_ProteinGym_substitutions/"
mean_standard_score = defaultdict(list)
median_robust_score = defaultdict(list)

for dms_id in dms_ids:
    site_wt = dict()
    site_to_effect = dict()

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
            else:
                sys.exit("Error: Duplicate DMS test for " + dms_file + " at position " + str(pos) + " for amino acid " + mut_aa)

        gene_standard_score = defaultdict(list)
        gene_robust_score = defaultdict(list)
        all_vals = []
        for site in site_to_effect.keys():
            site_scores = []
            for aa2 in site_to_effect[site].keys():
                site_scores.append(site_to_effect[site][aa2])

            if len(site_scores) < 15:
                continue
            else:
                all_vals = all_vals + site_scores
            
        mean_score = np.mean(all_vals)
        sd_score = np.std(all_vals)

        median_score = np.median(all_vals)
        mad_score = stats.median_abs_deviation(all_vals)

        if mad_score == 0:
            print("Warning: Median absolute deviation is zero for DMS file: " + dms_file + ". Will skip this gene.", file=sys.stderr)

        for site in site_to_effect.keys():
            ref_aa = site_wt[site]
            for aa2 in site_to_effect[site].keys():
                if aa2 == ref_aa:
                    sys.exit("Error: Unexpected reference amino acid in DMS file: " + dms_file + " for line " + str(site) + " " + ref_aa)
                else:
                    aa_combo = tuple(sorted([ref_aa, aa2]))
                    standard_score = (site_to_effect[site][aa2] - mean_score) / sd_score
                    gene_standard_score[aa_combo].append(standard_score)

                    if mad_score > 0:
                        robust_score = (site_to_effect[site][aa2] - median_score) / (1.4826 * mad_score)
                        gene_robust_score[aa_combo].append(robust_score)

        for aa_combo in gene_standard_score.keys():
            mean_standard_score[aa_combo].append(np.mean(gene_standard_score[aa_combo]))

            if mad_score > 0:
                median_robust_score[aa_combo].append(np.median(gene_robust_score[aa_combo]))

print('aa1\taa2\tmean_standard_score\tmedian_robust_score\tnum_genes_w_obs')
for aa_combo in sorted(list(mean_standard_score.keys())):
    print(aa_combo[0], aa_combo[1], np.mean(mean_standard_score[aa_combo]), np.median(median_robust_score[aa_combo]), len(mean_standard_score[aa_combo]), sep="\t")
