import os
import sys
from collections import defaultdict
import numpy as np
from scipy import stats
import csv

def compute_ranks_beta_cdf(ranks, alpha_opt=0.1848, beta_opt=0.3406):
    normalized_ranks = [(rank - 0.5) / (len(ranks) - 1) for rank in ranks]
    return stats.beta.cdf(normalized_ranks, alpha_opt, beta_opt)

if __name__ == "__main__":
    possible_aa = set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])

    dms_ids = []
    with open('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_score_type.csv', 'r') as f:
        csv_reader = csv.DictReader(f)
        for line in csv_reader:
            if line['DMS_id'] == 'NCAP_I34A1_Doud_2015':
                continue
            else:
                dms_ids.append(line['DMS_id'])

    mean_background_exposed_file = "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/stuctures/background_exposed_breakdown.tsv"
    mean_background_exposed = dict()
    with open(mean_background_exposed_file, 'r') as f:
        header = f.readline().strip().split("\t")
        if header != ['aa', 'yeast', 'ecoli', 'human', 'mean', 'median']:
            sys.exit("Error: Unexpected header in background exposed file: " + mean_background_exposed_file + " " + str(header))
        for line in f:
            line = line.strip().split("\t")
            aa = line[0]
            if aa not in possible_aa:
                sys.exit("Error: Unexpected amino acid in background exposed file: " + mean_background_exposed_file + " " + aa)
            mean_background_exposed[aa] = float(line[4])

    dms_folder = "/Users/gavin/Drive/research/aa_distance/data/proteinGym/DMS_ProteinGym_substitutions/"
    buried_gene_mean = defaultdict(list)
    exposed_gene_mean = defaultdict(list)
    buried_gene_median = defaultdict(list)
    exposed_gene_median = defaultdict(list)

    for dms_id in dms_ids:
        buried_raw = []
        exposed_raw = []
        buried_aa_combos = []
        exposed_aa_combos = []

        rsa_num_sites = 0
        buried_sites = set()
        exposed_sites = set()

        rsa_id = dms_id.split("_")[0] + '_' + dms_id.split("_")[1]
        rsa_file = os.path.join("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/stuctures/ProteinGym_AF2_structures_dssp_rsa", rsa_id + '.csv')
        with open(rsa_file, 'r') as rsa_fh:
            header = rsa_fh.readline().strip().split(",")
            if header != ['pdb_position', 'chain', 'structure', 'pdb_aa', 'rsa']:
                sys.exit("Error: Unexpected header in RSA file: " + rsa_file + " " + str(header))
            for line in rsa_fh:
                rsa_num_sites += 1
                line = line.strip().split(",")
                pos = int(line[0])
                rsa = float(line[4])

                if pos != rsa_num_sites:
                    sys.exit("Error: Position mismatch in RSA file: " + rsa_file + " at line " + str(pos) + " expected " + str(rsa_num_sites))

                if rsa < 0 or rsa > 10:
                    sys.exit("Error: Unexpected RSA value in RSA file: " + rsa_file + " at line " + str(pos) + " " + str(rsa))

                if rsa < 0.2:
                    buried_sites.add(pos)
                else:
                    exposed_sites.add(pos)

        # Make sure that all sites (except position 1) were tested for all possible derived amino acids.
        pos_to_wt_check = {}
        pos_to_mut_check = defaultdict(set)

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

                aa_combo = (wt_aa, mut_aa)
                
                if pos not in pos_to_wt_check:
                    pos_to_wt_check[pos] = wt_aa
                elif pos_to_wt_check[pos] != wt_aa:
                    sys.exit("Error: Different reference amino acid for same position in " + dms_file + " at position " + str(pos) + " expected " + pos_to_wt_check[pos] + " got " + wt_aa)

                if mut_aa not in pos_to_mut_check[pos]:
                    pos_to_mut_check[pos].add(mut_aa)
                else:
                    sys.exit("Error: Duplicate DMS test for " + dms_file + " at position " + str(pos) + " for amino acid " + mut_aa)

                if pos in buried_sites:
                    buried_aa_combos.append(aa_combo)
                    buried_raw.append(dms_score)
                elif pos in exposed_sites:
                    exposed_aa_combos.append(aa_combo)
                    exposed_raw.append(dms_score)

        if rsa_num_sites != seq_len:
            sys.exit("Error: Number of sites in RSA file does not match DMS file: " + rsa_file + " " + dms_file)

        # Get ranking of all values.
        buried_ranks = list(np.argsort(buried_raw) + 1)
        exposed_ranks = list(np.argsort(exposed_raw) + 1)
        buried_exp = compute_ranks_beta_cdf(buried_ranks)
        exposed_exp = compute_ranks_beta_cdf(exposed_ranks)

        # for i in range(len(buried_exp)):
        #     print(buried_exp[i], 'buried', ','.join(list(buried_aa_combos[i])))
        # for i in range(len(exposed_exp)):
        #     print(exposed_exp[i], 'exposed', ','.join(list(exposed_aa_combos[i])))

        # Check for extreme values
        for i in range(len(buried_exp)):
            if buried_exp[i] < 0 or buried_exp[i] > 10:
                print('wonky: ', dms_id, 'buried_rank', buried_ranks[i], 'out of ', len(buried_ranks), ' buried_exp', buried_exp[i], 'buried_raw', buried_raw[i], 'buried_aa_combo', buried_aa_combos[i], file=sys.stderr)
        working_gene_buried = defaultdict(list)
        for buried_i in range(len(buried_aa_combos)):
            working_gene_buried[buried_aa_combos[buried_i]].append(buried_exp[buried_i])
        for aa_combo in working_gene_buried.keys():
            buried_gene_mean[aa_combo].append(np.mean(working_gene_buried[aa_combo]))
            buried_gene_median[aa_combo].append(np.median(working_gene_buried[aa_combo]))
        
        working_gene_exposed = defaultdict(list)
        for exposed_i in range(len(exposed_aa_combos)):
            working_gene_exposed[exposed_aa_combos[exposed_i]].append(exposed_exp[exposed_i])
        for aa_combo in working_gene_exposed.keys():
            exposed_gene_mean[aa_combo].append(np.mean(working_gene_exposed[aa_combo]))
            exposed_gene_median[aa_combo].append(np.median(working_gene_exposed[aa_combo]))

    header = ['ref_aa',
              'mut_aa',
              'exposed_mean',
              'exposed_median',
              'N_exposed',
              'buried_mean',
              'buried_median',
              'N_buried',
              'weighted_mean',
              'weighted_median']

    print("\t".join(header))

    for aa1 in sorted(list(possible_aa)):
        for aa2 in sorted(list(possible_aa)):
            if aa1 == aa2:
                continue
            aa_combo = (aa1, aa2)

            if aa1 not in possible_aa or aa2 not in possible_aa:
                sys.exit("Error: Unexpected amino acid in DMS file: " + str(aa_combo))
            if aa1 == aa2:
                sys.exit("Error: Unexpected amino acid pair in DMS file: " + str(aa_combo))
            
            if len(exposed_gene_mean[aa_combo]) == 0:
                sys.exit(f"Error: No exposed data for {aa_combo}, cannot compute mean")

            if len(buried_gene_mean[aa_combo]) == 0:
                sys.exit(f"Error: No buried data for {aa_combo}, cannot compute mean")

            exposed_mean = np.mean(exposed_gene_mean[aa_combo])
            buried_mean = np.mean(buried_gene_mean[aa_combo])

            exposed_median = np.median(exposed_gene_mean[aa_combo])
            buried_median = np.median(buried_gene_mean[aa_combo])

            if exposed_mean <= 0 or buried_mean <= 0:
                sys.exit("Error: Unexpected mean value computed for AA combo: " + str(aa_combo) + " exposed mean: " + str(exposed_mean) + " buried mean: " + str(buried_mean))

            background_exposed = mean_background_exposed[aa1]
            background_buried = 1 - background_exposed

            weighted_mean = (background_exposed * exposed_mean + background_buried * buried_mean)
            weighted_median = (background_exposed * exposed_median + background_buried * buried_median)

            outline = [aa1, aa2,
                        f"{exposed_mean:.4f}",
                        f"{exposed_median:.4f}",
                        str(len(exposed_gene_mean[aa_combo])),
                        f"{buried_mean:.4f}",
                        f"{buried_median:.4f}",
                        str(len(buried_gene_mean[aa_combo])),
                        f"{weighted_mean:.4f}",
                        f"{weighted_median:.4f}"]
            print("\t".join(outline))
