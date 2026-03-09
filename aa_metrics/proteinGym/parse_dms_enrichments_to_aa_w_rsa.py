import os
import sys
from collections import defaultdict
import pandas as pd
import numpy as np
import csv
from scipy import stats

possible_aa = set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
unique_aa = sorted(list(possible_aa))

focal_ids = set()
with open('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_filt_dms_ids.txt', 'r') as f:
    for line in f:
        focal_ids.add(line.strip())

dms_to_transform = dict()
with open('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_score_type_unambig.csv', 'r') as f:
    csv_reader = csv.DictReader(f)
    for line in csv_reader:
        dms_to_transform[line['DMS_id']] = line['dms_score_transform']

dms_ids = sorted(list(dms_to_transform.keys()))

mean_per_gene = defaultdict(list)

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
buried_gene_mean_standard_score = defaultdict(list)
exposed_gene_mean_standard_score = defaultdict(list)
buried_gene_median_robust_score = defaultdict(list)
exposed_gene_median_robust_score = defaultdict(list)

for dms_id in dms_ids:
    if dms_id not in focal_ids:
        continue
    site_wt = dict()
    site_to_effect = dict()

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

            if dms_to_transform[dms_id] == "ln":
                dms_score = dms_score
            elif dms_to_transform[dms_id] == "log2":
                dms_score = np.log(2 ** dms_score)
            elif dms_to_transform[dms_id] == "raw ratio":
                if dms_score < 0.01:
                    dms_score = 0.01
                dms_score = np.log(dms_score)

            # Check if nan, and print line if so.
            if np.isnan(dms_score):
                sys.exit("Error: DMS score is NaN in DMS file: " + dms_file + " for line " + str(line))

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

    if rsa_num_sites != seq_len:
        sys.exit("Error: Number of sites in RSA file does not match DMS file: " + rsa_file + " " + dms_file)

    buried_vals = []
    exposed_vals = []
    for site in site_to_effect.keys():
        site_scores = []
        for aa2 in site_to_effect[site].keys():
            site_scores.append(site_to_effect[site][aa2])

        if len(site_scores) < 15:
            continue
        elif site in buried_sites:
            buried_vals = buried_vals + site_scores
        elif site in exposed_sites:
            exposed_vals = exposed_vals + site_scores
        else:
            sys.exit("Error: Site " + str(site) + " not found in either buried or exposed sites for DMS file: " + dms_file)

    if len(buried_vals) >= 10:
        buried_mean_score = np.mean(buried_vals)
        buried_sd_score = np.std(buried_vals)
        buried_median_score = np.median(buried_vals)
        buried_mad_score = stats.median_abs_deviation(buried_vals)
        buried_skip = False
        if buried_mad_score == 0:
            print("Warning: Buried median absolute deviation is zero for DMS file: " + dms_file + ". Will skip this gene.", file=sys.stderr)
            buried_skip = True

    else:
        buried_skip = True
        print("Warning: Not enough buried sites for DMS file: " + dms_file + ". Will skip this gene.", file=sys.stderr)
    if len(exposed_vals) >= 10:
        exposed_mean_score = np.mean(exposed_vals)
        exposed_sd_score = np.std(exposed_vals)
        exposed_median_score = np.median(exposed_vals)
        exposed_mad_score = stats.median_abs_deviation(exposed_vals)
        exposed_skip = False
        if exposed_mad_score == 0:
            print("Warning: Exposed median absolute deviation is zero for DMS file: " + dms_file + ". Will skip this gene.", file=sys.stderr)
            exposed_skip = True
    else:
        exposed_skip = True
        print("Warning: Not enough exposed sites for DMS file: " + dms_file + ". Will skip this gene.", file=sys.stderr)

    exposed_gene_standard_score = defaultdict(list)
    exposed_gene_robust_score = defaultdict(list)
    buried_gene_standard_score = defaultdict(list)
    buried_gene_robust_score = defaultdict(list)

    for site in site_to_effect.keys():
        ref_aa = site_wt[site]

        if site in exposed_sites:
            if exposed_skip:
                continue
            site_mean_score = exposed_mean_score
            site_sd_score = exposed_sd_score
            site_median_score = exposed_median_score
            site_mad_score = exposed_mad_score
            gene_standard_score = exposed_gene_standard_score
            gene_robust_score = exposed_gene_robust_score
        elif site in buried_sites:
            if buried_skip:
                continue
            site_mean_score = buried_mean_score
            site_sd_score = buried_sd_score
            site_median_score = buried_median_score
            site_mad_score = buried_mad_score
            gene_standard_score = buried_gene_standard_score
            gene_robust_score = buried_gene_robust_score

        for aa2 in site_to_effect[site].keys():
            if aa2 == ref_aa:
                sys.exit("Error: Unexpected reference amino acid in DMS file: " + dms_file + " for line " + str(site) + " " + ref_aa)
            else:
                aa_combo = (ref_aa, aa2)
                standard_score = (site_to_effect[site][aa2] - site_mean_score) / site_sd_score
                gene_standard_score[aa_combo].append(standard_score)

                if site_mad_score > 0:
                    robust_score = (site_to_effect[site][aa2] - site_median_score) / (1.4826 * site_mad_score)
                    gene_robust_score[aa_combo].append(robust_score)

    for aa_combo in exposed_gene_standard_score.keys():
        if len(exposed_gene_standard_score[aa_combo]) > 0:
            exposed_gene_mean_standard_score[aa_combo].append(np.mean(exposed_gene_standard_score[aa_combo]))

        if len(exposed_gene_robust_score[aa_combo]) > 0:
            exposed_gene_median_robust_score[aa_combo].append(np.median(exposed_gene_robust_score[aa_combo]))

    for aa_combo in buried_gene_standard_score.keys():
        if len(buried_gene_standard_score[aa_combo]) > 0:
            buried_gene_mean_standard_score[aa_combo].append(np.mean(buried_gene_standard_score[aa_combo]))
        if len(buried_gene_robust_score[aa_combo]) > 0:
            buried_gene_median_robust_score[aa_combo].append(np.median(buried_gene_robust_score[aa_combo]))

header = ['ref_aa', 'mut_aa',
          'exposed_mean_standard_score', 'N_exposed_mean_standard_score',
          'exposed_median_robust_score', 'N_exposed_median_robust_score',
          'buried_mean_standard_score', 'N_buried_mean_standard_score',
          'buried_median_robust_score', 'N_buried_median_robust_score',
          'weighted_mean_standard_score', 'weighted_median_robust_score',
           'weighted_median_standard_score']

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
        
        exposed_mean_score = np.mean(exposed_gene_mean_standard_score[aa_combo])
        exposed_median_robust_score = np.median(exposed_gene_median_robust_score[aa_combo])
        
        buried_mean_score = np.mean(buried_gene_mean_standard_score[aa_combo])
        buried_median_robust_score = np.median(buried_gene_median_robust_score[aa_combo])

        background_exposed = mean_background_exposed[aa1]
        background_buried = 1 - background_exposed

        weighted_mean_standard_score = (background_exposed * exposed_mean_score + background_buried * buried_mean_score)
        weighted_median_robust_score = (background_exposed * exposed_median_robust_score + background_buried * buried_median_robust_score)

        exposed_median_score = np.median(exposed_gene_mean_standard_score[aa_combo])
        buried_median_score = np.median(buried_gene_mean_standard_score[aa_combo])

        weighted_median_standard_score = (background_exposed * exposed_median_score + background_buried * buried_median_score)

        outline = [aa1, aa2,
                f"{exposed_mean_score:.4f}", str(len(exposed_gene_mean_standard_score[aa_combo])),
                    f"{exposed_median_robust_score:.4f}", str(len(exposed_gene_median_robust_score[aa_combo])),
                    f"{buried_mean_score:.4f}", str(len(buried_gene_mean_standard_score[aa_combo])),
                    f"{buried_median_robust_score:.4f}", str(len(buried_gene_median_robust_score[aa_combo])),
                    f"{weighted_mean_standard_score:.4f}", f"{weighted_median_robust_score:.4f}",
                    f"{weighted_median_standard_score:.4f}"]
        print("\t".join(outline))
