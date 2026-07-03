#!/usr/bin/python

# Compute background frequencies for all amino acids for % sites exposed (>= 0.2 RSA).
# Compute separately for humans, E. coli, and Drosophila, and compute mean and median as well.

import os
import sys
from collections import defaultdict

amino_acids = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}

amino_acids_ordered = sorted(list(amino_acids))

def compute_rsa_frequencies(subfolder):
    aa_totals = defaultdict(int)
    aa_exposed = defaultdict(int)
    
    rsa_files = [f for f in os.listdir(subfolder) if f.endswith('.csv')]
    for rsa_file in rsa_files:
        with open(os.path.join(subfolder, rsa_file), 'r') as rsa_fh:
            header = rsa_fh.readline().strip().split(',')
            if header[0] != 'pdb_position' or header[3] != 'pdb_aa' or header[4] != 'rsa':
                sys.exit(f"Error: Unexpected header in RSA file: {rsa_file} {header}")

            for line in rsa_fh:
                line = line.strip().split(',')
                pos = int(line[0])
                aa = line[3]
                rsa = float(line[4])

                # Skip start position and non-standard amino acids.
                if pos == 1 or aa not in amino_acids:
                    continue

                aa_totals[aa] += 1
                if rsa >= 0.2:
                    aa_exposed[aa] += 1

    aa_frequencies = {}
    for aa in amino_acids_ordered:
        aa_frequencies[aa] = aa_exposed[aa] / aa_totals[aa]

    return aa_frequencies


ecoli_percent_exposed = compute_rsa_frequencies('/Users/gavin/Downloads/alphafold_models/UP000000625_83333_ECOLI_v4_dssp_rsa')
human_percent_exposed = compute_rsa_frequencies('/Users/gavin/Downloads/alphafold_models/UP000005640_9606_HUMAN_v4_dssp_rsa')
yeast_percent_exposed = compute_rsa_frequencies('/Users/gavin/Downloads/alphafold_models/UP000002311_559292_YEAST_v4_dssp_rsa')

print("aa\tyeast\tecoli\thuman\tmean\tmedian")
for aa in amino_acids_ordered:
    yeast_freq = yeast_percent_exposed[aa]
    ecoli_freq = ecoli_percent_exposed[aa]
    human_freq = human_percent_exposed[aa]
    mean_freq = (yeast_freq + ecoli_freq + human_freq) / 3
    median_freq = sorted([yeast_freq, ecoli_freq, human_freq])[1]

    print(f"{aa}\t{yeast_freq:.4f}\t{ecoli_freq:.4f}\t{human_freq:.4f}\t{mean_freq:.4f}\t{median_freq:.4f}")
