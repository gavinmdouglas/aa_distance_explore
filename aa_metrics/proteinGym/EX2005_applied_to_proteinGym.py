import os
import sys
from collections import defaultdict
import numpy as np
from scipy.integrate import quad
from scipy import stats
import csv

# EX2005 (Yampolsky and Stoltzflus 2005. Genetics. https://doi.org/10.1534/genetics.104.039107)
# approach for inferring substitution impacts, but applied to the ProteinGym DMS dataset.

# Power law parameters (inferred in the EX2005 paper)
a = 0.278
b = 0.666

def cumul_power(t):
    """Cumulative distribution function"""
    return (t ** b) / (a + t ** b)


def inverse_cumul_power(p):
    """Inverse CDF (quantile function)"""
    return ((p * a) / (1 - p)) ** (1/b)


def pdf_power(t):
    """Probability density function"""
    return (a * b * (t ** (b-1))) / ((a + t ** b) ** 2)


def conditional_mean_in_range(L, U):
    """
    Calculate conditional mean within range [L, U]
    """
    # Numerator: integral of t * f(t)
    def integrand_numerator(t):
        return t * pdf_power(t)
    
    numerator, _ = quad(integrand_numerator, L, U)
    
    # Denominator: integral of f(t)
    denominator, _ = quad(pdf_power, L, U)
    
    return numerator / denominator


def ranks_to_conditional_means(ranks):
    """
    Convert ranks to expected values using conditional means within rank bins

    Expects ranks to be a list of integers representing ranks (1-based, meaning lowest rank is 1).
    """
    ranks = np.array(ranks)
    n = len(ranks)
    
    lower_percentiles = (ranks - 1) / n
    upper_percentiles = ranks / n
    
    # Steps added to make sure integration does not run into numerical issues
    max_percentile = 0.9999
    min_percentile = 0.0001
    lower_percentiles = np.clip(lower_percentiles, min_percentile, max_percentile)
    upper_percentiles = np.clip(upper_percentiles, min_percentile, max_percentile)
    mask = upper_percentiles < lower_percentiles
    upper_percentiles[mask] = lower_percentiles[mask] + 0.009
    
    # Clip again after the adjustment to prevent values > 1.0
    upper_percentiles = np.clip(upper_percentiles, min_percentile, max_percentile)

    # For values at edges, at max or min, manually set defined bounds.
    edge_gap = 0.009

    max_mask = upper_percentiles >= max_percentile
    upper_percentiles[max_mask] = max_percentile
    lower_percentiles[max_mask] = max_percentile - edge_gap

    min_mask = lower_percentiles <= min_percentile
    lower_percentiles[min_mask] = min_percentile
    upper_percentiles[min_mask] = min_percentile + edge_gap

    # Convert percentiles to threshold values
    lower_thresholds = np.array([inverse_cumul_power(p) for p in lower_percentiles])
    upper_thresholds = np.array([inverse_cumul_power(p) for p in upper_percentiles])
    
    # Calculate conditional mean for each rank bin
    conditional_means = np.array([
        conditional_mean_in_range(L, U) 
        for L, U in zip(lower_thresholds, upper_thresholds)
    ])
    
    return conditional_means


# Two functions used for testing.
def find_percentile_threshold(percentile):
    """Find threshold t corresponding to given percentile (0-1)"""
    return inverse_cumul_power(percentile)


def mean_below_percentile(percentile):
    """Calculate mean value from 0 to given percentile"""
    t_threshold = find_percentile_threshold(percentile)
    return conditional_mean_in_range(0, t_threshold)


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
        buried_exp = ranks_to_conditional_means(buried_ranks)
        exposed_exp = ranks_to_conditional_means(exposed_ranks)

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
