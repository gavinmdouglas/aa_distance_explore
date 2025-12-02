#!/usr/bin/python3

import gzip
import sys
import textwrap
import pandas as pd
import numpy as np
import os
from collections import defaultdict
import argparse

bases = set(['A', 'C', 'G', 'T'])

code11_codon_to_aa = {'TTT':'F', 'TTC':'F',
                       'TTA':'L', 'TTG':'L',
                       'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                       'TAT':'Y', 'TAC':'Y',
                       'TAA':'*', 'TAG':'*',
                       'TGT':'C', 'TGC':'C',
                       'TGA':'*',
                       'TGG':'W',
                       'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                       'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                       'CAT':'H', 'CAC':'H',
                       'CAA':'Q', 'CAG':'Q',
                       'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                       'ATT':'I', 'ATC':'I', 'ATA':'I',
                       'ATG':'M',
                       'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                       'AAT':'N', 'AAC':'N',
                       'AAA':'K', 'AAG':'K',
                       'AGT':'S', 'AGC':'S',
                       'AGA':'R', 'AGG':'R',
                       'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                       'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                       'GAT':'D', 'GAC':'D',
                       'GAA':'E', 'GAG':'E',
                       'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

# As listed here: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11
code11_possible_start_codons = set(['ATG', 'ATT', 'ATC', 'ATA', 'CTG', 'GTG', 'TTG'])

possible_aa = set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
sorted_aa = sorted(list(possible_aa))

ambig_aa_map = {"B": set(["D", "N"]),
                "Z": set(["E", "Q"]),
                "J": set(["I", "L"])}

def parse_str_arg_to_floats(value):
    try:
        floats = [float(x) for x in value.split(',')]
        return floats
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid float value: '{value}'")


def vespag_dir_to_prefs(in_dir, out_dir, all_sites_present=True, ref_aa_file=None, ref_fill="max", sum_scale=True):
    '''
    Given a VESPA-G output directory, output all amino acid preference tables.
    all_sites_present specifies whether all sites from min to max pos need to be present.
    ref_aa_file is an optional file that can be created with mappings of reference AA at each site. Can be useful for troubleshooting.
    ref_fill specifies the value to use for the reference AA. Must be one of "one", "max", "median", or "NA".
    '''

    if not os.path.exists(in_dir):
        sys.exit('Input directory does not exist.')

    # Read in all files in specified directory.
    files = os.listdir(in_dir)
    if len(files) == 0:
        sys.exit('No files found in input directory.')

    # Create output directory, if it does not exist.
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Optionally output a table of the reference AA at each site, which can be useful for sanity checks.
    ref_aa_tracker = []
    parsed_genes = set()
    for filename in files:
        if not filename.endswith('.csv.gz'):
            continue
        gene = filename.replace('.csv.gz', '')

        if gene in parsed_genes:
            sys.exit("Error: VESPA-G outfile for gene " + gene + " already parsed.")
        else:
            parsed_genes.add(gene)

        filename_path = in_dir + '/' + filename
        if not os.path.exists(filename_path):
            sys.exit('File ' + filename_path + ' does not exist.')

        # Read through file once, just to get the number of positions.
        unique_pos = set()
        with gzip.open(filename_path, 'rt') as vespag_in:
            vespag_in.readline()
            for line in vespag_in:
                line = line.rstrip().split(',')
                pos_info = line[0]
                pos = int(pos_info[1:-1])
                unique_pos.add(pos)

        sorted_unique_pos = sorted(list(unique_pos))

        if all_sites_present:
            max_pos = max(unique_pos)
            if len(unique_pos) != max_pos:
                sys.exit('Max pos of ' + str(max_pos) + ' does not match number of unique positions.')

        # Initialize preference table.
        prefs = pd.DataFrame(np.nan,
                             index=sorted_unique_pos,
                             columns=sorted_aa)

        ref_AAs = set()
        obs_pos = set()
        derived_AAs = []
        inferred_val = []

        with gzip.open(filename_path, 'rt') as vespag_in:
            vespag_in.readline()

            # Read in 19 lines at a time.
            for line in vespag_in:
                linesplit = line.rstrip().split(',')
                pos_info = linesplit[0]
                obs_pos.add(int(pos_info[1:-1]))
                ref_AAs.add(pos_info[0])
                derived_AAs += [pos_info[-1]]

                # Make sure value is between 0 and 1.
                if float(linesplit[1]) <= 0 or float(linesplit[1]) >= 1:
                    sys.exit('Preference value not between 0 and 1: ' + linesplit[1])
                inferred_val.append(float(linesplit[1]))

                if len(derived_AAs) == 19 and len(inferred_val) == 19:

                    # Sanity check that only one ref AA and one pos parsed.
                    if len(ref_AAs) != 1 or len(obs_pos) != 1:
                        sys.exit('Multiple ref AAs or positions parsed in one block.')

                    # Fill in preference table.
                    pos = list(obs_pos)[0]
                    ref_AA = list(ref_AAs)[0]
                    for i in range(19):
                        prefs.loc[pos, derived_AAs[i]] = inferred_val[i]

                    # Set ref. AA to median or max value of all derived AAs.
                    if ref_fill == "median":
                        prefs.loc[pos, ref_AA] = np.median(inferred_val)
                    elif ref_fill == "max":
                        prefs.loc[pos, ref_AA] = np.max(inferred_val)
                    elif ref_fill == "one":
                        prefs.loc[pos, ref_AA] = 1.0
                    elif ref_fill == "NA":
                        prefs.loc[pos, ref_AA] = np.nan
                    else:
                        sys.exit("Error: ref_fill option must be one, NA, median, or max.")

                    if ref_aa_file is not None:
                        ref_aa_tracker += [gene + "\t" + str(pos) + "\t" + ref_AA + "\n"]

                    ref_AAs = set()
                    obs_pos = set()
                    derived_AAs = []
                    inferred_val = []

        if sum_scale:
            # Total sum scale each row.
            prefs = prefs.div(prefs.sum(axis=1), axis=0)

        if ref_fill != "NA":
            # Throw error if any NaNs present.
            if prefs.isnull().values.any():
                print(prefs)
                sys.exit('NaNs present in preference table.')

        # Output preference table to file.
        prefs.index.name = 'site'
        prefs.to_csv(out_dir + '/' + gene + '.csv', sep=',', float_format='%.4f')

    if ref_aa_file is not None:
        with open(ref_aa_file, 'w') as ref_aa_out:
            ref_aa_out.write("gene\tsite\tref_AA\n")
            for line in ref_aa_tracker:
                ref_aa_out.write(line)

    return None


def seq_sub_prob_prefs(seq, effect_tab, outfile=None):
    '''
    Given a sequence and a dictionary of the prob. of each AA sub. type,
    return an amino acid preference table.
    '''
    prefs = pd.DataFrame(np.nan,
                        index=range(1, 1, len(seq) + 1),
                        columns=sorted_aa)

    for i in range(len(seq)):
        seq_aa = seq[i]
        if seq_aa not in possible_aa:
            continue
        for sub_aa in sorted_aa:
            prefs.loc[i + 1, sub_aa] = effect_tab[(seq_aa, sub_aa)]

    # Total sum scale each row.
    prefs = prefs.div(prefs.sum(axis=1), axis=0)

    prefs.index.name = 'site'

    if outfile is not None:
        prefs.to_csv(outfile, sep=',', float_format='%.4f', na_rep='')

    return prefs


def read_rad_vs_cons_effects(sub_file,
                             swap_rad_cons=False,
                             rad_pref=None,
                             set_identity_to_cons=True):
    '''
    Read through a sub. effect file and return a dictionary
    of radical vs conservative effects (where each sub pair is a tuple).
    Optionally swap in weights for radical vs. conservative values.
    '''

    if swap_rad_cons:
        if rad_pref is None:
            sys.exit('Error: rad not provided, but needed when swap_rad_cons=True.')
        else:
            cons_pref = 1 - rad_pref
    elif rad_pref is not None:
        sys.exit('Error: rad_pref provided, but not needed when swap_rad_cons=False.')

    sub_effects = {}
    with gzip.open(sub_file, 'rt') as sub_in:
        headerline = sub_in.readline()
        if headerline != "aa1\taa2\tsub_type\n":
            print(headerline, file=sys.stderr)
            sys.exit('Error: Sub. effect file does not have expected header.')

        for line in sub_in:
            aa1, aa2, sub_type = line.strip().split()
            if swap_rad_cons:
                if sub_type == 'Radical':
                    sub_effects[(aa1, aa2)] = rad_pref
                elif sub_type == 'Conservative':
                    sub_effects[(aa1, aa2)] = cons_pref
            else:
                sub_effects[(aa1, aa2)] = sub_type
    
    if set_identity_to_cons:
        for aa in sorted_aa:
            sub_effects[(aa, aa)] = cons_pref

    return sub_effects


def read_sim_matrix(sub_sim_file, identity_set='one'):
    '''
    Read through a substitution similarity file and return a dictionary
    of substitution matrix values (where each AA sub pair is a tuple).
    Can replace the self substitutions with 1.0 ("one"), or the "max" or the "median" of all other values (or to be "NA").
    '''

    if identity_set != 'max' and identity_set != 'median' and identity_set != 'NA' and identity_set != 'one':
        sys.exit('Error: "identity_set" must be "one", "max", "median", or "NA".')

    per_aa_sub_effects = defaultdict(list)
    sub_effects = {}
    with gzip.open(sub_sim_file, 'rt') as sub_in:
        headerline = sub_in.readline()
        if headerline[:8] != "aa1\taa2\t":
            print(headerline, file=sys.stderr)
            sys.exit('Error: Sub. effect file does not have expected header.')

        for line in sub_in:
            aa1, aa2, sim = line.strip().split()
            sim = float(sim)
            sub_effects[(aa1, aa2)] = sim
            per_aa_sub_effects[aa1].append(sim)

        if identity_set == 'max':
            for aa in possible_aa:
                sub_effects[(aa, aa)] = np.max(per_aa_sub_effects[aa])
        elif identity_set == 'median':
            for aa in possible_aa:
                sub_effects[(aa, aa)] = np.median(per_aa_sub_effects[aa])
        elif identity_set == 'NA':
            for aa in possible_aa:
                sub_effects[(aa, aa)] = np.nan
        elif identity_set == 'one':
            for aa in possible_aa:
                sub_effects[(aa, aa)] = 1.0

    return sub_effects


def read_genetic_code(code_file):
    '''
    Read in genetic code file in format on this repo: https://github.com/gavinmdouglas/codon_tables/.
    E.g., see https://github.com/gavinmdouglas/codon_tables/blob/main/genetic_codes/code1.tsv.
    
    Return tuple with:
        1. Dictionary with codon as key and amino acid as value.
        2. Set of possible start codons.
    '''
    with open(code_file, 'r') as code_in:
        codon_to_aa = {}
        possible_start_codons = set()
    
        for line in code_in:
            line = line.rstrip()
            linesplit = line.split('\t')

            if linesplit[0] == 'codon':
                continue

            codon_to_aa[linesplit[0]] = linesplit[1]

            if linesplit[2] == 'yes':
                possible_start_codons.add(linesplit[0])

        return (codon_to_aa, possible_start_codons)


def translate_DNA(in_dna: str, codon_to_aa: dict, possible_start_codons=None, check_start=True, silent=False, ambig_char='X'):
    '''
    Translate DNA sequence to amino acid sequence. Require genetic code info.
    '''

    # Check that input sequence is divisible by three.
    if not silent and len(in_dna) % 3 != 0:
        print('Warning: gene sequence is not perfectly divisible by three.',
              file = sys.stderr)
    
    # Initialize empty string to hold translated sequence.
    aa_seq = ''

    initial_i = 0

    if check_start:
        initial_i = 3
        if possible_start_codons is None:
            sys.exit('Error: possible_start_codons not provided.')

        # Do start codon separately.
        start_codon = in_dna[:3]
        
        if start_codon in codon_to_aa.keys():
            if start_codon in possible_start_codons:
                aa_seq += 'M'
            else:
                sys.exit('Error - first codon is not a possible start codon: ' + start_codon + 
                        ' in sequence: ' + in_dna)
        else:
            if not silent:
                print('Skipping first codon as it does not match any expected triplet (replacing it with ambig character, or gap if it was "---"): ' + start_codon, file=sys.stderr)
            if start_codon == '---':
                aa_seq += '-'
            else:
                codon_aa = ambig_char
            

    # Iterate over input DNA sequence, translating each codon to amino acid.
    for i in range(initial_i, len(in_dna), 3):
        codon = in_dna[i:i + 3]

        if codon in codon_to_aa.keys():
            codon_aa = codon_to_aa[codon]
        else:
            if not silent:
                print("Skipping this codon (and filling it in with an ambig character, or gap if it was '---') as it does not match any expected triplet: " + codon, file=sys.stderr)
            if codon == '---':
                codon_aa = '-'
            else:
                codon_aa = ambig_char
        aa_seq += codon_aa

    return aa_seq


def read_fasta(filename, cut_header=False, convert_upper=True):
    '''Read in FASTA file (gzipped or not) and return dictionary with each
    independent sequence id as a key and the corresponding sequence string as
    the value.'''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    def parse_fasta_lines(file_handle):
        for line in file_handle:

            line = line.rstrip()

            if len(line) == 0:
                continue

            if line[0] == ">":

                if cut_header:
                    name = line.split()[0][1:]
                else:
                    name = line[1:]

                name = name.rstrip("\r\n")

                # Make sure that sequence ID is not already in dictionary.
                if name in seq:
                    sys.stderr("Stopping due to duplicated ID in file: " + name)

                # Intitialize empty sequence with this ID.
                seq[name] = ""

            else:
                line = line.rstrip("\r\n")
                if convert_upper:
                    line = line.upper()
                seq[name] += line

    # Read in FASTA line-by-line.
    if filename[-3:] == ".gz":
        with gzip.open(filename, "rt") as fasta_in:
            parse_fasta_lines(fasta_in)
    else:
        with open(filename, "r") as fasta_in:
            parse_fasta_lines(fasta_in)

    return seq


def write_fasta(seq, outfile, id_order=None):
    out_fasta = open(outfile, "w")

    if id_order is None:
         # If not specified, thenn make IDs alphabetically ordered, so output is deterministic.
        id_order = sorted(seq.keys())
        
    for s in id_order:
        out_fasta.write(">" + s + "\n")
        out_fasta.write(textwrap.fill(seq[s], width=70) + "\n")

    out_fasta.close()


def write_fasta_single_seq(sequence, name, outfile):
    out_fasta = open(outfile, 'w')
    out_fasta.write('>' + name + '\n')
    out_fasta.write(textwrap.fill(sequence, width=70) + '\n')
    out_fasta.close()


def identify_potential_subs(in_dna, in_protein):
    '''
    Identify all possible non-synonymous substitutions that could occur at each codon in a protein-coding DNA sequence.
    Will skip the start codon by default. Will compare each translated codon to input protein sequence to address codons matching stop codons.
    '''

    # Replace all ambiguous IUPAC characters with N.
    in_dna = in_dna.replace('R', 'N').replace('Y', 'N').replace('S', 'N').replace('W', 'N').replace('K', 'N').replace('M', 'N').replace('B', 'N').replace('D', 'N').replace('H', 'N').replace('V', 'N')
 
    # Remove any gap characters.
    in_dna = in_dna.replace('-', '').replace('.', '')
    in_protein = in_protein.replace('-', '').replace('.', '')

    potential_aa_subs = list()

    seq_aa_num = 1

    for i in range(0, len(in_dna), 3):

        codon = in_dna[i:i + 3]

        if codon not in codon_to_aa.keys():
            seq_aa_num += 1
            print('Skipping codon as it does not match any expected triplet: ' + codon, file=sys.stderr)
            continue

        obs_aa = codon_to_aa[codon]
        exp_aa = in_protein[seq_aa_num - 1]

        if obs_aa != exp_aa:
            if seq_aa_num == 1 and codon in possible_start_codons and exp_aa == 'M':
                obs_aa = 'M'
            elif codon == 'TGA' and exp_aa == 'U':
                obs_aa = 'U'
         ##   elif codon == 'AAA' and (exp_aa == 'K' or exp_aa == 'R'):
          #      obs_aa = exp_aa
            else:
                print('Error: observed codon does not match expected amino acid: ' + codon + ' ' + obs_aa + ' ' + exp_aa, file=sys.stderr)
                print('Position:', seq_aa_num, file=sys.stderr)
                print('DNA:', in_dna, file=sys.stderr)
                print('Protein:', in_protein, file=sys.stderr)
                sys.exit(1)

        # Keep track of all possible substitutions.
        # Do so by looping over every possible nucleotide change.
        for j in range(3):
            for nt in ['A', 'C', 'G', 'T']:

                if codon[j] == nt:
                    continue

                codon_potential = codon[:j] + nt + codon[j + 1:]

                if codon_potential not in codon_to_aa.keys():
                    sys.exit('Why is potential codon not in codon_to_aa.keys()?')

                aa_potential = codon_to_aa[codon_potential]

                # Skip synonymous mutations and stop codons.
                if aa_potential == obs_aa or aa_potential == '*':
                    continue

                potential_aa_subs += [[str(seq_aa_num), codon, codon_potential, obs_aa, aa_potential]]

        seq_aa_num += 1

    return potential_aa_subs


def codon_ambig_check(codon):
    '''
    Quick check to see why codon did not match to an AA.
    Key interest is to know whether to intecrement the number of AAs or not in the AA alignment (which will be done if it's not all gaps).
    However, if there is a mix of gaps and nucleotides, then 'Other' will be returned, which warrants further investigation.
    '''

    if codon == '---':
        return 'all_gap'
    elif '_' not in codon and 'N' in codon:
        return 'no_gap_N'
    else:
        return 'Other'


def continuous_prefs_dict_sanity(prefs_dict):
    '''
    Make sure that the preferences dictionary includes all AA pairs,
    that the values are all floats, and that they are between 0 and 1.
    '''
    for aa1 in possible_aa:
        for aa2 in possible_aa:
            if (aa1, aa2) not in prefs_dict.keys():
                return False
            if not isinstance(prefs_dict[(aa1, aa2)], float):
                return False
            if not (0 <= prefs_dict[(aa1, aa2)] <= 1 or np.isnan(prefs_dict[(aa1, aa2)])):
                return False

    return True
