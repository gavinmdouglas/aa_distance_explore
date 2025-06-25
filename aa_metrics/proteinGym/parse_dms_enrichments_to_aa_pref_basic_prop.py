import os
import sys
from collections import defaultdict
import pandas as pd
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import read_fasta, write_fasta, possible_aa

# Suppress future warnings.
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Parse DMS enrichment values to amino acid preferences for input to phydms.
# Read through the DMS CSVs to get the tested positions.
# NOTE: This version of the script output to be a basic proportion.

dms_folder = "/home6/gmdougla/projects/folding_selection/proteinGym_parsing/proteinGym/DMS_ProteinGym_substitutions/"
dms_files = os.listdir(dms_folder)

# First read through EMBER3D test files to get tested positions to parse.
site_files = os.listdir("/home6/gmdougla/projects/folding_selection/proteinGym_parsing/ember3d_prepped_pos_files/")

for site_file in site_files:
    uniprot_id = site_file.replace('.tsv', '')
    if "Giacomelli" in uniprot_id or "Kotler" in uniprot_id:
        continue
    tested_sites = set()

    with open("/home6/gmdougla/projects/folding_selection/proteinGym_parsing/ember3d_prepped_pos_files/" + site_file, "r") as ember3d_fh:
        for site_line in ember3d_fh:
            site_line = site_line.strip().split("\t")
            tested_sites.add(int(site_line[0]))

    if len(tested_sites) == 0:
        continue

    uniprot_clean = uniprot_id.replace('.', '_')
    match_dms_files = []
    for f in dms_files:
        if f.startswith(uniprot_clean):
            match_dms_files.append(f)
    
    if len(match_dms_files) == 0:
        print("Skipping UniProt ID: No DMS files found for " + uniprot_clean, file=sys.stderr)
        continue

    clean_tab = []
    ref_preferred = []
    equal_preferred = []
    parsed_dms_file = []
    tested_seq = []
    tested_pos = []

    full_target_seq = read_fasta("/home6/gmdougla/projects/folding_selection/proteinGym_parsing/focal_seqs/individual_target_regions/" + uniprot_clean + ".faa")
    full_target_seq = list(full_target_seq.values())[0]

    pos_to_refaa = dict()

    for matched_dmsfile in match_dms_files:

        aa_pref_df = pd.DataFrame(index=sorted(list(tested_sites)), columns=sorted(list(possible_aa)))

        with open(dms_folder + matched_dmsfile, "r") as dms_fh:
            dms_fh.readline()
            for line in dms_fh:
                line = line.strip().split(",")
                info = line[0]

                # Skip multi-substitution tests.
                if ":" in info:
                    continue

                ref_aa = info[0]
                pos = int(info[1:-1])
                alt_aa = info[-1]

                if pos in pos_to_refaa.keys():
                    if ref_aa != pos_to_refaa[pos]:
                        sys.exit("Error: Different reference amino acid for same position in " + matched_dmsfile)
                else:
                    pos_to_refaa[pos] = ref_aa

                if pos not in tested_sites:
                    continue

                if ref_aa not in possible_aa or alt_aa not in possible_aa:
                    sys.exit("Error: Unrecognized amino acid in DMS file: " + ref_aa + " or " + alt_aa)

                # Check if cell value is non-NaN:
                if not pd.isna(aa_pref_df.loc[pos, ref_aa]):
                    print(dms_folder + matched_dmsfile, file=sys.stderr)
                    sys.exit("Error: Duplicate DMS test for " + uniprot_clean + " at position " + str(pos) + " for amino acid " + ref_aa)
                else:
                    aa_pref_df.loc[pos, alt_aa] = float(line[2])

            # Check that there is exactly one NaN value per row.
            if not aa_pref_df.isnull().sum(axis=1).eq(1).all():
                
                # Figure out the site indices with exactly one NaN value:
                passing_pos = aa_pref_df.index[aa_pref_df.isnull().sum(axis=1).eq(1)]

                # Skip if no rows remain.
                if len(passing_pos) == 0:
                    print("Skipping file: All tested sites do not have full test coverage for " + matched_dmsfile, file=sys.stderr)
                    continue

                aa_pref_df = aa_pref_df.loc[passing_pos]

            else:
                passing_pos = aa_pref_df.index

            # Get control table where the reference amino acid is greatly preferred compared to the others.
            ref_pref_tab = aa_pref_df.copy()
            ref_pref_tab = ref_pref_tab.mask(ref_pref_tab.notna(), 0.01)
            ref_pref_tab = ref_pref_tab.mask(ref_pref_tab.isna(), 0.81)
            ref_pref_tab.index = range(1, ref_pref_tab.shape[0] + 1)
            ref_pref_tab.index.name = "site"
            ref_preferred.append(ref_pref_tab)

            # Also, compare to a scenario where all amino acids are equally preferred.
            equal_pref_tab = ref_pref_tab.copy()
            equal_pref_tab = equal_pref_tab.mask(equal_pref_tab.notna(), 0.05)
            equal_pref_tab.index = range(1, equal_pref_tab.shape[0] + 1)
            equal_pref_tab.index.name = "site"
            equal_preferred.append(equal_pref_tab)

            # Replace NA value per row with the median value per row.
            aa_pref_df = aa_pref_df.apply(lambda x: x.fillna(x.median()), axis=1)

            # Then get sd per row, and add it plus the absolute min value to each row.
            row_sds = aa_pref_df.std(axis=1)
            row_mins = np.abs(aa_pref_df.min(axis=1))
            aa_pref_df = aa_pref_df.add(row_mins + row_sds, axis=0)

            # Then convert to proportions per row.
            aa_pref_df = aa_pref_df.div(aa_pref_df.sum(axis=1), axis=0)

            # Change index to go from 1 to number of rows, and write out.
            aa_pref_df.index = range(1, aa_pref_df.shape[0] + 1)
            aa_pref_df.index.name = "site"
            clean_tab.append(aa_pref_df)
            parsed_dms_file.append(matched_dmsfile)

            # Figure out tested residue positions and concatenate residues of tested sequence.
            if aa_pref_df.shape[0] > len(tested_sites):
                sys.exit("Error: More tested sites than expected for " + uniprot_clean)
            else:
                tested_seq_working = ''
                tested_pos_working = []
                for pos in passing_pos:
                    tested_pos_working.append(str(pos))
                    tested_seq_working += full_target_seq[pos - 1]
                    if full_target_seq[pos - 1] != pos_to_refaa[pos]:
                        sys.exit("Error: Reference amino acid does not match target sequence for " + uniprot_clean)

                tested_seq.append(tested_seq_working)
                tested_pos.append(','.join(tested_pos_working))

    if len(clean_tab) == 0:
        sys.exit("Error: No DMS could be successfully parsed for " + uniprot_clean)
    elif len(clean_tab) == 1:
        clean_tab[0].to_csv("/home6/gmdougla/projects/folding_selection/phydms_proteinGym_predict/dms_prefs_basic_prop/" + uniprot_clean + ".csv", sep=",")
        ref_preferred[0].to_csv("//home6/gmdougla/projects/folding_selection/phydms_proteinGym_predict/dms_prefs_controls/ref_preferred/" + uniprot_clean + ".csv", sep=",")
        equal_preferred[0].to_csv("/home6/gmdougla/projects/folding_selection/phydms_proteinGym_predict/dms_prefs_controls/equal_preferred/" + uniprot_clean + ".csv", sep=",")

    else:
        for i in range(len(clean_tab)):
            seq_id = uniprot_clean + "." + str(i + 1)
            clean_tab[i].to_csv("/home6/gmdougla/projects/folding_selection/phydms_proteinGym_predict/dms_prefs_basic_prop/" + seq_id + ".csv", sep=",")
            ref_preferred[i].to_csv("/home6/gmdougla/projects/folding_selection/phydms_proteinGym_predict/dms_prefs_controls/ref_preferred/" + seq_id + ".csv", sep=",")
            equal_preferred[i].to_csv("/home6/gmdougla/projects/folding_selection/phydms_proteinGym_predict/dms_prefs_controls/equal_preferred/" + seq_id + ".csv", sep=",")

