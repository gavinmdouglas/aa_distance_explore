import os
import sys
import pandas as pd
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from functions import read_fasta, possible_aa

# Suppress future warnings.
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Parse RaSP scores to amino acid preferences for input to phydms.
if not os.path.exists("/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/prefs/rasp/"):
    os.makedirs("/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/prefs/rasp/")

rasp_to_prot = {}
with open('/home6/gmdougla/projects/aa_selection/human_variants/gencode/rasp_gencode_map.tsv', 'r') as map_fh:
    map_fh.readline()
    for line in map_fh:
        line = line.strip().split("\t")
        rasp_to_prot[line[0].split(';')[1]] = line[1].split('|')[0]

prot_seqs = read_fasta('/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/gencode_w_snv.faa')

for rasp_id in rasp_to_prot.keys():
    rasp_predfile = '/home6/gmdougla/projects/aa_selection/prerun_struc_prefs/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2/' + 'rasp_pred_' + rasp_id + 'A.csv'
    
    prot_id = rasp_to_prot[rasp_id]
    if prot_id not in prot_seqs.keys():
        continue

    prot_seq = prot_seqs[prot_id]
    pos_set = list(range(1, len(prot_seq) + 1))
    pos_set = map(str, pos_set)
    max_pos = len(prot_seq) + 1

    pos_to_refaa = dict()
    aa_pref_df = pd.DataFrame(index=pos_set, columns=sorted(list(possible_aa)))

    rasp_predfile = '/home6/gmdougla/projects/aa_selection/prerun_struc_prefs/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2/' + 'rasp_pred_' + rasp_id + 'A.csv'
    with open(rasp_predfile, 'r') as rasp_pred_fh:
        pred_header = rasp_pred_fh.readline().rstrip().split(',')
        col_map = {}
        for i in range(len(pred_header)):
            col_map[pred_header[i]] = i
        for predline in rasp_pred_fh:
            predline = predline.rstrip().split(',')
            variant = predline[col_map['variant']]
            ref_aa = variant[0]
            pos = variant[1:-1]
            alt_aa = variant[-1]
            score = float(predline[col_map['score_ml']])

            if ref_aa not in possible_aa or alt_aa not in possible_aa:
                sys.exit("Error: Unrecognized amino acid in DMS file: " + ref_aa + " or " + alt_aa)

            if pos in pos_to_refaa.keys():
                if ref_aa != pos_to_refaa[pos]:
                        sys.exit("Error: Different reference amino acid for same position in " + rasp_predfile)
                else:
                    pos_to_refaa[pos] = ref_aa
            
            if int(pos) > max_pos:
                sys.exit("Error: Position greater than sequence length in " + rasp_predfile)
        
            aa_pref_df.at[pos, alt_aa] = score

    # Make sure there are non NaN values, otherwise throw error.
    if aa_pref_df.isnull().values.any():
        sys.exit("Error: NaN values in " + rasp_predfile)

    aa_pref_df.index.name = "site"

    aa_pref_df.to_csv('/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/prefs/rasp/' + prot_id + '.csv')
