#!/usr/bin/env python

import textwrap
import os

seqs = {}

for file in os.listdir("/home6/gmdougla/projects/aa_selection/prerun_struc_prefs/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2"):
    if not file.endswith('.csv'):
        continue
    with open("/home6/gmdougla/projects/aa_selection/prerun_struc_prefs/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2/" + file, 'r') as in_fh:
        headerline = in_fh.readline()
        line = in_fh.readline().strip().split(',')
        seqid = line[0] + ';' + line[1]
        seq = line[3]
        seqs[seqid] = seq

# Write the sequences to a fasta file.
with open("/home6/gmdougla/projects/aa_selection/prerun_struc_prefs/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.fasta", 'w') as out_fh:
    for seqid, seq in seqs.items():
        out_fh.write(">" + seqid + "\n")
        out_fh.write(textwrap.fill(seqs[seqid], width=70) + "\n")
