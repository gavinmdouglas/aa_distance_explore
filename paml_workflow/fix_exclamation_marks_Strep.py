#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import read_fasta, write_fasta

seq19 = read_fasta("/home6/gmdougla/projects/aa_distance/paml_workflow/concat_seqs/Strep/subsample19.fna")

for seq_id in list(seq19.keys()):
    for i in range(len(seq19[seq_id])):
        seq_i = seq19[seq_id][i]
        if seq_i not in ['A', 'C', 'G', 'T', '-']:
            print("Weird character: " + seq_i + " in sequence " + seq_id + ", will be replaced with '-'.")
            seq19[seq_id] = seq19[seq_id][:i] + '-' + seq19[seq_id][i + 1:]
write_fasta(seq19, "/home6/gmdougla/projects/aa_distance/paml_workflow/concat_seqs/Strep/subsample19.fna")

seq6 = read_fasta("/home6/gmdougla/projects/aa_distance/paml_workflow/concat_seqs/Strep/subsample6.fna")

for seq_id in list(seq6.keys()):
    for i in range(len(seq6[seq_id])):
        seq_i = seq6[seq_id][i]
        if seq_i not in ['A', 'C', 'G', 'T', '-']:
            print("Weird character: " + seq_i + " in sequence " + seq_id + ", will be replaced with '-'.")
            seq6[seq_id] = seq6[seq_id][:i] + '-' + seq6[seq_id][i + 1:]
write_fasta(seq6, "/home6/gmdougla/projects/aa_distance/paml_workflow/concat_seqs/Strep/subsample6.fna")
