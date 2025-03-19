#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from functions import read_fasta, write_fasta

# There are X residues in the sequences, which cause issues with VespaG.
# Replace all X residues with 'M' when they are at the start of the sequence (exclusively for VespaG!).
# These sites should not be actively analyzed though, and should be skipped.
# Throw an error if any X residues are found anywhere but the start of a sequence.

infasta = read_fasta('/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/gencode_w_snv.faa')

final_seqs = {}

for seqid, seq in infasta.items():
    if seq[0] == 'X':
        final_seqs[seqid] = 'M' + seq[1:]
    else:
        final_seqs[seqid] = seq
    if 'X' in final_seqs[seqid]:
        raise ValueError('X residues found in sequence: %s' % final_seqs[seqid])

write_fasta(final_seqs, '/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/gencode_w_snv_noX_temp.faa')
