#!/usr/bin/env python

import sys
import os

# Quick check that the FASTAs split into 32 chunks match exactly the original FASTA.

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from functions import read_fasta

orig = read_fasta('/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/gencode_w_snv.faa')

split_files = os.listdir('/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/gencode_w_snv_split')

parsed_ids = set()
num_confirmed = 0
for split_file in split_files:
    split = read_fasta(os.path.join('/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/gencode_w_snv_split', split_file))
    for id in split:
        if id in parsed_ids:
            sys.exit('Error: ' + id + ' found in multiple split files.')

        if id in orig and split[id] == orig[id]:
            num_confirmed += 1
        else:
            sys.exit('Error: ' + id + ' not found in original or sequences do not match.')
        parsed_ids.add(id)

if num_confirmed == len(orig):
    print('All sequences confirmed.')
else:
    sys.exit('Error: ' + str(len(orig) - num_confirmed) + ' sequences not confirmed.')
