#!/usr/bin/env python

import gzip
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from functions import read_fasta, write_fasta

# Identify protein sequences with RaSP predictions that are exact matches to proteins in Gencode.
rasp = read_fasta('/home6/gmdougla/projects/aa_selection/prerun_struc_prefs/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.fasta.gz')
gencode = read_fasta('/home6/gmdougla/projects/aa_selection/human_variants/gencode/gencode.v45.pc_translations.fa.gz')

rasp_seqs = set(list(rasp.values()))

rasp_seq_to_id = {}
for rasp_id, rasp_seq in rasp.items():
    # Note that there are some RaSP sequences that are identical to each other,
    # so just kept the last one parsed.
    #if rasp_seq in rasp_seq_to_id.keys():
        #print(rasp_id + ' is identical to another RaSP sequence.', file=sys.stderr)
    rasp_seq_to_id[rasp_seq] = rasp_id

intersecting = {}
id_map = dict()
for gencode_id, gencode_seq in gencode.items():
    if gencode_seq in rasp_seqs:
        rasp_id = rasp_seq_to_id[gencode_seq]
        if rasp_id in id_map.keys():
            # More than one Gencode sequence is an exact match for this RaSP ID.
            # So just keep one and skip others.
            continue
        else:
            id_map[rasp_id] = gencode_id

        intersecting[gencode_id] = gencode_seq

parsed_genes = set()
prot_to_exclude = set()
with open('/home6/gmdougla/projects/aa_selection/human_variants/gencode/rasp_gencode_map.tsv', 'w') as out_fh:
    out_fh.write('RaSP\tGencode\n')
    for rasp_id, gencode_id in id_map.items():
        gene_id = gencode_id.split('|')[2]
        if gene_id in parsed_genes:
            prot_to_exclude.add(gencode_id)
            continue
        else:
            parsed_genes.add(gene_id)
            out_fh.write(rasp_id + '\t' + gencode_id + '\n')

print('Excluding ' + str(len(prot_to_exclude)) + ' proteins from Gencode as they correspond to genes already represented.', file=sys.stderr)

final_out = {}
for seqid, seq in intersecting.items():
    if seqid in prot_to_exclude:
        continue
    final_out[seqid.split('|')[0]] = seq

write_fasta(final_out, '/home6/gmdougla/projects/aa_selection/human_variants/gencode/gencode.v45.pc_translations_rasp.fa')
