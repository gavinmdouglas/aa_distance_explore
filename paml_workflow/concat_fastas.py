#!/usr/bin/env python

import os
import sys
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import read_fasta, write_fasta

# Concatenate FASTAs based on randomly sampled IDs.
# Any species not found in one FASTA will be excluded from the output.

def read_subsample_sets(subsample_dir):
    subsample_sets = dict()
    for file in os.listdir(subsample_dir):
        if file.endswith('.txt'):
            subsample_num = file.replace('.txt', '')
            with open(os.path.join(subsample_dir, file), 'r') as f:
                lines = f.readlines()
                lines = [line.strip() for line in lines]
                subsample_sets[subsample_num] = lines
    return subsample_sets

def prep_concat_fasta(subsample_sets, fasta_dir, fasta_suffix, out_dir):
    for subsample_num in subsample_sets.keys():
        set_genes = subsample_sets[subsample_num]
        set_seqs = {}
        for gene in set_genes:
            gene_fasta = os.path.join(fasta_dir, gene + fasta_suffix)
            if not os.path.exists(gene_fasta):
                sys.exit('Gene FASTA file not found: ' + gene_fasta)
            set_seqs[gene] = read_fasta(gene_fasta)

        concat_seq = defaultdict(str)
        tip_tallies = defaultdict(int)
        for gene in set_genes:
            gene_seqs = set_seqs[gene]
            for seq_id in gene_seqs.keys():
                tip_tallies[seq_id] += 1
                concat_seq[seq_id] += gene_seqs[seq_id]
        to_ignore = set()
        for seq_id in tip_tallies.keys():
            if tip_tallies[seq_id] != 12:
                to_ignore.add(seq_id)
        
        for seq_id in to_ignore:
            del concat_seq[seq_id]
        
        out_fasta = os.path.join(out_dir, subsample_num + '.fna')
        write_fasta(concat_seq, out_fasta)


drosophila_sets = read_subsample_sets('/home6/gmdougla/projects/aa_distance/paml_workflow/concat_seqs/genes_to_concat/Drosophila')
prep_concat_fasta(subsample_sets = drosophila_sets,
                  fasta_dir = '/home6/gmdougla/projects/aa_distance/paml_workflow/Drosophila/all_cds',
                  fasta_suffix = '.DNA.afa.mask',
                  out_dir = '/home6/gmdougla/projects/aa_distance/paml_workflow/concat_seqs/Drosophila')

mammal_sets = read_subsample_sets('/home6/gmdougla/projects/aa_distance/paml_workflow/concat_seqs/genes_to_concat/Mammal')
prep_concat_fasta(subsample_sets = mammal_sets,
                  fasta_dir = '/home6/gmdougla/db/orthomam_v12/omm_filtered_NT_CDS_unzipped/',
                  fasta_suffix = '_NT_AL.fasta',
                  out_dir = '/home6/gmdougla/projects/aa_distance/paml_workflow/concat_seqs/Mammal')

strep_sets = read_subsample_sets('/home6/gmdougla/projects/aa_distance/paml_workflow/concat_seqs/genes_to_concat/Strep')
prep_concat_fasta(subsample_sets = strep_sets,
                  fasta_dir = '/home6/gmdougla/projects/aa_distance/paml_workflow/Streptococcus_working/macse_v2_workflow/fasta_links',
                  fasta_suffix = '_final_mask_align_NT.aln',
                  out_dir = '/home6/gmdougla/projects/aa_distance/paml_workflow/concat_seqs/Strep')
