#!/usr/bin/env python

import argparse
import sys
import os
import gzip
import csv

# Read in functions from script one directory up.
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import read_fasta, write_fasta

panta_dir = '/home6/gmdougla/projects/aa_distance/paml_workflow/Streptococcus_working/annotations/Streptococcus_panta/'
out_dir = '/home6/gmdougla/projects/aa_distance/paml_workflow/Streptococcus_working/annotations/Streptococcus_freq_genes/'
min_genomes = 68

genome_to_seqs = {}
samples_dir = os.path.join(panta_dir, 'samples')

# Loop through all genome IDs in sample_dir.
for genome_id in os.listdir(samples_dir):
    genome_dir = os.path.join(samples_dir, genome_id)
    if not os.path.isdir(genome_dir):
        sys.exit('Expected directory for genome ID: ' + genome_id + ' but found a file instead. Exiting.')
    for fasta_file in os.listdir(genome_dir):
        if fasta_file.endswith('.fna'):
            if genome_id in genome_to_seqs:
                sys.exit('Found multiple FASTA files for genome ID: ' + genome_id + ' Exiting.')
            fasta_path = os.path.join(genome_dir, fasta_file)
            genome_to_seqs[genome_id] = read_fasta(fasta_path, cut_header=True)
    if genome_id not in genome_to_seqs:
        sys.exit('Did not find FASTA file for genome ID: ' + genome_id + ' Exiting.')

with gzip.open(os.path.join(panta_dir, 'gene_presence_absence.csv.gz'), 'rt') as panta_fh:

    csv_reader = csv.reader(panta_fh, quotechar='"', delimiter=',')
    
    headerline = next(csv_reader)
    genome_ids = headerline[8:]

    for line_split in csv_reader:
        gene_seqs = {}
        gene = line_split[0]
        num_genomes = int(line_split[2])
        num_sequences = int(line_split[3])
        mean_seq_per_isolate = float(line_split[4])

        if mean_seq_per_isolate != 1.0 or num_genomes != num_sequences:
            continue
        elif num_genomes < min_genomes:
            continue

        genes = line_split[8:]
        for i in range(len(genome_ids)):
            gene_id = line_split[i + 8]
            if gene_id == '':
                continue
            elif ';' in gene_id or ',' in gene_id:
                sys.exit('Found multiple gene IDs for genome ID: ' + genome_ids[i] + ' and gene: ' + gene + ' Exiting.')
            genome_id = genome_ids[i]
            gene_seqs[genome_id] = genome_to_seqs[genome_id][gene_id]

        write_fasta(gene_seqs, os.path.join(out_dir, 'raw', gene + '.fna'))
