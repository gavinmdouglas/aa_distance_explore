#!/usr/bin/env python

import sys
import os
import gzip

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from functions import read_fasta, write_fasta, translate_DNA, read_genetic_code

codon_to_aa, possible_start_codons = read_genetic_code('/home6/gmdougla/db/codon_tables/genetic_codes/code1.tsv')

orig_fasta = read_fasta('/home6/gmdougla/projects/aa_selection/human_variants/gencode/gencode.v45.pc_translations_rasp.fa.gz')

cds_fasta = {}
with gzip.open('/home6/gmdougla/projects/aa_selection/human_variants/gencode/Homo_sapiens.GRCh38.cds.all.fa.gz', 'rt') as cds_fh:
    for fasta_line in cds_fh:
        if fasta_line[0] == '>':
            fasta_linesplit = fasta_line.split()
            transcript = fasta_linesplit[0][1:]
            cds_fasta[transcript] = ''
        else:
            cds_fasta[transcript] += fasta_line.rstrip().upper()

prot_to_transcript = {}
with open('/home6/gmdougla/projects/aa_selection/human_variants/gencode/gencode_header_split.tsv', 'r') as map_fh:
    for line in map_fh:
        line = line.strip().split('\t')
        if line[0] in prot_to_transcript.keys():
            sys.exit('Protein ID already present...')
        prot_to_transcript[line[0]] = line[1]

prot_w_snv = set()
snv_files = os.listdir('/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/by_chrom')
for snv_file in snv_files:
    if not snv_file.endswith('.tsv.gz'):
        continue
    with gzip.open('/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/by_chrom/' + snv_file, 'rt') as snv_fh:
        header = snv_fh.readline()
        for line in snv_fh:
            line = line.strip().split('\t')
            prot_w_snv.add(line[4])

final_prot_fasta = {}
final_cds_fasta = {}
for prot_id in sorted(list(prot_w_snv)):
    
    transcript_id = prot_to_transcript[prot_id]

    if transcript_id not in cds_fasta.keys():
        print('Missing transcript ID for ' + prot_id)
        continue

    final_cds_fasta[prot_id] = cds_fasta[transcript_id]

    final_prot_fasta[prot_id] = orig_fasta[prot_id]

    # Sanity check on translated gene.
    translated_prot = translate_DNA(in_dna=final_cds_fasta[prot_id],
                                    codon_to_aa=codon_to_aa,
                                    possible_start_codons=possible_start_codons,
                                    check_start=False,
                                    silent=True,
                                    ambig_char='X')

    # If last residue is stop codon, then remove, and also remove last three nucleotides from CDS.
    if translated_prot[-1] == '*' or (len(translated_prot) == len(final_prot_fasta[prot_id]) + 1 and translated_prot[-1] == 'X'):
        cds_trim_length = -3
        if len(cds_fasta[transcript_id]) % 3 != 0:
            cds_trim_length = (len(cds_fasta[transcript_id]) % 3) * -1
        translated_prot = translated_prot[:-1]
        final_cds_fasta[prot_id] = final_cds_fasta[prot_id][:cds_trim_length]
        translated_prot = translate_DNA(in_dna=final_cds_fasta[prot_id],
                                    codon_to_aa=codon_to_aa,
                                    possible_start_codons=possible_start_codons,
                                    check_start=False,
                                    silent=True,
                                    ambig_char='X')
    
    translated_prot = translated_prot.replace('*', 'X').replace('U', 'X')
    tmp_fasta = final_prot_fasta[prot_id]
    tmp_fasta = tmp_fasta.replace('*', 'X').replace('U', 'X')

    if translated_prot[1:] != tmp_fasta[1:]:
        print(transcript_id, prot_id)
        print(cds_fasta[transcript_id])
        print(translated_prot)
        print(final_prot_fasta[prot_id])
        sys.exit('ERROR: Translated protein does not match original protein')

write_fasta(final_prot_fasta, '/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/gencode_w_snv.faa')
write_fasta(final_cds_fasta, '/home6/gmdougla/projects/aa_selection/human_variants/nonsyn_snvs/gencode_w_snv.fna')
