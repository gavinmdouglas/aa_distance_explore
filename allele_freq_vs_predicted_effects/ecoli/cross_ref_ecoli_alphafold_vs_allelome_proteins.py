import os
import gzip
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

# Check that the protein sequences in the AlphaFold database for E. coli can be mapped 1:1 to the allelome proteins (which have Blattner numbers).

# This ID mapping file was downloaded (April 30, 2026) from:
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/ECOLI_83333_idmapping.dat.gz

uniprot_to_bnumber = {}
bnumber_to_uniprot = {}

uniprot_to_ignore = set()
bnumber_to_ignore = set()

with gzip.open('/Users/gavin/Drive/research/aa_distance/data/ECOLI_83333_idmapping.dat.gz', 'rt') as fh:
    for line in fh:
        uniprot, db, identifier = line.strip().split('\t')
        if db == 'Gene_OrderedLocusName' and identifier.startswith('b'):
            if uniprot in uniprot_to_bnumber:
                print(f'Multiple b-numbers for {uniprot}: {uniprot_to_bnumber[uniprot]} and {identifier}', file=sys.stderr)
                uniprot_to_ignore.add(uniprot)
            if identifier in bnumber_to_uniprot:
                print(f'Multiple UniProt accessions for {identifier}: {bnumber_to_uniprot[identifier]} and {uniprot}', file=sys.stderr)
                bnumber_to_ignore.add(identifier)
            uniprot_to_bnumber[uniprot] = identifier
            bnumber_to_uniprot[identifier] = uniprot

# Read in sequences linked to AlphaFold output.
af_seqs_raw = SeqIO.to_dict(SeqIO.parse('/Users/gavin/Desktop/alphafold_working/UP000000625_83333_ECOLI_v4.faa', 'fasta'))
af_seqs = {}
for seqid in af_seqs_raw.keys():
    clean_id = seqid.split('-')[1]
    af_seqs[clean_id] = af_seqs_raw[seqid].upper()
del af_seqs_raw

# Read in allelome sequences.
allelome_dir = '/Users/gavin/Drive/research/aa_distance/data/allelome_data'
allelome_seqs = {}
for filename in os.listdir(allelome_dir):
    if not filename.endswith('.csv.gz'):
        continue
    bnumber = filename.replace('_dfz.csv.gz', '')
    filepath = os.path.join(allelome_dir, filename)
    in_df = pd.read_csv(filepath, compression='gzip')
    aaseq = ''.join(in_df['consensus_AAseq'].tolist()).replace('-', '').replace('*', '').upper()
    allelome_seqs[bnumber] = aaseq

print(f'Number of sequences in AlphaFold dataset: {len(af_seqs)}')
print(f'Number of sequences in allelome dataset: {len(allelome_seqs)}')

matched_ids = set()
for bnumber in sorted(allelome_seqs.keys()):
    if bnumber in bnumber_to_ignore:
        continue
    if bnumber not in bnumber_to_uniprot:
        print(f'No UniProt accession found for {bnumber}')
        continue
    uniprot = bnumber_to_uniprot[bnumber]
    if uniprot in uniprot_to_ignore:
        continue
    if uniprot not in af_seqs:
        print(f'No AlphaFold sequence found for {uniprot} (b-number {bnumber})')
        continue
    af_seq = str(af_seqs[uniprot].seq)
    allelome_seq = allelome_seqs[bnumber]

    # Ignore start codon, as this will be ignored anyway.
    af_seq_trim = af_seq[1:]
    allelome_seq_trim = allelome_seq[1:]

    # And then comapre all residues, allowing for some differences (allowing for 1 difference).
    if len(af_seq_trim) == len(allelome_seq_trim):
        diffs = sum(a != b for a, b in zip(af_seq_trim, allelome_seq_trim))
        if diffs <= 1:
            matched_ids.add(bnumber)

print(f'Number of matching sequences: {len(matched_ids)}')

records = []
for bnumber in sorted(matched_ids):
    uniprot = bnumber_to_uniprot[bnumber]
    rec = SeqRecord(Seq(allelome_seqs[bnumber]), id=bnumber + '_' + uniprot, description='')
    records.append(rec)
SeqIO.write(records, '/Users/gavin/Drive/research/aa_distance/data/Ecoli_focal_seqs.faa', 'fasta')
