#!/usr/bin/env python

from Bio import SeqIO
import argparse

def main():
    parser = argparse.ArgumentParser(description='Rename sequences in a FASTA file to seqX format, and get mapfile of original to new IDs.')
    parser.add_argument('input_file', help='Path to input FASTA file')
    parser.add_argument('output_file', help='Path to output FASTA file')
    parser.add_argument('--map_file', default='id_mapping.tsv', 
                        help='Path to save ID mapping (default: id_mapping.tsv)')
    
    args = parser.parse_args()
    
    sequences = list(SeqIO.parse(args.input_file, "fasta"))
    id_mappings = []
    for i, record in enumerate(sequences):
        original_id = record.id
        new_id = f"seq{i + 1}"
        record.id = new_id
        record.description = new_id
        id_mappings.append((original_id, new_id))
    SeqIO.write(sequences, args.output_file, "fasta")
    
    with open(args.map_file, 'w') as f:
        f.write("Original_ID\tNew_ID\n")
        for original_id, new_id in id_mappings:
            f.write(f"{original_id}\t{new_id}\n")

if __name__ == "__main__":
    main()