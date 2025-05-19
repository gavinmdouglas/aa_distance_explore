#!/usr/bin/env python

from Bio import SeqIO
import argparse
import os
import sys

def convert_fasta_to_phylip(input_file, output_prefix):
    """
    Convert an alignment from FASTA format to PHYLIP format,
    renaming all sequences to "seqX" format to ensure unique IDs.
    Will write a map file of original to new IDs.
    
    Args:
        input_file (str): Path to input FASTA file
        output_prefix (str): Prefix for output files.
    """
    sequences = list(SeqIO.parse(input_file, "fasta"))
    
    # id_mappings = []
    # for i, record in enumerate(sequences):
    #     original_id = record.id
    #     new_id = f"seq{i + 1}"
    #     record.id = new_id
    #     record.description = new_id
    #     id_mappings.append((original_id, new_id))

    # with open(output_prefix + "_id_map.tsv", 'w') as f:
    #     f.write("Original_ID\tNew_ID\n")
    #     for original_id, new_id in id_mappings:
    #         f.write(f"{original_id}\t{new_id}\n")

    # Sanity check that all sequences are same length.
    seq_length = len(sequences[0].seq)
    for record in sequences:
        if len(record.seq) != seq_length:
            raise ValueError("All sequences must be of the same length.")

    with open(output_prefix + ".phylip", 'w') as out_fh:
        out_fh.write(f"     {len(sequences)} {seq_length}\n\n")
        for record in sequences:
            # Write in sequential format, with a space every 10 characters.
            seq_str = str(record.seq)

            # id_field = record.id
            # if len(id_field) > 10:
            #     sys.exit(f"ID '{id_field}' is too long for PHYLIP format. Please shorten it.")
            # elif len(id_field) < 10:
            #     id_field = id_field.ljust(10)

            out_fh.write(f"{record.id}                 {seq_str}\n")
        
        out_fh.write("\n// end of file\n")


def main():
    parser = argparse.ArgumentParser(description='Convert FASTA to PHYLIP format with renamed sequences.')
    parser.add_argument('input_file', help='Path to input FASTA file')
    parser.add_argument('output_prefix', help='Prefix for output files')
    
    args = parser.parse_args()
    convert_fasta_to_phylip(args.input_file, args.output_prefix)

if __name__ == "__main__":
    main()