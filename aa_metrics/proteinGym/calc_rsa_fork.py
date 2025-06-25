#!/usr/bin/python

# Original calc_rsa.py portions:
# Copyright (c) 2018 Benjamin R. Jack
# Licensed under MIT License - see https://opensource.org/license/mit

import argparse
import csv

# This block was taken from calc_rsa.py by Benjamin R. Jack (MIT license):
# ASA normalization constants were taken from:
# M. Z. Tien, A. G. Meyer, D. K. Sydykova, S. J. Spielman, C. O. Wilke (2013).
# Maximum allowed solvent accessibilities of residues in proteins. PLOS ONE
# 8:e80635.
RES_MAX_ACC = {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
                'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0, \
                'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0, \
                'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0, \
                'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0}

# Similarly, parse_dssp_line and parse_dssp functions were taken from calc_rsa.py by Benjamin R. Jack (MIT license).
def parse_dssp_line(line):
    '''
    Extract values from a single line of DSSP output and calculate RSA.
    '''
    solvent_acc = line[35:39].strip()  # record SA value for given AA
    amino_acid = line[13].strip()  # retrieve amino acid
    residue = line[6:10].strip()
    chain = line[11].strip()
    secondary_structure = line[16].strip() # secondary structure
    if amino_acid.islower():
        # if the AA is a lowercase letter, then it's a cysteine
        amino_acid = "C"
    if amino_acid in RES_MAX_ACC:
        max_acc = RES_MAX_ACC[amino_acid]  # Find max SA for residue
        rsa = float(solvent_acc) / max_acc # Normalize SA
    else:
        rsa = None
    return {'pdb_position': residue,
            'pdb_aa': amino_acid,
            'chain': chain,
            'rsa': rsa,
            'structure': secondary_structure}


def parse_dssp(raw_dssp_output):
    '''
    Parse a DSSP output file and return a dictionary of amino acids, PDB residue
    number, chain, and RSA.
    '''
    with open(raw_dssp_output, 'r') as dssp_file:
        lines = dssp_file.readlines()
        output = []  # list of dictionaries for output
    for line in lines[28:]:
        # Skip first 28 lines of header info
        output_line = parse_dssp_line(line)
        if output_line['rsa'] is not None:
            # Skip lines with no RSA value
            output.append(output_line)

    return output


def main():
    parser = argparse.ArgumentParser(

        description=
        '''
        Parse mkdssp output to get relative solvent accessibility (RSA) and other related info.

        This is a modification of calc_rsa.py by Benjamin R. Jack (MIT license): https://github.com/clauswilke/proteinER/blob/master/src/calc_rsa.py

        Difference is that the original code automatically ran mkdssp on input PDBs, while I prefer to run mkdssp separately and then parse the output with Python.

        The output is a CSV file with columns for PDB position, chain, secondary structure, PDB amino acid, and RSA value. Here is the description taken from calc_rsa.py:
            
            Column name   Description
            ===================================================================
            pdb_position  Residue number, extracted from the input PDB file.

            chain         PDB chain.

            structure     A letter indicating the secondary structure assigned 
                          to this residue. This column may be empty of no 
                          secondary structure can be assigned. Non-empty values 
                          are:

                          Code      Description
                          -----------------------------------------------------
                          H         Alpha Helix
                          B         Beta Bridge
                          E         Strand
                          G         Helix-3
                          I         Helix-5
                          T         Turn
                          S         Bend 
            
            pdb_aa        Single-letter amino acid.
            
            rsa           Relative solvent accessibility, normalized to the     
                          maximum possible solvent accessibility.
        ''',

        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input_dssp',
                        help='Path to input output (should be output of mkdssp).',
                        required=True,
                        metavar='INPUT',
                        type=str)
    
    parser.add_argument('-o', '--output_csv',
                        help='Path to output CSV file.',
                        required=True,
                        metavar='OUTPUT_CSV',
                        type=str)

    args = parser.parse_args()

    output_dict = parse_dssp(args.input_dssp)
    with open(args.output_csv, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['pdb_position', 
                                                     'chain', 'structure',
                                                     'pdb_aa', 'rsa'])
        writer.writeheader()
        writer.writerows(output_dict)


if __name__ == "__main__":
    main()
