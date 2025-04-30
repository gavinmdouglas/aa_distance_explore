#!/usr/bin/env python

import argparse
import sys
import os
import gzip

def read_fasta(filename, cut_header=False, convert_upper=True):
    seq = {}

    def parse_fasta_lines(file_handle):
        for line in file_handle:
            line = line.rstrip()
            if len(line) == 0:
                continue
            if line[0] == ">":
                if cut_header:
                    name = line.split()[0][1:]
                else:
                    name = line[1:]
                name = name.strip()
                if name in seq:
                    raise ValueError("Stopping due to duplicated ID in file: " + name)
                seq[name] = ''
            else:
                line = line.strip()
                if convert_upper:
                    line = line.upper()
                seq[name] += line

    if filename[-3:] == ".gz":
        with gzip.open(filename, "rt") as fasta_in:
            parse_fasta_lines(fasta_in)
    else:
        with open(filename, "r") as fasta_in:
            parse_fasta_lines(fasta_in)

    return seq


def main():
    parser = argparse.ArgumentParser(

        description=
        '''
        Summarizes proportion of codons in alignment that are variable (split by syn and non-syn), and the proportion of 
        codons that contain gaps or ambiguous codons.
        ''',

        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input',
                        help='A file of all FASTA paths.',
                        required=True,
                        metavar='PATH',
                        type=str)
    
    parser.add_argument('-g', '--genetic_code',
                        help='Integer specifying the genetic code to use (as defined by NCBI). Must be 1 or 11.',
                        required=True,
                        metavar='INT',
                        type=int,
                        choices=[1, 11])

    args = parser.parse_args()

    code1_codon_to_aa = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
    code1_possible_start_codons = {'TTG', 'CTG', 'ATG'}

    code11_codon_to_aa = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
    code11_possible_start_codons = {'TTG', 'ATA', 'CTG', 'GTG', 'ATG', 'ATC', 'ATT'}

    if args.genetic_code == 1:
        codon_to_aa = code1_codon_to_aa
        possible_start_codons = code1_possible_start_codons
    elif args.genetic_code == 11:
        codon_to_aa = code11_codon_to_aa
        possible_start_codons = code11_possible_start_codons

    columns = ['name', 'num_seqs', 'alignment_length_codons', 'num_codonsites_w_gap_or_ambig', 'codonsites_w_syn', 'codonsites_w_nonsyn', 'codonsites_invariant_nogap_or_ambig']

    print('\t'.join(columns))

    with open(args.input, 'r') as input_fh:
        for line in input_fh:
            fasta_file = line.strip()
            basename = os.path.basename(fasta_file)
            basename = os.path.splitext(basename)[0]
        seqs = read_fasta(fasta_file)

        num_seqs = len(seqs)
        alignment_length = len(seqs[list(seqs.keys())[0]])
        if alignment_length % 3 != 0:
            sys.exit('Alignment length is not divisible by three for ' + basename + ' Exiting.')

        for seq in seqs:
            if len(seqs[seq]) != alignment_length:
                sys.exit('Sequences are not all the same length for ' + basename + ' Exiting.')

        codon_alignment_length = alignment_length // 3

        num_variable_codon_sites_same_aa_aa = 0
        num_variable_codon_sites_diff_aa = 0
        num_codon_sites_with_gaps_or_ambig = 0
        num_codon_sites_invariant_nogap_or_ambig = 0

        # Identify stop codons.
        stop_codons = set()
        for codon in codon_to_aa:
            if codon_to_aa[codon] == '*':
                stop_codons.add(codon)

        for i in range(alignment_length // 3):
            codons = [seq[i * 3:i * 3 + 3] for seq in seqs.values()]
            allowable_codons = set()

            for codon in codons:
                if '-' in codon or codon not in codon_to_aa:
                    num_codon_sites_with_gaps_or_ambig += 1
                    break
                if codon not in stop_codons and codon in codon_to_aa:
                    allowable_codons.add(codon)
            if len(allowable_codons) > 1:
                aas = [codon_to_aa[codon] for codon in allowable_codons]
                if len(set(aas)) > 1:
                    num_variable_codon_sites_diff_aa += 1
                else:
                    num_variable_codon_sites_same_aa_aa += 1
            elif len(allowable_codons) == 1:
                num_codon_sites_invariant_nogap_or_ambig += 1

        outline = [basename, num_seqs, codon_alignment_length,
                   num_codon_sites_with_gaps_or_ambig,
                   num_variable_codon_sites_same_aa_aa,
                   num_variable_codon_sites_diff_aa,
                   num_codon_sites_invariant_nogap_or_ambig]
        
        outline = [str(x) for x in outline]
        print('\t'.join(outline))

if __name__ == '__main__':
    main()
