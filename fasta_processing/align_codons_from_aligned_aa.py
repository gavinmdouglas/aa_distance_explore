#!/usr/bin/python3

import argparse
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import read_fasta, write_fasta, read_genetic_code, translate_DNA


def main():

    parser = argparse.ArgumentParser(
        description='''
        Given a FASTA containing an amino acid alignment,
        and a separate FASTA with unaligned codons for the same sequences,
        align the codons to match the amino acid alignment
        and output a new FASTA with this alignment.

        Note that this is a really basic approach compared to MACSE2 and others: 
        it assumes that the codon sequences should agree perfectly with the amino acid alignments
        (and so would fail if applied to gene data with frameshifts and/or premature stop codons).
        ''',

        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-a', '--aa',
                        metavar='PATH',
                        type=str,
                        help='Path to file with aligned amino acids.',
                        required=True)

    parser.add_argument('-n', '--nucl',
                        metavar='PATH',
                        type=str,
                        help='Path to file with aligned amino acids.',
                        required=True)

    parser.add_argument('-o', '--output',
                        metavar='OUTPUT',
                        type=str,
                        help='Path to output FASTA with aligned codons.',
                        required=True)

    parser.add_argument('-g', '--genetic_code',
                        metavar='PATH',
                        type=str,
                        help='Path to genetic code file.',
                        required=False,
                        default='/home6/gmdougla/db/codon_tables/genetic_codes/code11.tsv')

    args = parser.parse_args()

    aa_seqs = read_fasta(args.aa)
    nucl_seqs = read_fasta(args.nucl)

    if len(aa_seqs) != len(nucl_seqs):
        sys.exit('Number of sequences in FASTAs do not match.')

    aa_seq_len = len(list(aa_seqs.values())[0])
    for seq_id in aa_seqs:
        if seq_id not in nucl_seqs:
            sys.exit('Sequence IDs do not match between FASTAs.')
        if len(aa_seqs[seq_id]) != aa_seq_len:
            sys.exit('Sequences are not the same length in amino acid alignment.')

    codon_to_aa, possible_start_codons = read_genetic_code(args.genetic_code)

    aligned_codons = dict()
    for seq_id in aa_seqs:
        aligned_codons[seq_id] = ''
        aa_seq = aa_seqs[seq_id]
        nucl_seq = nucl_seqs[seq_id]
        nucl_pos = 0
        for aa_pos in range(len(aa_seq)):
            aa = aa_seq[aa_pos]
            if aa == '-':
                aligned_codons[seq_id] += '---'
                continue

            codon = nucl_seq[nucl_pos:nucl_pos + 3]
            codon_aa = translate_DNA(codon, codon_to_aa, possible_start_codons, check_start=False, silent=True)

            if codon_aa != aa:
                sys.exit(f'Error: amino acid and codon sequences do not match at position {aa_pos} for sequence {seq_id}. ' + codon_aa + ' != ' + aa)
            else:
                aligned_codons[seq_id] += codon
                nucl_pos += 3

        # Check if any nucleotides remain.
        if nucl_pos < len(nucl_seq):
            final_nucls = nucl_seq[nucl_pos:]
            if len(final_nucls) != 3:
                sys.exit('Error: final nucleotides do not form a single complete codon.')
            elif final_nucls in codon_to_aa and codon_to_aa[final_nucls] == '*':
                print('All aligned except for stop codon at end', file=sys.stderr)
            else:
                sys.exit('Error: final three nucleotides do not form a stop codon.')

    # Check that, when translated, the sequences match each aligned amino acid sequence.
    for seq_id in aligned_codons:
        translated_from_codons = translate_DNA(aligned_codons[seq_id], codon_to_aa, possible_start_codons, check_start=False, silent=True)
        if translated_from_codons != aa_seqs[seq_id]:
            sys.exit(f'Error: translated codons do not match amino acid sequence for sequence {seq_id}.')

    write_fasta(aligned_codons, args.output)


if __name__ == '__main__':
    main()
