#!/usr/bin/env python

import argparse
import sys
import os
import pysam

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from functions import read_fasta

def main():
    parser = argparse.ArgumentParser(

        description=
        '''
        Parse gnomAD VCF to get all non-synonymous mutations.

        For each such mutation, return:
            - Name of protein
            - Position in protein
            - The reference sequence AA
            - The non-ref. AA
            - Non-ref. AA count
            - Number of samples covered at site
        ''',

        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-v', '--vcf',
                        help='Path to bgzipped and tabix indexed VCF.',
                        required=True,
                        metavar='VCF',
                        type=str)

    parser.add_argument('-f', '--fasta',
                        help='Path to FASTA file of proteins to consider.',
                        required=True,
                        metavar='FASTA',
                        type=str)

    parser.add_argument('-c', '--sample_coverage',
                        help='Minimum sample coverage to consider a site.',
                        required=False,
                        metavar='MIN_COV',
                        type=int,
                        default=100000)

    args = parser.parse_args()

    # Input file checks.
    if not os.path.exists(args.vcf):
        sys.exit('Error: VCF file does not exist.')
    elif not args.vcf.endswith('.vcf.bgz'):
        sys.exit('Error: VCF file must end in ".vcf.bgz".')
    elif not os.path.exists(args.vcf + '.tbi'):
        sys.exit('Error: VCF file must be tabix indexed (i.e., a corresponding .tbi must be present).')

    vcf_in = pysam.VariantFile(args.vcf)

    if not os.path.exists(args.fasta):
        sys.exit('Error: FASTA file does not exist.')
    
    fasta = read_fasta(args.fasta)
    protein_set = set(list(fasta.keys()))

    # Parse header to get indices of important information in VEP field.
    header = str(vcf_in.header)
    vep_index_map = {}
    for line in header.split('\n'):
        if line.startswith('##INFO=<ID=vep'):
            vep_fields = line.split('Format: ')[1].split('|')
            for i, field in enumerate(vep_fields):
                vep_index_map[field] = i

    if len(vep_index_map) == 0:
        sys.exit('Error: Could not find vep field in header.')

    # Print header for output.
    print('\t'.join(['chrom', 'genome_pos', 'ref_base', 'alt_base', 'protein_id', 'protein_pos', 'ref_aa', 'alt_aa', 'variant_count', 'sample_count']))

    for record in vcf_in:
        # Skip if below sample coverage threshold, if not a point mutation, (or if there are no alternative allele counts, which is weird but common).
        if record.info['AN'] < args.sample_coverage or len(record.ref) > 1 or len(record.alts[0]) > 1 or record.info['AC'][0] == 0:
            continue

        # Throw error if it is multi-allelic (as this script expects one alternative allele per line).
        # (Although the same site over multiple sites can occur)
        if len(record.alts) > 1:
            sys.exit('Error: Multi-allelic site found. This script expects one alternative allele per line.')

        chrom = record.chrom
        genome_pos = str(record.pos)
        ref_base = record.ref
        alt_base = record.alts[0]
        variant_count = str(record.info['AC'][0])
        sample_count = str(record.info['AN'])

        # Parse vep info to see if non-synonymous (missense).
        # Keep only one effect per variant (first one in a possible protein parsed).
        variant_effect = record.info['vep']
        variant_effect_cat = '\t'.join(variant_effect)
        if "missense_variant" in variant_effect_cat:
            print(record)
            sys.exit()
            for effect in variant_effect:
                effect = effect.split('|')
                if len(effect) < len(vep_index_map):
                    continue
                if effect[vep_index_map['Consequence']] != 'missense_variant' or effect[vep_index_map['VARIANT_CLASS']] != 'SNV' or effect[vep_index_map['BIOTYPE']] != 'protein_coding':
                    continue
                protein_id = effect[vep_index_map['HGVSp']].split(':')[0]
                if protein_id in protein_set:
                    aa_set = effect[vep_index_map['Amino_acids']].split('/')
                    ref_aa = aa_set[0]
                    alt_aa = aa_set[1]
                    protein_pos = effect[vep_index_map['Protein_position']]

                    # Run sanity check on protein sequence.
                    if fasta[protein_id][int(protein_pos) - 1] != ref_aa:
                        sys.exit('Error: Protein sequence does not match reference amino acid: ' + protein_id + ' ' + protein_pos + ' ' + ref_aa + ' ' + fasta[protein_id][int(protein_pos) - 1] + "\n" + fasta[protein_id])

                    print('\t'.join([chrom, genome_pos, ref_base, alt_base, protein_id, protein_pos, ref_aa, alt_aa, variant_count, sample_count]))
                    break

if __name__ == '__main__':
    main()
