#!/usr/bin/env python

import argparse
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import thermoMPNN_dir_to_prefs


def main():
    parser = argparse.ArgumentParser(

        description='''
            Parse directory of ThermoMPNN output files to amino acid preference tables.
            ''',

        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input',
                        help='Input folder of ThermoMPNN output CSVs (gzipped).',
                        required=True,
                        metavar='DIR',
                        type=str)

    parser.add_argument('-o', '--out_folder',
                        help='Output folder.',
                        required=True,
                        metavar='PATH',
                        type=str)

    parser.add_argument('-a', '--aa_ref',
                        help='Optional table to create of mappings from sites to reference AA. Useful for sanity checks/troubleshooting.',
                        required=False,
                        default=None,
                        metavar='PATH',
                        type=str)

    parser.add_argument('--all_sites_not_present',
                        help='Flag to indicate that not all sites from min to max pos need to be present.',
                        action='store_true')

    parser.add_argument('--ref_fill',
                        help='Value to fill in for reference AA. Must be one of "one", "NA", "max", or "median". Default is "max"/',
                        required=False,
                        default='max',
                        metavar='STR',
                        type=str)
    
    parser.add_argument('--no_sum_scale',
                         help='Flag to indicate that preference values should not be scaled to sum to 1 per site.',
                         action='store_true')

    args = parser.parse_args()

    thermoMPNN_dir_to_prefs(in_dir=args.input,
                        out_dir=args.out_folder,
                        all_sites_present=not args.all_sites_not_present,
                        ref_fill=args.ref_fill,
                        sum_scale=not args.no_sum_scale)


if __name__ == '__main__':
    main()
