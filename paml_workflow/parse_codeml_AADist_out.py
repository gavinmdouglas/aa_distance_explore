#!/usr/bin/env python

import argparse
import sys
import os


def main():
    parser = argparse.ArgumentParser(

        description=
        '''
        Parse PAML CodeML AADist output files to get information of interest about run.

        Note that rather than a general-purpose tool, this script assumes a certain folder structure (to parse out the AA dist definition used).
        ''',

        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input',
                        help='Path to CODEML results file.',
                        required=True,
                        metavar='PATH',
                        type=str)
    
    parser.add_argument('-H', '--header',
                        help='Print header line.',
                        action='store_true',
                        default=False)
    
    args = parser.parse_args()

    # Get directory of input file.
    input_dir = os.path.basename(os.path.dirname(os.path.abspath(args.input)))
    input_dir_split = input_dir.split(".")
    subsample_num = input_dir_split[0]
    distance_type = input_dir_split[1]

    with open(args.input, 'r') as input_fh:
        
        # Get first 6 lines, which should look like this:
        # CODONML (in paml version 4.10.7, June 2023)  ../../../concat_seqs/Drosophila/subsample9.fna
        # Model: One dN/dS ratio,
        # Codon frequency model: F3x4
        # grantham.dat, geometric
        # ns =   6  ls = 4751
        # 
        first6 = [next(input_fh) for _ in range(6)]

        if not first6[0].startswith("CODONML"):
            sys.exit("Error: First line of input file does not start with 'CODONML'.")
        elif not first6[1].startswith("Model: One dN/dS ratio,"):
            sys.exit("Error: Second line of input file does not start with 'Model: One dN/dS ratio,'.")
        elif not first6[2].startswith("Codon frequency model: F3x4"):
            sys.exit("Error: Third line of input file does not start with 'Codon frequency model: F3x4'.")
        elif not first6[3].startswith("grantham.dat, geometric"):
            sys.exit("Error: Fourth line of input file does not start with 'grantham.dat, geometric'.")
        elif not first6[4].startswith("ns ="):
            sys.exit("Error: Fifth line of input file does not start with 'ns ='.")
        
        ns_ls_linesplit = first6[4].split()
        if ns_ls_linesplit[0] == "ns" and ns_ls_linesplit[1] == "=" and ns_ls_linesplit[3] == "ls" and ns_ls_linesplit[4] == "=":
            ns_val = ns_ls_linesplit[2]
            ls_val = ns_ls_linesplit[5]

        tree_length_hits = []
        kappa_hits = []
        omega_hits = []
        b_hits = []
        a_hits = []
        tree_length_dN = []
        tree_length_dS = []
        lnL_hits = []
        np_hits = []

        for line in input_fh:
            linesplit = line.split()
            if len(linesplit) == 0:
                continue

            if line.startswith("tree length ="):
                tree_length_hits.append(linesplit[-1])
            elif line.startswith("kappa (ts/tv) ="):
                kappa_hits.append(linesplit[-1])
            elif line.startswith("omega (dN/dS) ="):
                omega_hits.append(linesplit[-1])
            elif line.startswith("b ="):
                b_hits.append(linesplit[-1])
            elif line.startswith("a ="):
                a_hits.append(linesplit[-1])
            elif line.startswith("tree length for dN:"):
                tree_length_dN.append(linesplit[-1])
            elif line.startswith("tree length for dS:"):
                tree_length_dS.append(linesplit[-1])
            elif line.startswith("lnL("):
                line = line.replace("ntime:", "ntime: ")
                line = line.replace("np:", "np: ")
                line = line.replace("):", "): ")
                linesplit=line.split()
                lnL_hits.append(linesplit[4])
                if linesplit[2] == "np:":
                    np_hits.append(linesplit[3].replace('):', ''))

        # Check that all hits are just one value.
        if len(tree_length_hits) != 1:
            print(args.input, file=sys.stderr)
            sys.exit("Error: Zero or more than one tree length hit found.")
        if len(kappa_hits) != 1:
            sys.exit("Error: Zero or more than one kappa hit found.")
        if len(omega_hits) != 1:
            sys.exit("Error: Zero or more than one omega hit found.")
        if len(b_hits) != 1:
            sys.exit("Error: Zero or more than one b hit found.")
        if len(a_hits) != 1:
            sys.exit("Error: Zero or more than one a hit found.")
        if len(tree_length_dN) != 1:
            sys.exit("Error: Zero or more than one tree length dN hit found.")
        if len(tree_length_dS) != 1:
            sys.exit("Error: Zero or more than one tree length dS hit found.")
        if len(lnL_hits) != 1:
            sys.exit("Error: Zero or more than one lnL hit found.")
        if len(np_hits) != 1:
            sys.exit("Error: Zero or more than one np hit found.")

        if args.header:
            print("subsample_num\tdistance_type\tlnL\tnum_param\ttree_length\tkappa\tomega\tb\ta\ttree_length_dN\ttree_length_dS\tns\tls")

        print(f"{subsample_num}\t{distance_type}\t{lnL_hits[0]}\t{np_hits[0]}\t{tree_length_hits[0]}\t{kappa_hits[0]}\t{omega_hits[0]}\t{b_hits[0]}\t{a_hits[0]}\t{tree_length_dN[0]}\t{tree_length_dS[0]}\t{ns_val}\t{ls_val}")


if __name__ == '__main__':
    main()
