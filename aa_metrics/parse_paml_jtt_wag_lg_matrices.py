import pandas as pd
import os
import sys
import gzip

# Parse the matrices from PAML's JTT, WAG, and LG files and write them to TSV format.

tri_to_single = {}
with gzip.open("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/aa_map.tsv.gz", 'rt') as f:
    f.readline()
    for line in f:
        line = line.strip().split('\t')
        tri_to_single[line[1]] = line[2]

datfiles = os.listdir("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/PAML_aa_matrices")

for datfile in datfiles:
    if not datfile.endswith('.dat'):
        continue

    datname = datfile.split('.')[0]

    # Rename jones to jtt, for consistency.
    if datname == "jones":
        datname = "jtt"

    datpath = os.path.join("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/PAML_aa_matrices", datfile)

    # Read through once to get AA order.
    aa_tri_order = None
    with open(datpath) as f:
        for line in f:
            line = line.strip().split()
            if not line:
                continue
            if line[0] == 'Ala':
                aa_tri_order = line
                break
    if aa_tri_order is None:
        sys.exit(f"Error: Could not find AA order in {datfile}")

    aa_order = [tri_to_single[tri] for tri in aa_tri_order]

    # Then read through again to fill in matrix values.
    df = pd.DataFrame(index=aa_order, columns=aa_order, dtype=float)
    with open(datpath) as f:
        current_row = 1
        for line in f:
            line = line.strip().split()
            if len(line) == 0:
                continue
            elif len(line) != current_row:
                sys.exit('Error: Unexpected line length in matrix values')
            
            for i in range(len(line)):
                df.iloc[current_row, i] = float(line[i])
                df.iloc[i, current_row] = float(line[i])
            current_row += 1

            if current_row == 20:
                break
    
    df.to_csv(f'/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/sub_matrices/{datname}.tsv.gz', sep='\t', compression='gzip', na_rep='NA')
