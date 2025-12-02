import os
import gzip

# Sanity check that the PAML models fit worse based on the complement of the distances.
# This script simply takes the complement of the prepped matrices and writes them out.

prepped_paml_dir = '/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/paml_format_dist_consistent'
outfolder = '/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/paml_format_dist_consistent_inverted_test'

prepped_files = os.listdir(prepped_paml_dir)
prepped_files = [f for f in prepped_files if f.endswith('.txt')]

for infile in prepped_files:

    with open(os.path.join(outfolder, infile), 'wt') as outfh:

        with open(os.path.join(prepped_paml_dir, infile), 'rt') as f:
            for line in f:
                if line.startswith('0'):
                    line_split = line.strip().split()
                    inverted_values = [str(1.0 - float(x)) for x in line_split]
                    outfh.write(' '.join(inverted_values) + '\n')
                else:
                    outfh.write(line)
                    continue
