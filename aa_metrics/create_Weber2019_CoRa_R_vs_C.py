import gzip
import sys

# Definitions of radical vs. conservative amino acids from CoRa GitHub.
partition=(1,2,3,3,3,2,2,1,2,4,4,2,4,4,1,3,3,2,2,1)
AAs='ARNDCQEGHILKMFPSTWYV'

# Create table of radical vs. conservative amino acids based on definition in CoRa paper (Weber et al. 2019).
def_file = '/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_RvC/RvC_polarity_volumeWeber2019.tsv.gz'
with gzip.open(def_file, 'wt') as def_rc_fh:
    def_rc_fh.write("aa1\taa2\tsub_type\n")
    for i in range(len(partition)):
        for j in range(len(partition)):
            if i == j:
                continue
        
            if partition[i] == partition[j]:
                def_rc_fh.write(f'{AAs[i]}\t{AAs[j]}\tConservative\n')
            else:
                def_rc_fh.write(f'{AAs[i]}\t{AAs[j]}\tRadical\n')
