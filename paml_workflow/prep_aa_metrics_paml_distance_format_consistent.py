import os
import gzip

outdir = '/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/paml_format_dist_consistent'

aa_order = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
            'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

short_to_tri = {}
with gzip.open('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/aa_map.tsv.gz', 'rt') as aa_map_fh:
    aa_map_fh.readline()
    for line in aa_map_fh:
        line = line.strip().split('\t')
        short_to_tri[line[2]] = line[1]

# Needed to uncomment the folders to run (and ran one at a time, as needed)
# sim_folder = "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity_consistent"
# sim_folder = "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_RvC/prepped_similarity_consistent"
sim_folder = "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/DISTATIS_working/prepped_similarity_consistent"

sim_files = os.listdir(sim_folder)
sim_files = [f for f in sim_files if f.endswith('.tsv.gz')]

for sim_file in sim_files:
    dist = {}
    with gzip.open(os.path.join(sim_folder, sim_file), 'rt') as sim_fh:
            sim_fh.readline()
            for line in sim_fh:
                line = line.strip().split('\t')
                aa1 = line[0]
                aa2 = line[1]
                dist_val = 1.0 - float(line[2])
                dist[(aa1, aa2)] = dist_val

    # Now write out the distance matrix
    with open(os.path.join(outdir, sim_file.replace('.tsv.gz', '.txt')), 'w') as out_fh:
        for aa1 in aa_order:
            outline = []
            for aa2 in aa_order:
                if aa1 == aa2:
                    if len(outline) > 0:
                        out_fh.write(' '.join(outline) + '\n')
                        break
                else:
                    outline.append(f'{dist[(aa1, aa2)]:.4f}')

        out_fh.write('\n')
        out_fh.write(' '.join(aa_order) + '\n')
        aa_order_tri = [short_to_tri[aa] for aa in aa_order]
        out_fh.write(' '.join(aa_order_tri) + '\n')
        out_fh.write('\n')
        out_fh.write('\n')
