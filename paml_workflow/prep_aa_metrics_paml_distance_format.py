import os
import sys
import gzip

# Get all AA metrics
# (mainly distances, but also some similarities and sub matrices)
# into distance matrix format required by PAML.

outdir = '/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/paml_format_dist'

aa_order = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
            'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

short_to_tri = {}
with gzip.open('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/aa_map.tsv.gz', 'rt') as aa_map_fh:
    aa_map_fh.readline()
    for line in aa_map_fh:
        line = line.strip().split('\t')
        short_to_tri[line[2]] = line[1]

dist_indir = '/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/distances'
dist_files = os.listdir(dist_indir)
for dist_file in dist_files:
    if dist_file[0] == '.':
        continue
    dist = {}
    with gzip.open(os.path.join(dist_indir, dist_file), 'rt') as dist_fh:
        aa_names = dist_fh.readline().strip().split('\t')
        if aa_names[0] == 'aa' or aa_names[0] == 'AA' or aa_names[0] == 'amino_acid':
            aa_names = aa_names[1:]
        if len(aa_names) != 20:
            raise ValueError(f'Expected 20 AA names, got {len(aa_names)}')
        for aa in aa_names:
            if aa not in aa_order:
                raise ValueError(f'AA {aa} not in AA order')

            for line in dist_fh:
                line = line.strip().split('\t')
                if len(line) != 21:
                    raise ValueError(f'Expected 21 columns, got {len(line)}')
                aa = line[0]
                if aa not in aa_order:
                    raise ValueError(f'AA {aa} not in AA order')
                vals = line[1:]
                if len(vals) != 20:
                    raise ValueError(f'Expected 20 values, got {len(vals)}')
                for i, val in enumerate(vals):
                    dist[aa, aa_names[i]] = float(val)

    # Now write out the distance matrix
    with open(os.path.join(outdir, dist_file.replace('.tsv.gz', '.txt')), 'w') as out_fh:
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

# And then had to prep a few other metrics a little differently, as they are not in pairwise distance format already.
# These include the substitution matrices BLOSUM62 and VTML200, as well as the similarity matrices csw, empar, and ex.
# I has previously transformed these into basic similarity format (ranging from 0 to 1), so can just take the inverse.

sim_folder = "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/prepped_similarity"
#sim_files = ["blosum62.prob.tsv.gz", "csw.maxscaled.tsv.gz", "empar.maxscaled.tsv.gz", "ex.maxscaled.tsv.gz", "vtml200.prob.tsv.gz",
#             "proteingym_rescale1.tsv.gz", "proteingym_rescale2.tsv.gz", "proteingym_ranked.tsv.gz", "demask.tsv.gz"]
#sim_files = ["proteingym_mean_standard.tsv.gz", "proteingym_median_robust.tsv.gz", "proteingym_median_robust_orderedNorm.tsv.gz"]

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
                dist[(aa1, aa2)] = 1.0 - float(line[2])
    
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
