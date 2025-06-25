rm(list = ls(all.names = TRUE))

demask_raw <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/literature_data/demask.txt',
                         header=TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(demask_raw) <- colnames(demask_raw)

write.table(demask_raw,
            file = gzfile('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/similarities/demask.tsv.gz'),
            sep = '\t', row.names = TRUE, col.names = NA, quote = FALSE)
