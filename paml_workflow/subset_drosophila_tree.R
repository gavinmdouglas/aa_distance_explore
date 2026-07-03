rm(list = ls(all.names = TRUE))

# Parse Drosophila tree.

library(ape)

full_tree <- ape::read.tree("/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/Drosophila/Thomas_Hahn_2017/drosophila-25spec-tree.tre")

long_to_short <- list()
for (tipname in full_tree$tip.label) {
  tmp_split <- strsplit(tipname, split = "\\.")[[1]]
  species_short <- substring(tmp_split[2], 1, 3)
  long_to_short[[tipname]] <- paste(tmp_split[1], species_short, sep = '')
}

for (long in names(long_to_short)) {
  short <- long_to_short[[long]]
  full_tree$tip.label[full_tree$tip.label == long] <- short
}

focal_short <- c('Dmel', 'Dsim', 'Dsec', 'Dyak', 'Dere', 'Dana')
other_tips <- setdiff(full_tree$tip.label, focal_short)

focal_tree <- ape::drop.tip(phy = full_tree, tip = other_tips, trim.internal = TRUE)

write.tree(phy = focal_tree, file = "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/Drosophila/focal_melgroup.tre")

# Unroot tree:
focal_tree <- ape::unroot(focal_tree)
write.tree(phy = focal_tree, file = "/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/PAML_workflow/Drosophila/focal_melgroup.unroot.tre")

