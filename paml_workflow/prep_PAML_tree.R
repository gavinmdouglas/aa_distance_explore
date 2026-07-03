library(ape)

write.tree.paml <- function(tree_file, file, ntrees = 1) {
  tree <- read.tree(tree_file)
  ntips <- length(tree$tip.label)
  con <- file(file, "w")
  cat(ntips, ntrees, "\n", file = con)
  cat(write.tree(tree), file = con)
  close(con)
  cat("Tree written to", file, "in PAML format\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript prep_PAML_tree.R <tree_file> <output_file>")
}

write.tree.paml(args[1], args[2])
