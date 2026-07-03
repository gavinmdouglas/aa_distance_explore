rm(list = ls(all.names = TRUE))

library(ggplot2)

# Highlight major differences across the DMS-EX buried vs exposed exchangeability metrics, which are asymmetric.
# This is meant to show the overall shift in exchangeability in buried vs. exposed sites, as well as the 
# differences in exchangeability for amino acid pairs depending on the direction.

combined_buried_exposed <- read.table('~/Drive/research/aa_distance/aa_distance_zenodo/aa_metrics/proteinGym_asymmetric_concat_rsa_group_similarity.tsv.gz',
                                      header=TRUE, sep="\t", stringsAsFactors = FALSE)

buried_tab <- combined_buried_exposed[combined_buried_exposed$rsa_group == "buried", ]
exposed_tab <- combined_buried_exposed[combined_buried_exposed$rsa_group == "exposed", ]

rownames(buried_tab) <- paste(buried_tab$ref_aa, buried_tab$mut_aa, sep = "_")
rownames(exposed_tab) <- paste(exposed_tab$ref_aa, exposed_tab$mut_aa, sep = "_")

combined_tab <- buried_tab[, "combined_buried_exposed_similiarity", drop = FALSE]
colnames(combined_tab) <- "dms.ex_buried"
combined_tab$dms.ex_exposed <- exposed_tab[rownames(combined_tab), "combined_buried_exposed_similiarity"]

# Split rownames back into the pair components
pairs_split <- do.call(rbind, strsplit(rownames(combined_tab), "_"))
combined_tab$ref_aa <- pairs_split[, 1]
combined_tab$mut_aa <- pairs_split[, 2]
combined_tab$pair  <- rownames(combined_tab)
combined_tab$rev_pair <- paste(combined_tab$mut_aa, combined_tab$ref_aa, sep = "_")

# Flag outliers by perpendicular distance from the 1:1 (fit) line
combined_tab$resid <- combined_tab$dms.ex_buried - combined_tab$dms.ex_exposed
cutoff <- 4.5 * sd(combined_tab$resid)
combined_tab$outlier <- abs(combined_tab$resid) > cutoff

# Flag the top 2 major outliers in addition to the cutoff
top_pos <- tail(order(combined_tab$resid), 2)
combined_tab$outlier[top_pos] <- TRUE

fwd <- combined_tab[combined_tab$outlier, ]

# Also keep track of the reverse substitutions, for reference.
rev <- combined_tab[fwd$rev_pair, ]

# build label "ref -> mut" for forward outliers (and reverse for matching reverse)
fwd$label <- paste(fwd$ref_aa, fwd$mut_aa, sep = " \u2192 ")
rev$label <- paste(rev$ref_aa, rev$mut_aa, sep = " \u2192 ")

buried_vs_exposed <- ggplot(combined_tab, aes(x = dms.ex_exposed, y = dms.ex_buried)) +
                            geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
                            geom_point(colour = "grey70", size = 1.5, alpha=0.5) +
                            geom_point(data = fwd, colour = "firebrick", size = 2) +
                            geom_point(data = rev, colour = "steelblue", size = 2) +
                            ggrepel::geom_text_repel(data = fwd, aes(label = label),
                                                     colour = "firebrick", size = 3, segment.colour = "grey60", force=7) +
                            ggrepel::geom_text_repel(data = rev, aes(label = label),
                                                     colour = "steelblue", size = 3, segment.colour = "grey60", force=1) +
                            coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
                            labs(x = "DMS-EX exposed site exchangeability", y = "DMS-EX buried site exchangeability") +
                            theme_bw() 

ggsave(plot = buried_vs_exposed,
       filename = '~/Drive/research/aa_distance/aa_distance_ms/display/Main_asymmetric_DMS-EX_by_burial_status.pdf',
       width = 5, height = 5, device=cairo_pdf)
