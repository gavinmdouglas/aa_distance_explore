# Definitions and useful objects for handling radical and conservative AA changes.
# Defined based on (1) charge, as well as (2) polarity and volume.
# As defined by Zhang 2000 and Dagan et al. 2002.

aa_charge_groups <- list(
  `positive` = c('R', 'H', 'K'),
  `negative` = c('D', 'E'),
  `uncharged` = c('A', 'N', 'C', 'Q', 'G', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
)

aa_polarity_volume_groups <- list(
  `special` = c('C'),
  `neutral and small` = c('A', 'G', 'P', 'S', 'T'),
  `polar and relatively small` = c('N', 'D', 'Q', 'E'),
  `polar and relatively large` = c('R', 'H', 'K'),
  `nonpolar and relatively small` = c('I', 'L', 'M', 'V'),
  `nonpolar and relatively large` = c('F', 'W', 'Y')
)

aa_to_groups <- list(
  `R` = c('positive', 'polar and relatively large'),
  `H` = c('positive', 'polar and relatively large'),
  `K` = c('positive', 'polar and relatively large'),
  `D` = c('negative', 'polar and relatively small'),
  `E` = c('negative', 'polar and relatively small'),
  `A` = c('uncharged', 'neutral and small'),
  `N` = c('uncharged', 'polar and relatively small'),
  `C` = c('uncharged', 'special'),
  `Q` = c('uncharged', 'polar and relatively small'),
  `G` = c('uncharged', 'neutral and small'),
  `I` = c('uncharged', 'nonpolar and relatively small'),
  `L` = c('uncharged', 'nonpolar and relatively small'),
  `M` = c('uncharged', 'nonpolar and relatively small'),
  `F` = c('uncharged', 'nonpolar and relatively large'),
  `P` = c('uncharged', 'neutral and small'),
  `S` = c('uncharged', 'neutral and small'),
  `T` = c('uncharged', 'neutral and small'),
  `W` = c('uncharged', 'nonpolar and relatively large'),
  `Y` = c('uncharged', 'nonpolar and relatively large'),
  `V` = c('uncharged', 'nonpolar and relatively small')
)


# Get mappings of radical and conservative AA changes based on the above definitions.
all_aa <- sort(names(aa_to_groups))

charge_sub_change <- list()
polarity_volume_sub_change <- list()
for (aa in all_aa) {
  charge_sub_change[[aa]] <- list()
  polarity_volume_sub_change[[aa]] <- list()
}

for (aa1 in all_aa[1:(length(all_aa) - 1)]) {

 for (aa2 in all_aa) {
  
   if (aa1 == aa2) { next }
  
   if (aa_to_groups[[aa1]][1] == aa_to_groups[[aa2]][1]) {
     charge_sub_change[[aa1]][[aa2]] <- 'Conservative'
     charge_sub_change[[aa2]][[aa1]] <- 'Conservative'
   } else {
     charge_sub_change[[aa1]][[aa2]] <- 'Radical'
     charge_sub_change[[aa2]][[aa1]] <- 'Radical'
   }

   if (aa_to_groups[[aa1]][2] == aa_to_groups[[aa2]][2]) {
     polarity_volume_sub_change[[aa1]][[aa2]] <- 'Conservative'
     polarity_volume_sub_change[[aa2]][[aa1]] <- 'Conservative'
   } else {
     polarity_volume_sub_change[[aa1]][[aa2]] <- 'Radical'
     polarity_volume_sub_change[[aa2]][[aa1]] <- 'Radical'
   }
 }
}
