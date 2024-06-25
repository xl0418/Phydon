#' @title PhyloPred
#' @description This function estimates the maximum growth rate of a species based on the phylogenetic history from the phylogeny.
#' @param rep_genomes_to_est The accession number of the genomes that represent a cluster of genomes shown in the phylogeny tree. This is used to calculate the phylogenetic distance between genomes.
#' @param GTDB_tax_trait_repGenome_in_tree_expanded The data set in the package that contains the trait, accession number, species names, etc. information.
#' @param phydis_tree The phylogenetic tree pulled from GTDB database that contains user's genome and the genomes in the data set of the package.
#'
#' @noRd

PhloPred <- function(rep_genomes_to_est,
                     GTDB_tax_trait_repGenome_in_tree_expanded,
                     phydis_tree) {
  train_rep_genomes <- setdiff(
    unique(
      GTDB_tax_trait_repGenome_in_tree_expanded$Representative.genome
    ),
    unique(rep_genomes_to_est)
  )
  train_rep_genomes <- intersect(
    train_rep_genomes,
    phydis_tree$tip.label
  )


  known_traits <-
    GTDB_tax_trait_repGenome_in_tree_expanded[which(
      GTDB_tax_trait_repGenome_in_tree_expanded$Representative.genome %in% train_rep_genomes
    ), c("Representative.genome", "doubling_h")]
  # unique train_rep_genomes by Representative.genome and doubling_h
  known_traits <- unique(known_traits, by = c("Representative.genome", "doubling_h"))

  known_traits_vec <- known_traits$doubling_h
  names(known_traits_vec) <- known_traits$Representative.genome
  phylopred <-
    picante::phyEstimate(phydis_tree, known_traits_vec, method = "pic")
  est_rep_genomes <- row.names(phylopred)


  temp_sub_phylopred_df <- data.frame()
  for (rep_genome in est_rep_genomes) {
    temp_sub_phylopred_df <- rbind(
      temp_sub_phylopred_df,
      data.frame(rep_genome = rep_genome, phylopred = phylopred[rep_genome, "estimate"])
    )
  }



  return(temp_sub_phylopred_df)
}
