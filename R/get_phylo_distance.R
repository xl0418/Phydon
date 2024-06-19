#' Compute the phylogenetic distance of the user's species to the database and generate a subtree.
#'
#' @param genomes_to_est The species name of the user. It should be the in the same format as in the GTDB database.
#' @param gtdb_tree The phylogenetic tree of the GTDB database.
#' @param GTDB_tax_trait_repGenome_in_tree_expanded The dataframe that contains the representative genomes in the tree.
#' @param sp_clusters The dataframe that contains the clustered genomes.
#' @param save_tree The name of the file to save the subtree. Default is NULL and not save the tree.
#' @return A list of the dataframe that contains the information of the user's species and its phylogenetic distance to the database and the subtree.
#' @noRd


get_phylo_distance <- function(genomes_to_est,
                               gtdb_tree,
                               GTDB_tax_trait_repGenome_in_tree_expanded,
                               sp_clusters,
                               save_tree = NA) {
  rep_genomes_database <- GTDB_tax_trait_repGenome_in_tree_expanded$Representative.genome

  ### This can be very computational expensive; To be parallelized
  species_genome_repgenome_df <- data.frame()
  for (genome_to_est in genomes_to_est) {
    temp_df <- data.frame(
      species = NA,
      rep_genome = NA,
      genome = genome_to_est
    )
    for (i in 1:nrow(sp_clusters)) {
      # if genome_to_est is in the string of Clustered.genomes
      clustered_genomes <- sp_clusters[i, "Clustered.genomes"]
      if (genome_to_est %in% strsplit(clustered_genomes, ",")[[1]]) {
        temp_df$species <- sp_clusters[i, "GTDB.species"]
        temp_df$rep_genome <- sp_clusters[i, "Representative.genome"]
        break
      }
    }
    species_genome_repgenome_df <- rbind(
      species_genome_repgenome_df,
      temp_df
    )

  }

  unknown_genomes <- species_genome_repgenome_df[which(is.na(species_genome_repgenome_df$species)), ]
  known_genomes <- species_genome_repgenome_df[which(!is.na(species_genome_repgenome_df$species)), ]

  rep_genome <- known_genomes$rep_genome
  # # if rep_genome in rep_genomes_database, return
  # if (all(rep_genome %in% rep_genomes_database)) {
  #   print("The representative genome is in the database, just run the training model.")
  #   return(list())
  # }

  all_rep_genomes <- unique(c(rep_genomes_database, rep_genome))


  rep_genomes_intree <- gtdb_tree$tip.label
  # get the intersection of matched_genomes and rep_genomes_intree
  matched_genomes <- intersect(all_rep_genomes, rep_genomes_intree)


  sub.gtdb.tree <-
    castor::get_subtree_with_tips(gtdb_tree, only_tips = matched_genomes)
  gtdb.subtree <- sub.gtdb.tree$subtree

  info_df <- data.frame()
  for (i in 1:nrow(known_genomes)) {
    single_rep_genome <- known_genomes$rep_genome[i]
    min_dist_test_train <-
      min(ape::cophenetic.phylo(gtdb.subtree)[single_rep_genome, rep_genomes_database])

    # the index of the minimum distance
    index_min_dist <-
      which(ape::cophenetic.phylo(gtdb.subtree)[single_rep_genome, rep_genomes_database] == min_dist_test_train)

    neighbor_repgenome_train <-
      unique(rep_genomes_database[index_min_dist])

    info_df <- rbind(
      info_df,
      data.frame(
        species = known_genomes$species[i],
        rep_genome = single_rep_genome,
        genome = known_genomes$genome[i],
        neighbor_repgenome_train = neighbor_repgenome_train,
        phy_distance = min_dist_test_train
      )
    )

  }

  if(nrow(unknown_genomes)>0){
    unknown_genomes$neighbor_repgenome_train <- NA
    unknown_genomes$phy_distance <- Inf
  }

  info_df <- rbind(info_df, unknown_genomes)

  if (!is.na(save_tree)) {
    # check if the folder output exists, if not, create it
    if (!dir.exists("output")) {
      dir.create("output")
    }
    ape::write.tree(gtdb.subtree, file = paste0("output/", save_tree, ".tree"))
  }



  return(list(info_df, gtdb.subtree))
}
