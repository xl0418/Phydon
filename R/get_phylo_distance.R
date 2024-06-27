#' Compute the phylogenetic distance of the user's species to the database and generate a subtree.
#'
#' @param genomes_to_est The species name of the user. It should be the in the same format as in the GTDB database.
#' @param inputtree The phylogenetic tree of the GTDB database or the user's tree getting from gtdbtk.
#' @param GTDB_tax_trait_repGenome_in_tree_expanded The dataframe that contains the representative genomes in the tree.
#' @param sp_clusters The dataframe that contains the clustered genomes.
#' @param usertree The user's phylogenetic tree is provided. Default is FALSE.
#' @param save_tree The name of the file to save the subtree. Default is NULL and not save the tree.
#' @return A list of the dataframe that contains the information of the user's species and its phylogenetic distance to the database and the subtree.
#' @noRd


get_phylo_distance <- function(genomes_to_est,
                               inputtree,
                               GTDB_tax_trait_repGenome_in_tree_expanded,
                               sp_clusters,
                               usertree=FALSE,
                               save_tree = NA) {
  if (!usertree){
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
          temp_df$species = sp_clusters[i, "GTDB.species"]
          temp_df$rep_genome = sp_clusters[i, "Representative.genome"]

          break
        }

      }
      species_genome_repgenome_df <- rbind(
        species_genome_repgenome_df,
        temp_df
      )


    }

    unknown_genomes <- species_genome_repgenome_df[is.na(species_genome_repgenome_df$species),]

    species_genome_repgenome_df <- species_genome_repgenome_df[!is.na(species_genome_repgenome_df$species),]
    info_df <- data.frame()
    if(nrow(species_genome_repgenome_df) > 0){
      rep_genome <- species_genome_repgenome_df$rep_genome

      all_rep_genomes <- unique(c(rep_genomes_database, rep_genome))


      rep_genomes_intree <- inputtree$tip.label
      # get the intersection of matched_genomes and rep_genomes_intree
      matched_genomes <- intersect(all_rep_genomes, rep_genomes_intree)


      sub.gtdb.tree <-
        castor::get_subtree_with_tips(inputtree, only_tips = matched_genomes)
      gtdb.subtree <- sub.gtdb.tree$subtree


      for (i in 1:nrow(species_genome_repgenome_df)) {
        single_rep_genome <- species_genome_repgenome_df$rep_genome[i]
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
            species = species_genome_repgenome_df$species[i],
            rep_genome = single_rep_genome,
            genome = species_genome_repgenome_df$genome[i],
            neighbor_repgenome_train = neighbor_repgenome_train,
            phy_distance = min_dist_test_train
          )
        )

      }
    }


    # info_df for unknown genomes
    unknown_info_df <- data.frame(
      species = NA,
      rep_genome = NA,
      genome = unknown_genomes$genome,
      neighbor_repgenome_train = NA,
      phy_distance = NA
    )

    info_df <- rbind(info_df, unknown_info_df)

    if (!is.na(save_tree)) {
      # check if the folder output exists, if not, create it
      if (!dir.exists("output")) {
        dir.create("output")
      }
      ape::write.tree(gtdb.subtree, file = paste0("output/", save_tree, ".tree"))
    }
  } else {
    training_tree_tips <- GTDB_tax_trait_repGenome_in_tree_expanded$Representative.genome

    training_testing_tree_tips <- c(training_tree_tips, genomes_to_est)

    sub.tree.list <-
      castor::get_subtree_with_tips(inputtree, only_tips = training_testing_tree_tips)


    gtdb.subtree <- sub.tree.list$subtree

    sub_training_tree_tips <- intersect(training_tree_tips, gtdb.subtree$tip.label)

    info_df <- data.frame()
    for (test_genome in genomes_to_est) {
      min_dist_test_train <-
        min(ape::cophenetic.phylo(gtdb.subtree)[test_genome, sub_training_tree_tips])

      # the index of the minimum distance
      index_min_dist <-
        which(ape::cophenetic.phylo(gtdb.subtree)[test_genome, sub_training_tree_tips] == min_dist_test_train)

      neighbor_repgenome_train <-
        unique(sub_training_tree_tips[index_min_dist])

      info_df <- rbind(
        info_df,
        data.frame(
          species = NA,
          rep_genome = test_genome,
          genome = test_genome,
          neighbor_repgenome_train = neighbor_repgenome_train,
          phy_distance = min_dist_test_train
        )
      )

    }

  }




  return(list(info_df, gtdb.subtree))
}
