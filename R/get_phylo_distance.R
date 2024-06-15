library(ape)
library(readr)
library(castor)
## read tree 


get_phylo_distance <- function(spe){
  gtdb_tree <- read.tree("data/bac120_r214.tree")
  gtdb_species <- read.table("data/sp_clusters_r214.tsv", header = T, sep = "\t")
  rep_genomes_database <- read.table("data/GTDB_tax_trait_repGenome_in_tree_expanded.csv", header = T, sep = ",")
  rep_genomes_database <- rep_genomes_database$Representative.genome
  
  # representive genome
  rep_genome <- gtdb_species[gtdb_species$GTDB.species == spe, "Representative.genome"]
  if (length(rep_genome) == 0) {
    return("The species is not in the database.")
  }
  
  # if rep_genome in rep_genomes_database, return 
  if (rep_genome %in% rep_genomes_database) {
    return("The representative genome is in the database, just run the training model.")
  }
  
  all_rep_genomes <- c(rep_genomes_database, rep_genome)
  
  
  rep_genomes_intree <- gtdb_tree$tip.label
  # get the intersection of matched_genomes and rep_genomes_intree
  matched_genomes <- intersect(all_rep_genomes, rep_genomes_intree)

  
  sub.gtdb.tree <-
    get_subtree_with_tips(gtdb_tree, only_tips = matched_genomes)
  gtdb.subtree <- sub.gtdb.tree$subtree
  
  min_dist_test_train <-
    min(cophenetic(gtdb.subtree)[rep_genome, rep_genomes_database])
  
  # the index of the minimum distance
  index_min_dist <-
    which(cophenetic(gtdb.subtree)[rep_genome, rep_genomes_database] == min_dist_test_train)
  
  neighbor_repgenome_train <-
    unique(rep_genomes_database[index_min_dist])

  
  
  info_df <- data.frame(
    species = spe,
    rep_genome = rep_genome,
    neighbor_repgenome_train = neighbor_repgenome_train,
    min_dist_test_train = min_dist_test_train
  )
  # check if the folder output exists, if not, create it
  if (!dir.exists("output")) {
    dir.create("output")
  }
  write.tree(gtdb.subtree, file = "output/subtree.tree")
  
  
  return(info_df)
}


# Example
spe <- "s__Methylobacterium radiotolerans"
get_phylo_distance(spe)

