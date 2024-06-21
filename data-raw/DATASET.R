## code to prepare `DATASET` dataset goes here

## load the phylogenetic tree
gtdb_tree <- ape::read.tree(paste0("data-raw/bac120_r220.tree"))

# usethis::use_data(gtdb_tree, overwrite = TRUE)


## load the clustered genome information
sp_clusters <- read.table(
  paste0("data-raw/sp_clusters_r220.tsv"),
  header = T,
  sep = "\t"
)
sp_clusters <- sp_clusters[,c("Representative.genome", "GTDB.species", "Clustered.genomes")]
# usethis::use_data(sp_clusters, overwrite = TRUE)


## load the trait information
GTDB_tax_trait_repGenome_in_tree_expanded <- read.table(
  "data-raw/GTDB_tax_trait_repGenome_in_tree_expanded.csv",
  header = T,
  sep = ","
)
# usethis::use_data(GTDB_tax_trait_repGenome_in_tree_expanded, overwrite = TRUE)


## load the regression model
load("data-raw/reg_model.RData")
load("data-raw/reg_model_tmp.RData")
# usethis::use_data(reg_model, overwrite = TRUE)


usethis::use_data(gtdb_tree, sp_clusters, GTDB_tax_trait_repGenome_in_tree_expanded, reg_model, reg_model_tmp, internal = TRUE, overwrite = TRUE)
