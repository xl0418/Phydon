## code to prepare `DATASET` dataset goes here

## load the phylogenetic tree
gtdb_tree <- ape::read.tree(paste0("data-raw/bac120_r220.tree"))
gtdb_tree_archaea <- ape::read.tree(paste0("data-raw/ar53_r220.tree"))

# usethis::use_data(gtdb_tree, overwrite = TRUE)


## load the clustered genome information
sp_clusters_bac <- read.table(paste0("data-raw/sp_clusters_r220_bac.tsv"),
                          header = T,
                          sep = "\t")
sp_clusters_arc <- read.table(paste0("data-raw/sp_clusters_r220_arc.tsv"),
                              header = T,
                              sep = "\t")
# sp_clusters <- sp_clusters[, c("Representative.genome", "GTDB.species", "Clustered.genomes")]
# usethis::use_data(sp_clusters, overwrite = TRUE)


## load the trait information
GTDB_tax_trait_repGenome_in_tree_expanded <- read.table(
  "data-raw/GTDB_tax_trait_repGenome_in_tree_expanded.csv",
  header = T,
  sep = ","
)
GTDB_tax_trait_repGenome_in_tree_expanded_archaea <- read.table(
  "data-raw/GTDB_tax_trait_repGenome_in_tree_expanded_archaea.csv",
  header = T,
  sep = ","
)
# usethis::use_data(GTDB_tax_trait_repGenome_in_tree_expanded, overwrite = TRUE)


## load the regression model
### bacteria full temperature and not temperature
load("data-raw/reg_model_tmp_bacteria_full.RData")
load("data-raw/reg_model_bacteria_full.RData")
### bacteria metagenome temperature and not temperature
load("data-raw/reg_model_tmp_bacteria_metagenome.RData")
load("data-raw/reg_model_bacteria_metagenome.RData")

### archaea full temperature and not temperature
load("data-raw/reg_model_tmp_archaea_full.RData")
load("data-raw/reg_model_archaea_full.RData")
### archaea metagenome temperature and not temperature
load("data-raw/reg_model_tmp_archaea_metagenome.RData")
load("data-raw/reg_model_archaea_metagenome.RData")


usethis::use_data(
  gtdb_tree,
  sp_clusters_bac,
  sp_clusters_arc,
  GTDB_tax_trait_repGenome_in_tree_expanded,
  reg_model_bacteria_full,
  reg_model_tmp_bacteria_full,
  reg_model_bacteria_metagenome,
  reg_model_tmp_bacteria_metagenome,
  gtdb_tree_archaea,
  GTDB_tax_trait_repGenome_in_tree_expanded_archaea,
  reg_model_archaea_full,
  reg_model_tmp_archaea_full,
  reg_model_archaea_metagenome,
  reg_model_tmp_archaea_metagenome,
  internal = TRUE,
  overwrite = TRUE
)
