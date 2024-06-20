#' Estimate maximum growth rates/doubling time of species based on genomic data and phylogenetic history.
#'
#' @param data_info_df The data frame that contains the information of the genomes to be estimated. The first column should be the directory of the genomic data, i.e., xxx.ffn. The second column should be the accession number of genomes. The third column should be the optimal growth temperature. Default is 20.
#' @param opt_temp The optimal growth temperature of the species. Default is 20.
#' @return A data frame that contains the species name and the estimated maximum growth rates by three methods.
#' @export


Phydon <- function(data_info_df,
                       opt_temp = 20) {

  ## load the necessary data
  print("Loading internal data ...")
  GTDB_tax_trait_repGenome_in_tree_expanded <- get0("GTDB_tax_trait_repGenome_in_tree_expanded", envir = asNamespace("Phydon"))
  gtdb_tree <- get0("gtdb_tree", envir = asNamespace("Phydon"))
  sp_clusters <- get0("sp_clusters", envir = asNamespace("Phydon"))
  reg_model <- get0("reg_model", envir = asNamespace("Phydon"))
  reg_model_tmp <- get0("reg_model_tmp", envir = asNamespace("Phydon"))

  ## check the length of species and the gene files
  genomes_to_est <- data_info_df$accession_no

  print(
    paste0(
      length(genomes_to_est),
      " genomes are found in the gene file. Start estimating maximum growth rates ..."
    )
  )

  ## check if the third column exists
  if (ncol(data_info_df) == 3) {
    opt_temps <- data_info_df$opt_temp
    print("The optimal growth temperatures are provided.")
  } else {
    print("The optimal growth temperatures are not provided. Default is 20.")
  }

  if (nrow(data_info_df) == 0) {
    return("The gene file is empty.")
  }



  representative.genomes <- GTDB_tax_trait_repGenome_in_tree_expanded$Representative.genome


  est_gRodon_df <- data.frame()
  ## Future feature to be implemented: Can be parallelized

  ## gRodon estimates
  for (i in 1:nrow(data_info_df)) {
    gene_loc <- data_info_df$gene_loc[i]
    genome <- data_info_df$accession_no[i]
    if (exists("opt_temps")) {
      opt_temp <- opt_temps[i]
    }
    Est_gRodon <- gRodonpred(gene_loc, opt_temp = opt_temp)
    est_gRodon_df <- rbind(est_gRodon_df,
                           data.frame(genome = genome, gRodonpred = Est_gRodon))
  }



  ## phylopred estimates
  prepare_for_genomes_to_est <- get_phylo_distance(
    genomes_to_est,
    gtdb_tree,
    GTDB_tax_trait_repGenome_in_tree_expanded,
    sp_clusters
  )

  rep_genomes_df <- prepare_for_genomes_to_est[[1]]
  # deal with unkown species, only gRodon estimates are provided
  unknown <- rep_genomes_df[which(is.na(rep_genomes_df$species)), ]


  # for known species, both gRodon and phylopred estimates are provided
  rep_genomes_df <- rep_genomes_df[which(!is.na(rep_genomes_df$species)), ]
  rep_genomes_to_est <- unique(rep_genomes_df$rep_genome)
  phydis_tree <- prepare_for_genomes_to_est[[2]]
  Est_phylopred <- PhloPred(rep_genomes_to_est,
                            GTDB_tax_trait_repGenome_in_tree_expanded,
                            phydis_tree)

  phylopred_df <- merge(rep_genomes_df, Est_phylopred, by = "rep_genome")



  ## combopred estimates
  input_df <- merge(est_gRodon_df, phylopred_df, by = "genome")
  combopred_df <- combopred(input_df, reg_model)

  ## rbind unknown species
  if(nrow(unknown) > 0){
    unknown <- merge(est_gRodon_df, unknown, by = "genome")
    unknown$phylopred <- NA
    unknown$combopred <- NA
    combopred_df <- rbind(combopred_df, unknown)
  }

  ## rearrange the columns
  combopred_df <- combopred_df[, c(
    "genome",
    "rep_genome",
    "species",
    "neighbor_repgenome_train",
    "phy_distance",
    "gRodonpred",
    "phylopred",
    "combopred"
  )]

  return(combopred_df)


}
