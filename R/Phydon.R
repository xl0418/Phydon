#' Estimate maximum growth rates/doubling time of species based on genomic data and phylogenetic history.
#'
#' @param data_info_df The data frame that contains the information of the genomes to be estimated. The first column should be the directory of the genomic data, i.e., xxx.ffn. The second column should be the accession number of genomes. Or genome names if user provides phylogenetic tree. The third column should be the optimal growth temperature.
#' @param user_tree The user-provided phylogenetic tree. It can be obtained by running gtdbtk on users' genomes and generated under the foler "/output/classify/". If not provided (user_tree = NULL), the user is supposed to provide the accession number of the genomes. Default is NULL.
#' @param tax The taxon of the species. It can be "bacteria" or "archaea". Default is "bacteria".
#' @param regression_mode The regression model to be used. It can be "arithmetic_mean" or "geometric_mean". Default is "geometric_mean".
#' @return A data frame that contains the species name and the estimated maximum growth rates by three methods.
#' @export


Phydon <- function(data_info_df, user_tree = NULL, tax = "bacteria", regression_mode="geometric_mean") {
  ## load the necessary data
  if (tax == "bacteria") {
    print("Estimating maximum growth rates for bacteria ...")
    GTDB_tax_trait_repGenome_in_tree_expanded <- get0("GTDB_tax_trait_repGenome_in_tree_expanded",
                                                      envir = asNamespace("Phydon"))
    gtdb_tree <- get0("gtdb_tree", envir = asNamespace("Phydon"))
    sp_clusters <- get0("sp_clusters_bac", envir = asNamespace("Phydon"))
    reg_model <- get0("reg_model", envir = asNamespace("Phydon"))
    reg_model_tmp <- get0("reg_model_tmp", envir = asNamespace("Phydon"))
  } else if (tax == "archaea") {
    print("Estimating maximum growth rates for archaea ...")
    GTDB_tax_trait_repGenome_in_tree_expanded <- get0("GTDB_tax_trait_repGenome_in_tree_expanded_archaea",
                                                      envir = asNamespace("Phydon"))
    gtdb_tree <- get0("gtdb_tree_archaea", envir = asNamespace("Phydon"))
    sp_clusters <- get0("sp_clusters_arc", envir = asNamespace("Phydon"))
    reg_model <- get0("reg_model_archaea", envir = asNamespace("Phydon"))
    reg_model_tmp <- get0("reg_model_tmp_archaea", envir = asNamespace("Phydon"))

  } else {
    return("The taxon is not supported.")
  }
  print("Loading internal data ...")


  if (!is.null(user_tree)) {
    print("Phylogenetic tree is provided. Phylogenetic distance of genomes will be calculated based on the tree ...")
    ## build the tree
    inputtree <- user_tree
    usertree <- TRUE
  } else{
    print("Phylogenetic tree is not provided. Will search genomes by accession number in the database ...")
    inputtree <- gtdb_tree
    usertree <- FALSE
  }

  ## check if gene_location and genome_name are in the column names
  if (!("gene_location" %in% colnames(data_info_df))) {
    return("The column name of the gene location is not correct. Please provide the gene locations in the column named 'gene_location'.")
  }

  if (!("genome_name" %in% colnames(data_info_df))) {
    return("The column name of the genome names is not correct. Please provide the accession numbers in the column named 'genome_name'.")
  }

  ## check the length of species and the gene files
  genomes_to_est <- data_info_df$genome_name

  if (length(genomes_to_est) == 0) {
    return("The gene file is empty.")
  } else {
    print(
      paste0(
        length(genomes_to_est),
        " genomes are found in the input file. Start estimating maximum growth rates ..."
      )
    )
  }


  ## check if the temperature column exists
  if ("temperature" %in% colnames(data_info_df)) {
    temps <- data_info_df$temperature
    print("The growth temperatures are provided.")
  } else {
    print("The growth temperatures are not provided.")
  }




  ## phylopred estimates
  print("Start phylopred predictions ...")
  prepare_for_genomes_to_est <- get_phylo_distance(
    genomes_to_est,
    inputtree,
    GTDB_tax_trait_repGenome_in_tree_expanded,
    sp_clusters,
    usertree
  )

  rep_genomes_df <- prepare_for_genomes_to_est[[1]]
  if ("temperature" %in% colnames(data_info_df)) {
    temp_info <- data_info_df[, c("genome_name", "temperature")]
    colnames(temp_info) <- c("genome", "temperature")
    rep_genomes_df <- merge(rep_genomes_df, temp_info, by = "genome")
  }

  # for known species, both gRodon and phylopred estimates are provided
  rep_genomes_to_est <- unique(rep_genomes_df$rep_genome)
  rep_genomes_to_est <- rep_genomes_to_est[!is.na(rep_genomes_to_est)]
  if(length(rep_genomes_to_est) == 0){
    phylopred_df <- rep_genomes_df
    phylopred_df$phylopred <- NA
    print("Cannot find any representative genomes of user's genomes in the data base. ")
    print("Please check whether the kingdom is correct. Or please provide user's tree.")
    print("Only the gRodon predictions will be provided.")
  } else {
    if (length(rep_genomes_to_est) < length(genomes_to_est)) {
      print(paste0("Only ", length(rep_genomes_to_est), " genomes are matched to the representative genomes in the database."))
    }

    phydis_tree <- prepare_for_genomes_to_est[[2]]
    Est_phylopred <- PhloPred(rep_genomes_to_est,
                              GTDB_tax_trait_repGenome_in_tree_expanded,
                              phydis_tree)

    phylopred_df <- merge(rep_genomes_df, Est_phylopred, by = "rep_genome", all = TRUE)
  }



  ## Future feature to be implemented: Can be parallelized
  print("Start gRodon prediction ...")
  ## gRodon estimates
  est_gRodon_df <- data.frame()

  for (i in 1:nrow(data_info_df)) {
    gene_location <- data_info_df$gene_location[i]
    genome <- data_info_df$genome_name[i]
    gff_file <- gsub(".ffn", ".gff", gene_location)

    # check if the required ffn and gff files exist
    if (!file.exists(gene_location)) {
      print(paste0("The ffn file ", gene_location, " does not exist."))
      next
    }


    if (exists("temps")) {
      temp <- temps[i]
    } else {
      temp <- NA
    }
    Est_gRodon <- gRodonpred(gene_location, temp = temp)
    Est_gRodon$genome <- genome
    colnames(Est_gRodon)[which(colnames(Est_gRodon) == "d")] <- "gRodonpred"
    est_gRodon_df <- rbind(est_gRodon_df,
                           Est_gRodon)
  }




  ## combopred estimates
  print("Start Phydon regression predictions ...")
  input_df <- merge(est_gRodon_df, phylopred_df, by = "genome")

  combopred_df <- combopred(input_df, regression_mode)


  ## rearrange the columns

  combopred_df <- combopred_df[, c(
    "genome",
    "rep_genome",
    "species",
    "neighbor_repgenome_train",
    "phy_distance",
    "OGT",
    "gRodonpred",
    "phylopred",
    "combopred",
    "CUBHE",
    "GC",
    "GCdiv",
    "ConsistencyHE",
    "CUB",
    "CPB",
    "FilteredSequences",
    "nHE",
    "dCUB",
    "LowerCI",
    "UpperCI"
  )]


  print("Estimation is done.")

  return(combopred_df)


}
