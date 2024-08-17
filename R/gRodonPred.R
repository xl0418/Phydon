#' @title gRodonpred
#' @description This function estimates the maximum growth rate of a species based on the gene sequences using gRodon package.
#' @param gene_file The file that contains the gene sequences in fasta format.
#' @param temp The growth temperature of the species.
#'
#' @noRd

gRodonpred <- function(gene_file, temp) {
  ### format in xxx.ffn
  genes <- Biostrings::readDNAStringSet(gene_file)

  ### replace ".ffn" with "_CDS_names.txt"
  CDS_file <- gsub(".ffn", "_CDS_names.txt", gene_file)
  gff_file <- gsub(".ffn", ".gff", gene_file)

  ### check if the file exists or the size of file is zero
  if (!file.exists(CDS_file) || file.size(CDS_file) == 0) {
    print("The CDS file does not exist. Trying to create it using sed ... ")
    # This part is not working for now as it requests the system to have sed and awk installed.
    sed_lind <- utils::capture.output(
      cat(
        "sed -n '/##FASTA/q;p'",
        gff_file,
        "| awk '$3==\"CDS\"' | awk '{print $9'} | awk 'gsub(\";.*\",\"\")' | awk 'gsub(\"ID=\",\"\")' > ",
        CDS_file
      )
    )
    system(sed_lind, intern = TRUE)
    if (!file.exists(CDS_file)) {
      print(
        "Creating failed. It might be due to that sed cannot be used, particularly on Windows. Please check https://microbialgamut.com/gRodon-vignette to mannually generate CDS files."
      )
      stop("The CDS file does not exist.")
    }
  }

  # Subset your sequences to those that code for proteins
  CDS_IDs <- readLines(paste0(CDS_file))
  gene_IDs <- gsub(" .*", "", names(genes)) #Just look at first part of name before the space
  genes <- genes[gene_IDs %in% CDS_IDs]

  #Search for genes annotated as ribosomal proteins
  highly_expressed <- grepl("ribosomal protein", names(genes), ignore.case = T)

  # Deal with errors from gRodon
  tryCatch({
    maxg <- gRodon::predictGrowth(genes, highly_expressed, temperature = temp)
    # list to data frame
    maxg <- as.data.frame(maxg)
    maxg$OGT <- temp
  }, error = function(e) {
    print(paste("Error in predicting the growth rate. Please check the gene sequences at", gene_file))
    maxg <<- data.frame(CUBHE = NA,
                        GC = NA,
                        GCdiv = NA,
                        ConsistencyHE = NA,
                        CUB = NA,
                        CPB = NA,
                        FilteredSequences = NA,
                        nHE = NA,
                        dCUB = NA,
                        OGT = temp,
                        d = NA,
                        LowerCI = NA,
                        UpperCI = NA
                        )


  })


  return(maxg)
}
