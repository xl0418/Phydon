
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Phydon

<!-- badges: start -->
<!-- badges: end -->

Estimating maximum growth rates of bacteria from genomic data with
**Phydon**.

# Introduction

Phydon is a R package that estimates maximum growth rates of bacteria
from genomic data. This package implements a method detailed in [Our
paper](), which leverages phylogenetic signals and mechanistic
statistics, specifically [codon usage bias
(CUB)](https://github.com/jlw-ecoevo/gRodon2) (see
[gRodon](https://github.com/jlw-ecoevo/gRodon2) for details of this
mechanistic method), to enhance the accuracy of growth rate estimations.
By integrating these advanced techniques, Phydon provides a robust tool
for microbial growth rate analysis based on genomic information.

# Installation

You can install Phydon from GitHub with:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("xl0418/Phydon")
```

> \[!IMPORTANT\]  
> Phydon has several essential dependencies, including `gRodon2`, `ape`,
> `Biostrings`, `picante`, etc. Please make sure you have installed
> these packages before using Phydon.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

BiocManager::install("Biostrings")
BiocManager::install("coRdon")
install.packages("matrixStats") 

devtools::install_github("jlw-ecoevo/gRodon2")
```

# Usage

## Minimal example

Phydon is super user friendly. With an internal data, a minimal example
is as simple as:

``` r
library(Phydon)
# load the example data
data_info <- read.csv(system.file("extdata", "data_info.csv", package = "Phydon"))

#### To get the directory of the genomic data  ####
#### This is only for internal data as directory of genomic data is ad hoc per user's system ####
#### For user's input data, see below for details ####
genomes <- data_info$accession_no
gene_loc <- c()
for (i in 1:length(genomes)) {
  gene_loc[i] <- system.file("extdata", paste0("genefiles/",genomes[i],"/",genomes[i],".ffn"), package = "Phydon")
}
data_info$gene_loc <- gene_loc


# estimate the maximum growth rate
result <- Phydon(data_info)
```

## Input data structure

The input data requires a data frame that contains the following
columns:

The columns are:

- `gene_loc`: the directory of the annotated genomic data. The fasta
  file should be named as the accession number of the genome. For
  example, the fasta file of the genome with the accession number
  `RS_GCF_000024245.1` should be named as `RS_GCF_000024245.1.ffn`.

- `accession_no`: the accession number of the genome that users want to
  predict the growth rate.

- `Opt_temp`: **Optional**. This is used for gRodon predictions the
  growth temperature of the genome. If not provided, the growth
  temperature of the genome will be set to 20 degree Celsius.

<details>
<summary>
What the input file looks like
</summary>

| gene_loc                                            | accession_no       |
|-----------------------------------------------------|--------------------|
| genefiles/RS_GCF_002749895.1/RS_GCF_002749895.1.ffn | RS_GCF_002749895.1 |
| genefiles/RS_GCF_002849855.1/RS_GCF_002849855.1.ffn | RS_GCF_002849855.1 |
| genefiles/RS_GCF_002906655.1/RS_GCF_002906655.1.ffn | RS_GCF_002906655.1 |
| genefiles/RS_GCF_003026105.1/RS_GCF_003026105.1.ffn | RS_GCF_003026105.1 |
| genefiles/RS_GCF_003026475.1/RS_GCF_003026475.1.ffn | RS_GCF_003026475.1 |
| genefiles/RS_GCF_003026815.1/RS_GCF_003026815.1.ffn | RS_GCF_003026815.1 |
| genefiles/RS_GCF_003144035.1/RS_GCF_003144035.1.ffn | RS_GCF_003144035.1 |
| genefiles/RS_GCF_003544875.1/RS_GCF_003544875.1.ffn | RS_GCF_003544875.1 |
| genefiles/RS_GCF_003716875.1/RS_GCF_003716875.1.ffn | RS_GCF_003716875.1 |
| genefiles/RS_GCF_900130105.1/RS_GCF_900130105.1.ffn | RS_GCF_900130105.1 |

</details>

## Output

The output of Phydon is a data frame that contains the following
columns:

- `genome`: the accession number of the genome that users want to
  predict the growth rate

- `rep_genome`: the representative genome of the species that the genome
  belongs to based on GTDB database; If it is NA, no matched
  representative genomes are found. The genome will be treated as a
  single species outside of the GTDB database.In this case, only gRodon
  predictions will be used. (See example result)

- `species`: the species names of the genome based on GTDB database; If
  it is NA, no matched species are found.

- `neighbor_repgenome_train`: the accession number of the representative
  genome that is most phylogeneitcally close to the genome in the data
  of the package. If it is NA, no matched representative genomes are
  found. The phylogenetic distance of the user’s genome to the database
  is then set to Inf.

- `phy_distance`: the phylogenetic distance of the user’s genome to the
  database. If it is Inf, no matched representative genomes are found.

- `gRodonpred`: the gRodon predictions of the genome based on the
  [gRodon2](https://github.com/jlw-ecoevo/gRodon2).

- `phylopred`: the phylogenetic predictions of the genome based on the
  [picante](https://www.rdocumentation.org/packages/picante/versions/1.8.2).

- `combopred`: the combined predictions of the genome based on the
  gRodon and phylogenetic predictions. The regression model is trained
  based on the data of the package. See details in the [paper]().

<details>
<summary>
What ‘result’ looks like
</summary>

| genome             | rep_genome         | species                         | neighbor_repgenome_train | phy_distance | gRodonpred | phylopred | combopred |
|--------------------|--------------------|---------------------------------|--------------------------|--------------|------------|-----------|-----------|
| RS_GCF_002749895.1 | RS_GCF_002749895.1 | s\_\_Vibrio fujianensis         | RS_GCF_000621645.1       | 0.12220675   | 1.1346118  | 0.4697509 | 0.4697509 |
| RS_GCF_002849855.1 | RS_GCF_002849855.1 | s\_\_Vibrio azureus             | RS_GCF_900460535.1       | 0.05166449   | 0.4913568  | 0.2783264 | 0.2783264 |
| RS_GCF_002906655.1 | RS_GCF_002906655.1 | s\_\_Vibrio hyugaensis          | RS_GCF_900460535.1       | 0.03447402   | 0.7658604  | 0.2783264 | 0.2783264 |
| RS_GCF_003026105.1 | RS_GCF_003026105.1 | s\_\_Photobacterium sp003026105 | RS_GCF_003025615.1       | 0.16585578   | 1.5986925  | 2.0439385 | 2.0439385 |
| RS_GCF_003026475.1 | RS_GCF_003026475.1 | s\_\_Photobacterium lipolyticum | RS_GCF_003025615.1       | 0.08080289   | 1.4924976  | 1.8805498 | 1.8805498 |
| RS_GCF_003026815.1 | RS_GCF_003026815.1 | s\_\_Photobacterium damselae    | RS_GCF_003025615.1       | 0.15489320   | 0.7407403  | 2.0439385 | 2.0439385 |
| RS_GCF_003144035.1 | RS_GCF_003144035.1 | s\_\_Vibrio albus               | RS_GCF_900460535.1       | 0.15592277   | 0.6097861  | 0.8919112 | 0.8919112 |
| RS_GCF_003544875.1 | RS_GCF_003544875.1 | s\_\_Vibrio alfacsensis         | RS_GCF_900460535.1       | 0.03506859   | 0.6234453  | 0.2783264 | 0.2783264 |
| RS_GCF_003716875.1 | RS_GCF_003716875.1 | s\_\_Vibrio zhugei              | RS_GCF_000621645.1       | 0.18306980   | 0.6364132  | 0.4697509 | 0.4697509 |
| RS_GCF_900130105.1 | RS_GCF_024346755.1 | s\_\_Vibrio aerogenes           | RS_GCF_900460535.1       | 0.17884206   | 1.6145463  | 0.5057288 | 0.5057288 |
| RS_GCF_102749895.1 | NA                 | NA                              | NA                       | Inf          | 1.1346118  | NA        | NA        |

</details>

> \[!Note\]  
> The last genome `RS_GCF_102749895.1` was made up on purpose to show
> the case when the genome is not in the GTDB database. In this case,
> the representative genome and species will be NA, and the phylogenetic
> distance will be Inf. The gRodon predictions will be used for this
> genome.

# Before you start, annotating genomes

> \[!IMPORTANT\]  
> As gRodon requires annotated genomes, users need to annotate the
> genomes before using Phydon. We recommend using
> [Prokka](https://github.com/tseemann/prokka) for genome annotation.
> Thus, the directory structure of the genomic data should be like this:

- `genefiles/`
  - `genome1/`
    - `genome1.ffn`
    - `genome1.gff`
    - `genome1_CDS_names.txt`
  - `genome2/`
    - `genome2.ffn`
    - `genome2.gff`
    - `genome2_CDS_names.txt`
  - `genome3/`
    - `genome3.ffn`
    - `genome3.gff`
    - `genome3_CDS_names.txt`
  - …

> \[!NOTE\]  
> The `genome1.ffn` is the fasta file of the genome, and the
> `genome1.gff` is the annotation file of the genome. These two files
> can be obtained when annotating genomes. The `genome1_CDS_names.txt`
> is the file that contains the names of the CDSs in the genome. Phydon
> will try to generate this file using `sed` command. If it fails, users
> need to generate this file manually. `sed` is usually available in
> Linux and MacOS. For Windows users, you can install `sed` from
> [here](https://sourceforge.net/projects/gnuwin32/files/sed/4.2.1/).
> See more details in the
> [gRodon2](https://github.com/jlw-ecoevo/gRodon2) package. We recommend
> run Phydon in Linux or MacOS to avoid the potential issues in
> generating the `genome1_CDS_names.txt` file.

# The basic workflow

1.  Annotate the genomes using Prokka and obtain `.ffn` and `.gff`
    files. Optionally, generate the `genome1_CDS_names.txt` file. Not
    necessary if you run Phydon in Linux or MacOS.

2.  Organize the genomic data in the `genefiles/` directory and generate
    `data_info.csv` file.

3.  Run Phydon with the `data_info.csv` file.

# Citation