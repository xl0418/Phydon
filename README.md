
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Phydon

<!-- badges: start -->
<!-- badges: end -->

The goal of Phydon is to …

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

You can install the development version of Phydon from GitHub with:

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

BiocManager::install("Biostrings")
BiocManager::install("coRdon")
install.packages("matrixStats") 

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("jlw-ecoevo/gRodon2")
```

# Usage

First, load the package and the example data.

``` r
library(Phydon)
# load the example data
# devtools::load_all()
data_info <- read.csv(system.file("extdata", "data_info.csv", package = "Phydon"))
```

The data_info is a data frame that contains the following columns:

    #>                                               gene_loc       accession_no
    #> 1  gene_file/RS_GCF_002749895.1/RS_GCF_002749895.1.ffn RS_GCF_002749895.1
    #> 2  gene_file/RS_GCF_002849855.1/RS_GCF_002849855.1.ffn RS_GCF_002849855.1
    #> 3  gene_file/RS_GCF_002906655.1/RS_GCF_002906655.1.ffn RS_GCF_002906655.1
    #> 4  gene_file/RS_GCF_003026105.1/RS_GCF_003026105.1.ffn RS_GCF_003026105.1
    #> 5  gene_file/RS_GCF_003026475.1/RS_GCF_003026475.1.ffn RS_GCF_003026475.1
    #> 6  gene_file/RS_GCF_003026815.1/RS_GCF_003026815.1.ffn RS_GCF_003026815.1
    #> 7  gene_file/RS_GCF_003144035.1/RS_GCF_003144035.1.ffn RS_GCF_003144035.1
    #> 8  gene_file/RS_GCF_003544875.1/RS_GCF_003544875.1.ffn RS_GCF_003544875.1
    #> 9  gene_file/RS_GCF_003716875.1/RS_GCF_003716875.1.ffn RS_GCF_003716875.1
    #> 10 gene_file/RS_GCF_900130105.1/RS_GCF_900130105.1.ffn RS_GCF_900130105.1

The columns are:

- `gene_loc`: the location of the gene file (.ffn) in the package

- `accession_no`: the accession number of the genome

Now, you are ready to estimate the maximum growth rate.

``` r
# estimate the maximum growth rate
result <- Phydon(data_info)
#> [1] "Loading internal data ..."
#> [1] "10 genomes are found in the gene file and for estimation of the maximum growth rates ..."
```

    #>                genome         rep_genome                       species
    #> 1  RS_GCF_002749895.1 RS_GCF_002749895.1         s__Vibrio fujianensis
    #> 2  RS_GCF_002849855.1 RS_GCF_002849855.1             s__Vibrio azureus
    #> 3  RS_GCF_002906655.1 RS_GCF_002906655.1          s__Vibrio hyugaensis
    #> 4  RS_GCF_003026105.1 RS_GCF_003026105.1 s__Photobacterium sp003026105
    #> 5  RS_GCF_003026475.1 RS_GCF_003026475.1 s__Photobacterium lipolyticum
    #> 6  RS_GCF_003026815.1 RS_GCF_003026815.1    s__Photobacterium damselae
    #> 7  RS_GCF_003144035.1 RS_GCF_003144035.1               s__Vibrio albus
    #> 8  RS_GCF_003544875.1 RS_GCF_003544875.1         s__Vibrio alfacsensis
    #> 9  RS_GCF_003716875.1 RS_GCF_003716875.1              s__Vibrio zhugei
    #> 10 RS_GCF_900130105.1 RS_GCF_024346755.1           s__Vibrio aerogenes
    #>    neighbor_repgenome_train phy_distance gRodonpred phylopred combopred
    #> 1        RS_GCF_000621645.1   0.12220675  1.1346118 0.4697509 0.4697509
    #> 2        RS_GCF_900460535.1   0.05166449  0.4913568 0.2783264 0.2783264
    #> 3        RS_GCF_900460535.1   0.03447402  0.7658604 0.2783264 0.2783264
    #> 4        RS_GCF_003025615.1   0.16585578  1.5986925 2.0439385 2.0439385
    #> 5        RS_GCF_003025615.1   0.08080289  1.4924976 1.8805498 1.8805498
    #> 6        RS_GCF_003025615.1   0.15489320  0.7407403 2.0439385 2.0439385
    #> 7        RS_GCF_900460535.1   0.15592277  0.6097861 0.8919112 0.8919112
    #> 8        RS_GCF_900460535.1   0.03506859  0.6234453 0.2783264 0.2783264
    #> 9        RS_GCF_000621645.1   0.18306980  0.6364132 0.4697509 0.4697509
    #> 10       RS_GCF_900460535.1   0.17884206  1.6145463 0.5057288 0.5057288
