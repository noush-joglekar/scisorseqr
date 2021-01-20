
<!-- README.md is generated from README.Rmd. Please edit that file -->

## scisorseqr - a comparative analysis of alternative splicing patterns

<img src="man/figures/scisorseqr.png" width="25%" style="float:left; padding:20px" />

scisorseqr is a linux based R-package for analyzing differential isoform
expression in single cells. The methods are based on our recent
[publication](https://www.nature.com/articles/s41467-020-20343-5) and
the [scISOrSeq
workflow](https://www.nature.com/articles/nbt.4259?draft=marketing)

Any comparative studies of alternative splicing can be performed with
scisorseqr. The package includes functions for barcode deconvolution
from fastqs, integration with long read alignment tools, mapping and
filtering of high confidence, full-length spliced reads, and some handy
tools to conduct differential expression analysis.

The tools are also applicable to long-read spatial transcriptomics, and
can be used to resolve exon expression at the
[spatial](https://isoformatlas.com/) level

-----

## Hardware / software requirements

The package has only been tested on a CentOS x86\_64 machine. For
alignment and mapping, we recommend

  - [STARlong](https://github.com/alexdobin/STAR/) software installation
    for PacBio reads
  - [Minimap2](https://github.com/lh3/minimap2) installation for Oxford
    Nanopore (or PacBio) reads
  - samtools
  - bedtools
  - python version 3.7

## Installation

The easiest way to install scisorseqr is through
[Github](https://github.com)
with:

``` r
devtools::install_github('noush-joglekar/scisorseqr',build_vignettes = TRUE)
```

## Workflow

<img src="man/figures/README-flow-2.png" width="60%" />

These steps are available as functions in the package. For example,
barcode deconvolution can be done using the following command

``` r
library(scisorseqr)
GetBarcodes('FastqFiles/','userInput/BarcodeCluster_Assignments', concatenate=TRUE, 
  filterReads=FALSE, numProcesses=24)
```

A step-by-step outline of the various functions is available as a
vignette. To access it, run

``` r
browseVignettes("scisorseqr")
```

## Support

We appreciate any and all inputs for improving scisorseqr. Feel free to
send us an [email](mailto:anj2026@med.cornell.edu).
