#' Perform differential splicing analysis
#' @aliases DiffIsoTest
#' @description This function allows you to test alternative
#' splicing differences between any two groups
#' based on the full splice-form, exons,
#' the transcription start site (TSS) or end site (polyA).
#' @param configFile File in a .tsv format with four columns
#' specifying the comparisons to be made
#' @param numIsoforms Cap on the number of sites to test
#' per gene. Defaults to 10
#' @param minNumReads Parameter to filter out lowly expressed
#' genes. Defaults to 25
#' @param typeOfTest Exon, TSS, PolyA or Isoform test conducted
#' depending on input. Defaults to Iso
#' @param is.hier OPTIONAL logical indicating hierarchical structure
#' @param region1 string containing comparison1 for parsing if is.hier
#' is TRUE
#' @param region2 string containing comparison2 for parsing if is.hier
#' is TRUE
#' @param yamlFile .yaml file indicating hierarchical structure
#' @return Output folder with results for each line of config
#' file, i.e. for each set of comparisons
#' @return Result file with gene, uncorrected pval, deltaPI, FDR
#' @return handy .Robj per comparison to import if need be
#' @return file with significantly differentially spliced genes
#' @return file with numbers per gene and isoform tested in each group
#' @usage DiffSplicingAnalysis('config',10,25,'Exon',is.hier=FALSE)
#'
#' @seealso \code{\link{IsoQuant}}
#'
#' @export
DiffSplicingAnalysis <- function(configFile, numIsoforms = 10, minNumReads = 25, typeOfTest = "Iso",
                                 is.hier = FALSE, region1 = NULL, region2 = NULL, yamlFile = NULL) {
  R_file <- system.file("RScript", "IsoformTest.R", package = "scisorseqr")
  exonFile <- system.file("RScript", "ExonTest.R", package = "scisorseqr")

  if (typeOfTest == "Iso"){
    inputFile = "IsoQuantOutput/NumIsoPerCluster"
  } else if (typeOfTest == "TSS"){
    inputFile = "IsoQuantOutput/NumTSSPerCluster"
  } else if (typeOfTest == "PolyA"){
    inputFile = "IsoQuantOutput/NumPolyAPerCluster"
  } else if (typeOfTest == "Exon"){
    inputFile = "ExonQuantOutput/InclusionExclusionCounts.tsv"
  }

  if (typeOfTest != "Exon"){
    if (is.hier == FALSE){
      call1 <- paste("while read -r line; do echo $line | Rscript", R_file,
                     inputFile, "$(awk \'{split($0,a,\"\\t\"); {print a[1],a[2],a[3],a[4]}}\')",
                     numIsoforms, minNumReads, typeOfTest, FALSE, "; done <", configFile)
    } else {
      call1 <- paste("while read -r line; do echo $line | Rscript ", R_file,
                     inputFile, "$(awk \'{split($0,a,\"\\t\"); {print a[1],a[2],a[3],a[4]}}\')",
                     numIsoforms, minNumReads, typeOfTest, TRUE ,
                     yamlFile, region1, region2, "; done <", configFile)
    }
  } else {
    if (is.hier == FALSE){
      call1 <- paste("while read -r line; do echo $line | Rscript", exonFile,
                     inputFile, "$(awk \'{split($0,a,\"\\t\"); {print a[1],a[2],a[3],a[4]}}\')",
                     FALSE, "; done <", configFile)
    } else {
      call1 <- paste("while read -r line; do echo $line | Rscript ", exonFile,
                     inputFile, "$(awk \'{split($0,a,\"\\t\"); {print a[1],a[2],a[3],a[4]}}\')",
                     TRUE, numIsoforms, minNumReads, yamlFile,
                     region1, region2, "; done <", configFile)
    }
  }

  system(call1)
}
