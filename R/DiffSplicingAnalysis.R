#' Perform differential splicing analysis
#' @aliases DSA DiffTest DiffIsoTest
#' @description This function allows you to test alternative
#' splicing differences between any two groups
#' based on the full splice-form, exons,
#' the transcription start site (TSS) or end site (polyA).
#' @param config File in a .tsv format with four columns
#' specifying the comparisons to be made
#' @param numIsoforms Cap on the number of sites to test
#' per gene. Defaults to 10
#' @param minNumReads Parameter to filter out lowly expressed
#' genes. Defaults to 25
#' @param typeOfTest Exon, TSS, PolyA or Isoform test conducted
#' depending on input. Defaults to Iso
#' @return Output folder with results for each line of config
#' file, i.e. for each set of comparisons
#' @usage DiffSplicingAnalysis('config',10,25,'Exon')
#'
#' @export
DiffSplicingAnalysis <- function(configFile, numIsoforms = 10, minNumReads = 25, typeOfTest = "Iso") {
  R_file <- system.file("RScript", "DiffExp_Stretches_v00.R", package = "scisorseqr")
  call1 <- paste0("while read -r line; do echo $line | Rscript ", R_file,
                  " IsoQuantOutput/NumIsoPerCluster $(awk \'{split($0,a,\"\\t\"); {print a[1],a[2],a[3],a[4]}}\') ",
                  numIsoforms, " ", minNumReads," ", typeOfTest, " ; done < ", configFile)
  system(call1)
}
