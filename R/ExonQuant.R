#' Count exon inclusion and exclusion according to grouping factor
#'
#' @description This function groups reads by either the celltype
#' or cell-barcode to count inclusion and exclusion of exons. Inputs
#' to this function are based on previous preprocessing steps
#'
#' For ONT data we recomment \url{https://github.com/ablab/IsoQuant.git}
#' which allows for non-exact splice-site matching
#' @param allInfoFile file containing barcode, celltype, and exon information
#' per read. Defaults to output of the InfoPerLongRead function
#' @param groupingFactor "Celltype" or "Barcode" to group reads by for
#' counting inclusion levels. Defaults to Celltype
#' @param threshold minimum number of reads per grouping factor in order
#' to consider that exon to be sufficiently expressed. Defaults to 10
#' @param threads number of threads to parallelize the function. Defaults
#' to 4
#' @seealso \code{\link{MapAndFilter}}
#' @seealso \code{\link{InfoPerLongRead}}
#' @return ExonQuantOutput/InclusionExclusionCounts.tsv
#' 
#' @import dplyr
#' @import parallel
#' @export

ExonQuant <- function(allInfoFile = 'LongReadInfo/AllInfo_IncompleteReads.gz',
                           groupingFactor = "Celltype",threshold = 10, numThreads = 4) {

  if(!dir.exists('ExonQuantOutput')){dir.create("ExonQuantOutput/")}
  if(file.exists('ExonQuantOutput/InclusionExclusionCounts.tsv')){
    file.remove('ExonQuantOutput/InclusionExclusionCounts.tsv')}

  R_file <- system.file("RScript", "ExonCounting.R", package = "scisorseqr")

  countExons <- paste("Rscript", R_file, allInfoFile, groupingFactor, numThreads, threshold)
  system(countExons)

}
