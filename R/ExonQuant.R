#' Count exon inclusion and exclusion according to grouping factor
#'
#' @description This function groups reads by either the celltype
#' or cell-barcode to count inclusion and exclusion of exons. Inputs
#' to this function are based on previous preprocessing steps
#' @param allInfoFile file containing cell barcode and celltype information
#' per read. Defaults to output of the InfoPerLongRead function
#' @param exon.gff gzipped file containing exon mapping information.
#' Defaults to output from MapAndFilter
#' @param groupingFactor "Celltype" or "Barcode" to group reads by for
#' counting inclusion levels. Defaults to celltype
#' @param threshold minimum number of reads per grouping factor in order
#' to consider that exon to be sufficiently expressed. Defaults to 10
#' @seealso \code{\link{MapAndFilter}} \code{\link{InfoPerLongRead}}
#' @return ExonQuantOutput/InclusionExclusionCounts.tsv
#' @export
#'

ExonQuant <- function(allInfoFile = 'LongReadInfo/AllInfo',
                           exon.gff = 'LRProcessingOutput/mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz',
                           groupingFactor = "Celltype",threshold = 10) {

  if(!dir.exists('ExonInfo')){dir.create("ExonInfo/")}

  R_file <- system.file("RScript", "ExonCounting.R", package = "scisorseqr")

  countExons <- paste("Rscript", allInfoFile, exon.gff, threshold, groupingFactor)
  system(countExons)

}
