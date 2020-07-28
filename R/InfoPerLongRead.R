#' Combine all information per read to a flat file
#' @aliases LongReadInfo
#' @description Function to concatenate all the information per read, i.e
#' gene-name, cellular barcode, UMI, cell-type information, and isoform
#' information into one file. Also outputs basic stats such as number of
#' reads/ genes/ UMIs per cellular barcode
#' @seealso \code{\link{MapAndFilter}}
#' @seealso \code{\link{GetBarcodes}}
#' @param barcodeOutput .csv file containing barcode and cell-type information
#' per read from the output of \code{\link{GetBarcodes}}
#' @param mapAndFilterOut output directory of the mapping function. If full-length
#' reads have been filtered using CAGE and PolyA site peaks, then it defaults
#' to that output, else it uses the canonically spliced full-length reads
#' @param minTimesIsoObserve minimum number of times an isoform is observed in the
#' dataset. Defaults to 5
#' @export
InfoPerLongRead <- function(barcodeOutput, mapAndFilterOut, minTimesIsoObserve = 5) {

  longReadInfo_sh <- system.file("bash", "longReadInfo.sh", package = "scisorseqr")

  longReadInfoFolder <- "LongReadInfo/"
  dir.create(longReadInfoFolder)

  if(file.exists(file.path(mapAndFilterOut,'newIsoforms_vs_Anno_ignoreAnno/CagePolyA.complete.stretches.gz'))){
    stretchesFile <- file.path(mapAndFilterOut,'newIsoforms_vs_Anno_ignoreAnno/CagePolyA.complete.stretches.gz')
  } else {file.path(mapAndFilterOut,'newIsoforms_vs_Anno_ignoreAnno/stretches.gz')}

  if(file.exists(file.path(mapAndFilterOut,'CagePolyA.complete.mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz'))){
    geneFile <- file.path(mapAndFilterOut,'CagePolyA.complete.mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz')
  } else {geneFile <- file.path(mapAndFilterOut,'mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz')}

  longReadInfoComm <- paste(longReadInfo_sh, perfectMatchFile, geneFile, stretchesFile, minTimesIsoObserve)
  system(longReadInfoComm)

}
