#' Filter unbarcoded reads
#'
#' @description This function filters out all the reads that don't have a barcode
#' and outputs them into a separate folder.
#' Option to concatenate output from several files if they belong to
#' the same sample.
#' @param outputFolder automatically sourced from GetBarcodes script
#' @param raw_Output unfiltered output of GetBarcodes script
#' @param concatenate Combines output of all fastq.gz files if they belong
#' to the same sample. Defaults to TRUE
#' @export

FilterBCOutput <- function(outputFolder, raw_Output,  concatenate, filterReads, fqFolder) {
  filterOut <- system.file("bash", "FilterBCReads.sh", package = "scisorseqr")

  filt_Output <- paste0(outputFolder, "/OutputFiltered/")
  dir.create(filt_Output)

  # Maybe can add feature to randomize selection of one read per barcode-UMI pair
  # code in KLa's folder
  if (concatenate == TRUE){
    filterAndConcat <- paste("sh", filterOut, raw_Output, filt_Output, TRUE, filterReads, fqFolder)
    system(filterAndConcat)
  } else {
    filterOnly <- paste("sh", filterOut, raw_Output, filt_Output, FALSE, filterReads, fqFolder)
    system(filterOnly)
  }
}
