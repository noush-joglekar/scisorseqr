% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InfoPerLongRead.R
\name{InfoPerLongRead}
\alias{InfoPerLongRead}
\alias{LongReadInfo}
\title{Combine all available information per read to a flat file}
\usage{
InfoPerLongRead(
  barcodeOutputFile,
  mapAndFilterOut,
  minTimesIsoObserve = 5,
  rmTmpFolder = TRUE
)
}
\arguments{
\item{barcodeOutputFile}{.csv file containing barcode and cell-type information
per read from the output of \code{\link{GetBarcodes}}}

\item{mapAndFilterOut}{output directory of the mapping function. If full-length
reads have been filtered using CAGE and PolyA site peaks, then it defaults
to that output, else it uses the canonically spliced full-length reads}

\item{minTimesIsoObserve}{minimum number of times an isoform is observed in the
dataset. Defaults to 5}

\item{rmTmpFolder}{Logical indicating whether you want to delete contents of the 
pre-processing folder. Defaults to TRUE}
}
\description{
Function to concatenate all the information per read, i.e
gene-name, cellular barcode, UMI, cell-type information, and isoform
information into one file. Isoform information includes a string separated list of 
introns, Cage and PolyA peak information if available, and a 
string separated list of exons. 
AllInfo contains the list of complete and full-length mapped, spliced, and barcoded reads.
AllInfo_Incomplete is a superset of the above, and also contains reads that are not
classified as complete based on Cage/PolyA data. 
The function also outputs basic stats 
such as number of reads/ genes/ UMIs per cellular barcode
}
\seealso{
\code{\link{MapAndFilter}}

\code{\link{GetBarcodes}}
}
