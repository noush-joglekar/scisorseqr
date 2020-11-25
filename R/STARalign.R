#' Align fastq.gz files with STARlong
#'
#' @aliases STAR
#' @description This function is a wrapper for STARlong
#' with stringent parameters optimized for PacBio sequenced
#' data. It performs sam to bam conversion to reduce
#' space, and .bam output is necessary for downstream analysis
#' @seealso \code{\link{MapAndFilter}}
#' @param fqFolder fastq.gz files from a single sample or replicate,
#' or barcoded output
#' @param starProgPath path to STARlong aligner
#' @param refGenomeIndex location of STAR reference genome index.
#' @param numThreads number of parallel threads to use, Defaults to 8
#' @return STARoutput Folder containing star output files and reports
#' @usage STARalign('FastqFiles/','~/STARlong','~/starIndex_gencode10_sequins/',32)
#' @export

STARalign <- function(fqFolder, starProgPath, refGenomeIndex, numThreads=8) {

  # Check that bedtools and samtools are correctly installed and loaded.
  checkFile <- system.file("bash", "toolCheck.sh", package = "scisorseqr")
  if(system(paste("sh", checkFile)) == 127) {
    stop(paste0("Error: samtools necessary for conversion to .bam format \n",
                "Check that both bedtools and samtools are installed and loaded before moving forward."))
  }

  starComm <- system.file("bash", "STARcomm.sh", package = "scisorseqr")
  miscFolder <- 'Misc/'
  starOut <- 'STARoutput/'

  if(!dir.exists(miscFolder)){dir.create(miscFolder)}
  if(!dir.exists(starOut)){dir.create(starOut)}

  print("Aligning with STAR")
  runIt <- paste("sh", starComm, fqFolder, starOut, starProgPath, refGenomeIndex, numThreads)
  system(runIt)

}
