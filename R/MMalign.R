#' Align fastq.gz files with Minimap2
#' @aliases Minimap2
#' @description This function is a wrapper for Minimap2
#' with stringent parameters optimized for ONT/PacBio sequenced
#' data. It performs sam to bam conversion to reduce
#' space, and .bam output is necessary for downstream analysis
#' @seealso \code{\link{MapAndFilter}}
#' @param fqFolder fastq.gz files from a single sample or replicate,
#' or barcoded output
#' @param mmProgPath path to Minimap2 aligner
#' @param refGenome path to reference genome.fa
#' @param numThreads number of parallel threads to use, Defaults to 8
#' @return MMoutput folder containing minimap2 aligned files and report
#' @usage MMalign('FastqFiles/','~/minimap2',
#' '~/mm10.fa',32)
#' @export

MMalign <- function(fqFolder, mmProgPath, refGenome, numThreads=8) {

  # Check that bedtools and samtools are correctly installed and loaded.
  checkFile <- system.file("bash", "toolCheck.sh", package = "scisorseqr")
  if(system(paste("sh", checkFile)) == 127) {
    stop(paste0("Error: samtools necessary for conversion to .bam format \n",
                "Check that both bedtools and samtools are installed and loaded before moving forward."))
  }

  mmComm <- system.file("bash", "minimap2comm.sh", package = "scisorseqr")
  miscFolder <- 'Misc/'
  mmOut <- 'MMoutput/'

  if(!dir.exists(miscFolder)){dir.create(miscFolder)}
  if(!dir.exists(mmOut)){dir.create(mmOut)}

  print("Aligning with minimap2")
  runIt <- paste("sh", mmComm, fqFolder, mmOut, mmProgPath, refGenome, numThreads)
  system(runIt)

}
