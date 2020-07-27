#' Align fastq.gz files with Minimap2
#' @aliases Minimap2
#' @description This function is a wrapper for Minimap2
#' with stringent parameters optimized for ONT/PacBio sequenced
#' data. It performs sam to bam conversion to reduce
#' space, and .bam output is necessary for downstream analysis
#' @seealso MapAndAlign
#' @param fqFolder fastq.gz files from a single sample or replicate or barcoded output
#' @param mmProgPath path to Minimap2 aligner
#' @param annoGZ annotation.gtf.gz file. Defaults to v21 mm10
#' @param numThreads number of parallel threads to use, Defaults to 8
#' @return MMoutput folder containing minimap2 aligned files and report
#' @usage MMalign('FastqFiles/','~/minimap2',
#' '/athena/tilgnerlab/store/hut2006/data/annotations/M.Musculus/mm10/gencode.vM21.annotation.gtf.gz',
#' 32)
#' @export

MMalign <- function(fqFolder, mmProgPath, annoGZ, numThreads=8) {

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


  runIt <- paste("sh", mmComm, fqFolder, mmOut, mmProgPath, annoGZ, numThreads)
  system(runIt)

}
