#' Align fastq.gz files with STARlong
#' @aliases STAR
#' @description This function is a wrapper for STARlong
#' with stringent parameters optimized for PacBio sequenced
#' data. It also performs sam to bam conversion to reduce
#' space and .bam output is necessary for downstream analysis
#' @seealso MapAndQuantify
#' @param fqFolder fastq.gz files from a single sample or replicate
#' @param starProgPath path to STARlong aligner
#' @param refGenome location of reference genome. Defaults to mouse mm10
#' @param numThreads number of parallel threads to use, Defaults to 8
#' @return STARoutput Folder containing star output files and reports
#' @usage STARalign('FastqFiles/',
#' '/athena/tilgnerlab/store/hut2006/soft/src/star-mapper/2.5.2b/STAR-2.5.2b/bin/Linux_x86_64/STARlong',
#' /athena/tilgnerlab/store/hut2006/data/seq/genomes/M.musculus/mm10/wholeGenomeUnzipped/starIndex_gencode10_sequins/,
#' 32)
#' @export

STARalign <- function(fqFolder, starProgPath, refGenome, numThreads=8) {

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


  runIt <- paste("sh", starComm, fqFolder, miscFolder, starOut, starProgPath, refGenome, numThreads)
  system(runIt)

}
