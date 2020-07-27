#' Extract isoform information from aligned reads
#'
#' @description This function follows the STARalign function.
#' It is a wrapper for several scripts that impose checks and
#' balances on the mapping quality, consensus nature of splice
#' sites etc
#' @seealso STARalign
#' @param outputDir location within the working directory for files
#' that will be created during this process
#' @param annoGZ annotation.gtf.gz file. Defaults to v21 mm10
#' @param numThreads number of threads for parallel processes. Defaults to 12
#' @param filterFullLength logical indicating whether reads should be filtered
#' for having start and end sites falling into annotated cage peaks or polyA
#' sites. Defaults to FALSE
#' @param cageBed bed.gz file containing annotated CAGE peaks
#' @param polyABed bed.gz file containing annotated polyA sites
#' @param cp_distance distance from annotated CAGE peak or polyA site
#' for a read to be considered full-length. Defaults to 50
#' @return directory containing several bam, gff.gz, and flat files
#' necessary for downstream analysis
#' @usage MapAndAlign('LRoutput','gencode.vM21.annotation.gtf.gz',16)
#' @export
#'
MapAndAlign <- function(outputDir = 'LRProcessingOutput/', annoGZ, numThreads = 12,
                        filterFullLength = FALSE, cageBed = NULL, polyABed = NULL,
                        cp_distance = 50){

  checkFile <- system.file("bash", "toolCheck.sh", package = "scisorseqr")
  if(system(paste("sh", checkFile)) == 127) {
    stop(paste0("Error: samtools necessary for conversion to .bam format \n",
                "Check that both bedtools and samtools are installed and loaded before moving forward."))
  }

  if(filterFullLength == TRUE){
    if(is.null(cageBed) || is.null(polyABed)){
      stop("Please provide bed.gz files for CAGE and PolyA sites")
    } else if (grepl(".gz",cageBed) == FALSE || grepl(".gz",polyABed) == FALSE){
      stop("Please provide bed.gz files for CAGE and PolyA sites")
    }
  }

  miscFolder <- 'Misc/'
  mainFile <- system.file("bash", "mapAndAlignReads.sh", package = "scisorseqr")

  fastqGuide <- paste0(miscFolder,fastqGuide)
  bamGuide <- paste0(miscFolder,bamGuide)

  tmpFolder <- 'tmpDir'
  if(!dir.exists(tmpFolder)){dir.create(tmpFolder)}

  annoGZ <- system.file("extdata/", "gencode.vM21.annotation.gtf.gz", package = "scisorseqr")
  chromFaDir <- system.file("extdata/", "chromFa/", package = "scisorseqr")

  awkScriptDir <- system.file("bash", package = "scisorseqr")

  runCommand <- paste("sh", mainFile, fastqGuide, outputDir, tmpFolder, numThreads,
                      chromFaDir, annoGZ, bamGuide, awkScriptDir)

  if(filterFullLength == TRUE){
    cageGZ <- cageBed
    polyaGZ <- polyABed
    filterFL_script <- system.file("bash", "cagePolyA.sh", package = "scisorseqr")
    filterFLcomm <- paste("sh", filterFL_script, "zcat", outputDir,
                          cageGZ, polyaGZ, awkScriptDir, cp_distance)
  }

}
