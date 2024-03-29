#' Extract isoform information from aligned reads
#'
#' @description This function follows the alignment function.
#' It is a wrapper for several scripts that impose checks and
#' balances on the mapping quality, consensus nature of splice
#' sites, and annotated start and end sites per transcript.
#' @seealso \code{\link{STARalign}} \code{\link{MMalign}}
#' @param outputDir location within the working directory for files
#' that will be created during this process
#' @param annoGZ annotation.gtf.gz file. Defaults to v21 mm10
#' @param numThreads number of threads for parallel processes. Defaults to 12
#' @param seqDir directory containing chromosome.fa.gz files
#' @param filterFullLength OPTIONAL logical indicating whether reads should be filtered
#' for having start and end sites falling into annotated cage peaks or polyA
#' sites. Defaults to FALSE
#' @param cageBed OPTIONAL bed.gz file containing annotated CAGE peaks
#' @param polyABed OPTIONAL bed.gz file containing annotated polyA sites
#' @param cp_distance OPTIONAL distance from annotated CAGE peak or polyA site
#' for a read to be considered full-length. Defaults to 50
#' @param genomeVersion genome and version for GFF output file. Defaults to mm10
#' @param onlyFullLength OPTIONAL If you have already run this command without
#' annotated CAGE and PolyA sites and now want to just run that portion, set to
#' TRUE. Defaults to FALSE
#' @return directory containing several bam, gff.gz, and flat files
#' necessary for downstream analysis
#' @usage MapAndFilter('LRoutput','gencode.vM21.annotation.gtf.gz',16)
#' @export
#'
MapAndFilter <- function(outputDir = 'LRProcessingOutput/', annoGZ = NULL, numThreads = 12,
                        seqDir, filterFullLength = FALSE, cageBed = NULL, polyABed = NULL,
                        cp_distance = 50, genomeVersion = NULL, onlyFullLength = FALSE){

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

  fastqGuide <- paste0(miscFolder,'fastqGuide')
  bamGuide <- paste0(miscFolder,'bamGuide')

  tmpFolder <- 'tmpDir'

  if(is.null(annoGZ)){
    annoGZ <- system.file("extdata/", "gencode.vM21.annotation.gtf.gz", package = "scisorseqr")
  }

  if(is.null(genomeVersion)){
    genomeVersion <- "mm10"
  }

  awkScriptDir <- system.file("bash", package = "scisorseqr")
  pyScriptDir <- system.file("python", package = "scisorseqr")

  runCommand <- paste("sh", mainFile, fastqGuide, outputDir, tmpFolder, numThreads,
                      seqDir, annoGZ, bamGuide, awkScriptDir, pyScriptDir, genomeVersion)

  if(!dir.exists(outputDir)){dir.create(outputDir)}

  if(onlyFullLength == FALSE){
    if(!dir.exists(tmpFolder)){dir.create(tmpFolder)}
    system(runCommand)
  }

  if(filterFullLength == TRUE){
    cageGZ <- cageBed
    polyaGZ <- polyABed
    filterFL_script <- system.file("bash", "cagePolyA.sh", package = "scisorseqr")
    filterFLcomm <- paste("sh", filterFL_script, "zcat", outputDir,
                          cageGZ, polyaGZ, awkScriptDir, cp_distance)
    system(filterFLcomm)
  }

}
