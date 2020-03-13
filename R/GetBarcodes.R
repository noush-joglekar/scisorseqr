#' Get single-cell barcodes from long read files
#'
#' This is an R-wrapper for a python function. It parallelizes the reading
#' and processing of fastq.gz files, and uses a sliding window approach to
#' identify cell barcodes, assign cluster labels and output useful statistics
#' in one .csv file and one _summary file per fastq.gz input file
#' @param fqFolder Path to a folder containing fastq.gz files
#' @param BCClustAssignFile tab separated file containing barcode in column 1
#' and cluster in column 2
#' @param outputFolder OPTIONAL Path to working directory. Defaults to current
#' @param numProcesses Number of chunks to split processing into. Defaults to 10
#' @param concatenate OPTIONAL Concatenate output from all files in the folder
#' @return OutputRaw folder containing a .csv file and summary stats for each
#' input fastq.gz file
#' @return OutputFiltered Filtered file with one line for each read containing
#' a barcode. Option to concatenate into one file per input folder
#' @usage GetBarcodes("FastqFolder","Barcode-Clust","~/MyDir/",6,concatenate = FALSE)
#' @export
GetBarcodes <- function(fqFolder, BCClustAssignFile, outputFolder = ".", numProcesses = 10, concatenate = TRUE) {
  if(!dir.exists(outputFolder))
    dir.create(outputFolder)
  # Checks for working directory as stated in input, by default will use
  # current directory
  raw_Output <- paste0(outputFolder, "/OutputRaw/")
  # Python and bash files are included in package
  py_file <- system.file("python", "polyA_barcode_UMI_detection.py", package = "scisorseqr")
  simConcat <- system.file("bash", "concat_ParallelizedBC.sh", package = "scisorseqr")
  filterOut <- system.file("bash", "FilterBCReads.sh", package = "scisorseqr")

  # Creates subdirectory for unfiltered output and summary files
  dir.create(raw_Output)
  files <- list.files(fqFolder)[grep(pattern = "q.gz",list.files(fqFolder))]

  # Iterates through each file, activating necessary scripts for
  # analysis and concatenation of parallel processes
  for (file in files) {
    input_file <- paste0(fqFolder, file)
    f_name <- unlist(strsplit(file,".fastq|.fq"))[1]
    run_bc_deConv <- paste("python3", py_file, input_file, BCClustAssignFile,
                           "--outDir", raw_Output, "--numProc", numProcesses)
    run_concatSimProc <- paste("bash", simConcat, f_name, raw_Output)
    system(run_bc_deConv)
    system(run_concatSimProc)
  }

  # Calls next function automatically once all files have been analyzed and
  # .csv files created.
  MaptoCellType(outputFolder, raw_Output, concatenate = TRUE)
}

#' @export
#' This function filters out all the reads that don't have a barcode
#' and outputs them into a separate folder.
#' Option to concatenate output from several files if they belong to
#' the same sample.
MaptoCellType <- function(outputFolder, raw_Output,  concatenate = TRUE) {
  filt_Output <- paste0(outputFolder, "/OutputFiltered/")
  dir.create(filt_Output)

  # Maybe can add feature to randomize selection of one read per barcode-UMI pair
  # code in KLa's folder
  if (concatenate == "TRUE"){
    filterAndConcat <- paste("sh", filterOut, raw_Output, filt_Output, TRUE)
    system(filterAndConcat)
  } else {
    filterOnly <- paste("sh", filterOut, raw_Output, filt_Output)
    system(filterOnly)
  }
}
