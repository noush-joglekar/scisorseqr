#' Conducts gene set enrichment analysis of significantly spliced genes
#' @description An easy way to look for enriched genes in the differential
#' splicing analysis. Uses the ClusterProfiler package to compare geneset
#' against the reference genome and outputs an .Robj file along with .pdf
#' network graphs
#' @seealso \url{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}
#' @param working_directory where the OutputGO folder should be created. Defaults to current
#' working directory
#' @param comparisonsDir location of TreeTraversal folder containing output
#' files from differential isoform expression analysis
#' @return OutputGO containing folders for Molecular Function, Biological
#' Process, and Cellular Component
#' @return .Robj files for MF, BP, and CC
#' @usage GO.Reference('./','TreeTraversal_Iso/')
#'
#' @export
GO.Reference <- function(working_directory = './',comparisonsDir = 'TreeTraversal_Iso/') {

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package \"clusterProfiler\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  go_folder <- 'OutputGO/'
  if(!dir.exists(paste0(working_directory,'/',go_folder))){
    dir.create(paste0(working_directory,'/',go_folder))
  }

  R_file <- system.file("RScript", "GOAnalysisWithPlots.R", package = "scisorseqr")
  run_prog <- paste("Rscript", R_file, working_directory, comparisonsDir)
  system(run_prog)
}
