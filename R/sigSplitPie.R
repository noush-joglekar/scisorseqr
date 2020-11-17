#' Creates a piechart of number of isoforms needed to cross deltaPI threshold
#' @description A quick way to summarize DIE for a single comparison
#' @seealso \code{\link{DiffSplicingAnalysis}}
#' @param compDir output directory containing DiffSplicingAnalysis for a
#' specific comparison
#' @return Directory labelled 'Visualizations' containing a pie chart containing
#' the breakup of number of isoforms needed to cross deltaPI threshold
#' of pairwise comparisons containing cell groups provided in input
#' @usage sigSplitPie('compDir/')
#'
#' @import dplyr
#' @import ggplot2
#' @import scales
#' @importFrom magrittr %>%
#'
#' @export

sigSplitPie <- function(compDir) {

  `%>%` <- magrittr::`%>%`

  R_file <- system.file("RScript", "PieChart_Viz.R", package = "scisorseqr")

  outDir <- "Visualizations/"
  dir.create(outDir)

  run_prog <- paste("Rscript", R_file, compDir)

  system(run_prog)
}
