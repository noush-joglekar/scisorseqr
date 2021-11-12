#' Creates a heatmap of DIE abundance across comparisons
#' @description A quick way to summarize DIE across multiple
#' pairwise comparisons
#' @seealso \code{\link{DiffSplicingAnalysis}}
#' @param treeDir output directory containing DiffSplicingAnalysis
#' @param comparisonList List of comparisons to be included in plot
#' @param typeOfTest type of test being visualized, e.g. "Exon" or "Iso".
#' Defaults to "Iso".
#' @param outName plot label
#' @return Directory labelled 'Visualizations' containing a triangular heatmap
#' of pairwise comparisons containing cell groups provided in input
#' @usage triHeatmap('TreeTraversal_Iso/','neuronal_cellTypes','NeuronSubtype')
#' @import dplyr
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import cowplot
#' @export

triHeatmap <- function(treeDir, comparisonList, typeOfTest = "Iso", outName = "triHeatmap") {

  `%>%` <- magrittr::`%>%`

  R_file <- system.file("RScript", "triangleHeatmap_Viz.R", package = "scisorseqr")

  outDir <- "Visualizations/"
  dir.create(outDir)

  run_prog <- paste("Rscript", R_file, treeDir, comparisonList, typeOfTest, outName)

  system(run_prog)
}
