% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/triHeatmap.R
\name{triHeatmap}
\alias{triHeatmap}
\title{Creates a heatmap of DIE abundance across comparisons}
\usage{
triHeatmap('TreeTraversal_Iso/','neuronal_cellTypes','NeuronSubtype')
}
\arguments{
\item{treeDir}{output directory containing DiffSplicingAnalysis}

\item{comparisonList}{List of comparisons to be included in plot}

\item{typeOfTest}{type of test being visualized, e.g. "Exon" or "Iso".
Defaults to "Iso".}

\item{outName}{plot label}
}
\value{
Directory labelled 'Visualizations' containing a triangular heatmap
of pairwise comparisons containing cell groups provided in input
}
\description{
A quick way to summarize DIE across multiple
pairwise comparisons
}
\seealso{
\code{\link{DiffSplicingAnalysis}}
}
