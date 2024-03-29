% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MMalign.R
\name{MMalign}
\alias{MMalign}
\alias{Minimap2}
\title{Align fastq.gz files with Minimap2}
\usage{
MMalign('FastqFiles/','~/minimap2',
'~/mm10.fa',32)
}
\arguments{
\item{fqFolder}{fastq.gz files from a single sample or replicate,
or barcoded output}

\item{mmProgPath}{path to Minimap2 aligner}

\item{refGenome}{path to reference genome.fa}

\item{numThreads}{number of parallel threads to use, Defaults to 8}
}
\value{
MMoutput folder containing minimap2 aligned files and report
}
\description{
This function is a wrapper for Minimap2
with stringent parameters optimized for ONT/PacBio sequenced
data. It performs sam to bam conversion to reduce
space, and .bam output is necessary for downstream analysis
}
\seealso{
\code{\link{MapAndFilter}}
}
