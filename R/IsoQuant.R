#' Quantifies full-length isoforms in sample(s)
#'
#' @description Depending on whether you want to do downstream
#' differential expression analysis of full-length isoforms,
#' TSS, or PolyA sites, this function annotates and counts up
#' transcripts expressed in the given set of samples
#'
#' @param AllInfo Merged file of Gene, Barcode, UMI, Cluster,
#' and transcript across samples and replicates. Major isoform
#' is calculated based on the entire input.
#' @param Iso Logical indicating whether to conduct analysis for full-length
#' isoforms. Defaults to TRUE
#' @param TSS OPTIONAL logical indicating whether to conduct analysis for transcription
#' start sites (TSS). Defaults to FALSE
#' @param PolyA OPTIONAL logical indicating whether to conduct analysis for end sites
#' (PolyA). Defaults to FALSE
#' @return IsoQuantOutput folder containing all useful files for downstream
#' analysis
#'
#' @usage IsoQuant('AllInfo',Iso=TRUE, TSS=FALSE, PolyA=FALSE)
#' @seealso \code{\link{DiffSplicingAnalysis}}
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @export
IsoQuant <- function(AllInfoFile, Iso=TRUE, TSS=FALSE, PolyA=FALSE) {

  R_file <- system.file("RScript", "IsoQuant.R", package = "scisorseqr")
  sh_file <- system.file("bash","isoQuant_ds.sh",package = "scisorseqr")

  iso_quant_folder <- "IsoQuantOutput/"
  dir.create(iso_quant_folder)

  run_prog <- paste("Rscript", R_file, AllInfoFile, Iso, TSS, PolyA, iso_quant_folder)
  DS_iso <- paste("sh", sh_file, iso_quant_folder, "Gene_Cluster_Cell_IsoIDs_PerRead.csv", Iso, TSS, PolyA)

  system(run_prog)
  if(file.exists(file.path(iso_quant_folder, "Gene_Cluster_Cell_IsoIDs_PerRead.csv"))){
    system(DS_iso)}

}
