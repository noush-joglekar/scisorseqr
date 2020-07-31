#' Generate matrices for clustering analysis
#'
#' @description Converts the counts for isoforms
#' into matrices clustered either by cell-type
#' or by barcode.
#' @seealso \code{\link{IsoQuant}}
#' @param isoQuantOutDir output directory of IsoQuant
#' function. Defaults to IsoQuantOutput
#' @param groupBy String indicating group for matrix output.
#' Current options involve "Cell" or "Celltype", and defaults
#' to "Cell"
#' @param ensemblToClear A mapping of gene-names from
#' ENSEMBL ids to human readable gene names for downstream
#' analysis. Defaults to mm10 gencode v21 annotation.
#' @param convertTo10XoutputFormat OPTIONAL If grouped by cell,
#' logical inputt indicating whether
#' or not matrix should be converted into a format similar to 10X
#' cellranger output for downstream analysis. Defaults to TRUE
#'
#' @export
MakeMatrices <- function(isoQuantOutDir = 'IsoQuantOutput/', groupBy = "Cell",
                         ensemblToClear = NULL, convertTo10XoutputFormat = TRUE) {

  usethis::use_package('utils')
  usethis::use_package('plyr')
  usethis::use_package('tidyr')
  usethis::use_package('methods')
  usethis::use_package('DropletUtils','Suggests')


  py_file1 <- system.file("python", "AllIsoformsXCell_sparseMatrix.py", package = "scisorseqr")
  py_file2 <- system.file("python", "Get_IsoformXClusterMat.py", package = "scisorseqr")

  matricesFolder <- "Matrices/"
  if(!dir.exists(matricesFolder)){dir.create(matricesFolder)}

  if(groupBy == "Cell"){
    isoXcell <- file.path(isoQuantOutDir,'IsoXNumInCell')
    run_cell <- paste("python3", py_file1, isoXcell)
    system(run_cell)

    if(convertTo10XoutputFormat == TRUE){

      if (!requireNamespace("DropletUtils", quietly = TRUE)) {
        stop("Package \"DropletUtils\" needed for this function to work. Please install it.",
             call. = FALSE)
      }

      isoMat <- utils::read.table("Matrices/AllIsoformsXAllCells_matrix.csv")
      geneNames <- data.frame("gene"=unlist(strsplit(rownames(isoMat),"_"))[c(TRUE,FALSE)])
      isoID <- as.vector(unlist(strsplit(rownames(isoMat),"_"))[c(FALSE,TRUE)])

      if(is.null(ensemblToClear)){
        ensemblToClear <- system.file("extdata", "ENSEMBLE_v21_ClearGeneNames", package = "scisorseqr")}

      clearGenes <- as.data.frame(utils::read.table(ensemblToClear,sep="\t"))
      colnames(clearGenes) <- c("gene","clear")

      isoClear <- as.vector(plyr::join(geneNames,clearGenes)$clear)
      isoClearID <- data.frame("clear"=isoClear,"id"=isoID)
      isoClearwithID <- tidyr::unite(isoClearID, "clearID", sep = "_")

      DropletUtils::write10xCounts("Matrices/SeuratFriendly/",
                                   methods::as(as.matrix(isoMat),"dgCMatrix"),
                                   barcodes=colnames(isoMat),
                                   gene.id=rownames(isoMat),
                                   gene.symbol=isoClearwithID$clearID)
    }

  } else if (groupBy == "Celltype") {
    isoXcluster <- file.path(isoQuantOutDir,'NumIsoPerCluster')
    run_clust <- paste("python3", py_file2, isoXcluster)
    system(run_clust)
  }
}
