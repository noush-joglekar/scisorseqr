## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----makeBCfile, eval=FALSE, echo=TRUE----------------------------------------
#  
#  library(Seurat)
#  library(tidyr)
#  library(dplyr)
#  
#  bc <- data.frame("Barcode"=rownames(my_object@meta.data),
#                   "Celltype"=my_object@meta.data$celltype) %>%
#    separate(Barcode,into=c("BC","suffix"),sep="-") %>%
#    select(-suffix) %>% as.data.frame()
#  
#  write.table(bc, file = "PATH_TO_FILE", sep="\t", row.names = FALSE,
#              col.names = FALSE, quote = FALSE)

## ---- echo=FALSE--------------------------------------------------------------
bc_clust <- read.table(system.file("extdata/", "bc_celltype_assignments", package = "scisorseqr"), 
                       header = FALSE, sep = "\t")
knitr::kable(head(bc_clust))

## ----getBarcodes, eval=FALSE, echo=TRUE---------------------------------------
#  bc_clust <- system.file("extdata/", "bc_celltype_assignments", package = "scisorseqr")
#  fastqFolder <- system.file("extdata/", "FastqFiles", package = "scisorseqr")
#  
#  ## run command
#  GetBarcodes(fqFolder = fastqFolder, BCClustAssignFile = bc_clust,
#              chemistry = "v2", concatenate = TRUE, filterReads = FALSE)

## ----STARalign, eval = FALSE, echo = TRUE-------------------------------------
#  STARalign('FastqFiles/','~/bin/Linux_x86_64/STARlong',
#            refGenome='~/starIndex_gencode10_sequins/',
#            numThreads=24)

## ----MapAndFilter, eval = FALSE, echo = TRUE----------------------------------
#  
#  MapAndFilter(annoGZ='[PATH_TO_annotation.gtf.gz]',numThreads=24,
#               seqDir = '[PATH_TO_DIR/chr*.fa.gz]',
#               filterFullLength=TRUE,
#               cageBed = '[PATH_TO_Cage.bed.gz]',
#               polyABed = '[PATH_TO_PolyA.bed.gz]',
#               cp_distance=50)
#  

## ----InfoPerLongRead, eval=FALSE, echo=TRUE-----------------------------------
#  
#  InfoPerLongRead(barcodeOutputFile =
#                    'OutputFiltered/FilteredDeconvBC_AllFiles.csv',
#                  mapAndFilterOut = 'LRProcessingOutput/',
#                  minTimesIsoObserve = 3)
#  

## ----IsoQuant, echo = TRUE, eval= FALSE---------------------------------------
#  IsoQuant(AllInfoFile = 'LongReadInfo/AllInfo',
#            Iso = TRUE, TSS = TRUE, PolyA = TRUE)
#  

## ----ExonQuant, echo = TRUE, eval= FALSE--------------------------------------
#  ExonQuant(AllInfoFile = 'LongReadInfo/AllInfo',
#            groupingFactor = "Celltype", threshold = 10)

## ----DiffSplicing, eval=FALSE, echo=TRUE--------------------------------------
#  config <- system.file("extdata/", "config",
#                        package = "scisorseqr")
#  
#  DiffSplicingAnalysis(configFile = config,
#                       typeOfTest = 'Iso',
#                       minNumReads = 25, is.hier = FALSE)
#  
#  ## or for exons
#  DiffSplicingAnalysis(configFile = config,
#                       typeOfTest = 'Exon',
#                       minNumReads = 25,
#                       is.hier = FALSE)

## ----TriHeatmap, echo = TRUE, eval= FALSE-------------------------------------
#  triHeatmap(treeDir = 'TreeTraversal_Iso/',
#             comparisonList ='condensedCellTypes',
#             outName ="condensedCellTypes")

## ----heatmap, out.width = '60%', echo=F---------------------------------------
knitr::include_graphics("../man/figures/condensedCellTypes.png")

## ----PieChart, echo = TRUE, eval= FALSE---------------------------------------
#  sigSplitPie(compDir = 'Uniq_TreeTraversal_Iso/Astro_ExcitNeuron_10/')

## ----pie, out.width = '60%', echo=F-------------------------------------------
knitr::include_graphics("../man/figures/examplePie.png")

