#' By Anoushka Joglekar 2019. Edited 04/2020
#' Added layer to regular exon level analysis script where hierarchical testing is allowed.
#' Input args: 1. file.tsv with Exon-Gene-Inclusion-Exclusion counts
#' args 2,3,4,5 (tsv from each line of config file):
#' --- Shorthand comparison1 of interest,
#' --- comma separated cell-types included,
#' --- Shorthand comparison2 of interest,
#' --- comma separated cell-types included

devtools::use_package('dplyr')
devtools::use_package('data.table')

args <- commandArgs(trailingOnly=TRUE)


inclusionFile <- data.table::fread(args[1])

comps <- c(args[2],args[4])

type1 <- unlist(strsplit(args[3],","))
type2 <- unlist(strsplit(args[5],","))

inclusionFile$V3[inclusionFile$V3 %in% type1] = args[2]	## rename cell-subtype to comparison name
inclusionFile$V3[inclusionFile$V3 %in% type2] = args[4]
comparisons <- comps

is.hier <- args[6]

if(is.hier == TRUE){
  if(length(args) < 11){
    print("Please input path to hierarchy.yaml file")
    print("Indicate names of case-control corresponding to config file")
    print("Aborting")
    stop()
  } else {
    devtools::use_package('yaml')
    inum <- as.integer(args[7])
    threshold = as.integer(args[8])
    h <- yaml::yaml.load_file(args[9])
    hier <- as.Node(h)
    region <- c(args[10],args[11])
    if(!dir.exists('TreeTraversal_Hier_Exon')){
      dir.create('TreeTraversal_Hier_Exon/')
      out_dir <- paste0('TreeTraversal_Hier_Exon/',comps[1],"_",comps[2],"/")
      dir.create(out_dir)
    } else {
      out_dir <- paste0('TreeTraversal_Hier_Exon/',comps[1],"_",comps[2],"/")
      dir.create(out_dir)}
  }
} else {
  if(!dir.exists('TreeTraversal_Exon')){
    dir.create('TreeTraversal_Exon/')
    out_dir <- paste0('TreeTraversal_Exon/',comps[1],"_",comps[2],"/")
    dir.create(out_dir)
  } else {
    out_dir <- paste0('TreeTraversal_Exon/',comps[1],"_",comps[2],"/")
    dir.create(out_dir)}
}



inclusionFile <- inclusionFile[inclusionFile$V3 %in% comparisons]
colnames(inclusionFile) <- c("Exon","Gene","Celltype","inclusion","exclusion")

##### Check if hierarchical analysis is true
if(is.hier == TRUE){
	cellGroup <- unlist(strsplit(comps[1],region[1]))[2]
	if (cellGroup != hier$root$name){
		parentCG <- FindNode(hier, cellGroup)$parent$name
		sig <- read.table(paste0("TreeTraversal_Hier_Iso/",
                        region[1],parentCG,"_",region[2],parentCG,"_",inum,"/",
			'sigGenes_',region[1],parentCG,"_",region[2],parentCG,threshold),header=TRUE)
		inclusionFile <- inclusionFile[inclusionFile$Gene %in% sig$x, ]
	}
}


sink(paste0(out_dir,"REPORT.txt"), append=TRUE, split=TRUE)
Sys.time()
cat("Starting exon level analysis in",comps[1],"vs",comps[2],"\n")


### Group by same cell-type and add counts (eg: hippEN has 3 contributing cell-subtypes)
condensedDF <- as.data.frame(inclusionFile %>% group_by(Exon,Gene,Celltype) %>%
                               dplyr::summarise(Inclusion = sum(inclusion),Exclusion= sum(exclusion)) %>%
                               dplyr::group_by(Exon) %>% dplyr::filter( n() > 1 ))

numExons <- dim(condensedDF)[1]/2

pval <- geneL <- exonL <- dpsi <- psi1 <- psi2 <- c()

## dataframe is sorted so two lines correspond to our comparisons of interest
## Filter on counts: chi-sq test criterion + PSI of matrix should be >=0.1 and <=0.9
checkAndCompute <- function(inputMat,ix){
	mat<- as.matrix(inputMat[c((2*ix)-1,2*ix),c(4,5)])
	if (sum(mat) > 0 && max(colSums(mat))/sum(mat) <= 0.9 &&
	min(rowSums(mat))*min(colSums(mat))/sum(mat) > 5){
		exonL <- c(exonL,as.character(inputMat[2*ix,1]))
		geneL <- c(geneL,as.character(inputMat[2*ix,2]))
		psis <- mat[,1]/rowSums(mat)
		psi1 <- c(psi1,psis[1])
		psi2 <- c(psi2,psis[2])
		pval <- c(pval,chisq.test(mat)$p.value)
		dpsi <- psis[1]-psis[2]
	}
	return(list(exonL,geneL,pval,dpsi,psi1,psi2))}

## Run function for the full dataset and convert to dataframe
res <- parallel::mclapply(1:numExons, function(ix) checkAndCompute(condensedDF,ix),mc.cores=24)
res <- unlist(plyr::compact(res))
res <- as.data.frame(matrix(res, ncol = 6,  byrow = TRUE), stringsAsFactors = FALSE)
colnames(res) <- c("Exon","Gene","Pval","dPSI","psi1","psi2")

res$FDR <- p.adjust(res$Pval,method="BY")


## summarize results for REPORT
sigCor <- which(res$FDR <= 0.05 & abs(as.numeric(res$dPSI)) >= 0.1)
sig <- which(res$FDR <= 0.05)
dPSI <- which(abs(as.numeric(res$dPSI)) >= 0.1)
inputList <- nrow(inclusionFile %>% dplyr::group_by(Gene) %>%
                    dplyr::do(data.frame(nrow=nrow(.))))
genesTested <- length(unique(res$Gene))
genesSig <- length(unique(res[sigCor,"Gene"]))

cat("Total exons tested:",length(res$FDR),"\n")
cat("Number exons with |deltaPSI| >= 0.1:",length(dPSI),"\n")
cat("Number of significant exons by BY:",length(sig),"\n")
cat("Number of BY sig genes with deltaPSI >= 0.1:",length(sigCor),"\n")
cat('\n')
cat("Input geneset:",inputList,"\n")
cat("Total genes tested:",genesTested,"\n")
cat("Significantly differentially spliced genes:",genesSig,"\n")
cat('\n')
Sys.time()
sink()

testedDF <- condensedDF[condensedDF$Exon %in% res$Exon,]
write.table(testedDF,file=paste0(out_dir,"TestedExons_",comps[1],"_",comps[2],".csv"),
	sep="\t",quote=FALSE,row.names=FALSE)

write.table(res,
	file=paste0(out_dir,comps[1],'_',comps[2],'_results.csv'),
	sep="\t",quote = FALSE,row.names=FALSE)
write.table(unique(res[sigCor,"Gene"]),file = paste0(out_dir,'sigGenes_',comps[1],'_',comps[2]),
	quote=FALSE,row.names = FALSE)

