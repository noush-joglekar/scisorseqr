## By Anoushka Joglekar 2019. Edited 04/2020, 07/2020
#' @import dplyr
#' @importFrom magrittr %>%

`%>%` <- magrittr::`%>%`

args <- commandArgs(trailingOnly=TRUE)

### Initial set-up
numIsoPerCluster <- data.table::fread(args[1])	## Need to make the package seamless
					# enough that it takes the correct
					#file depending on iso,tss,polyA

comps <- c(args[2],args[4])

type1 <- unlist(strsplit(args[3],","))
type2 <- unlist(strsplit(args[5],","))

numIsoPerCluster$V2[numIsoPerCluster$V2 %in% type1] = "Group1"
numIsoPerCluster$V2[numIsoPerCluster$V2 %in% type2] = "Group2"
comparisons <- c("Group1","Group2")

inum <- as.integer(args[6])
threshold <- as.integer(args[7])

typeOfTest <- args[8]
is.hier <- args[9]

if(is.hier == TRUE){
  if(length(args) < 12){
    print("Please input path to hierarchy.yaml file")
    print("Indicate names of case-control corresponding to config file")
    print("Aborting")
    stop()
  } else {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Packages \"yaml\" and \"data.tree\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

    h <- yaml::yaml.load_file(args[10])
    hier <- data.tree::as.Node(h)
    region <- c(args[11],args[12])
    if(!dir.exists(paste0('TreeTraversal_Hier_',typeOfTest))){
      dir.create(paste0('TreeTraversal_Hier_',typeOfTest))
      out_dir <- paste0('TreeTraversal_Hier_',typeOfTest,'/',comps[1],"_",comps[2],"_",inum,"/")
      dir.create(out_dir)
    } else {
      out_dir <- paste0('TreeTraversal_Hier_',typeOfTest,'/',comps[1],"_",comps[2],"_",inum,"/")
      dir.create(out_dir)
    }
  }
} else {
  if(!dir.exists(paste0('TreeTraversal_',typeOfTest))){
    dir.create(paste0('TreeTraversal_',typeOfTest))
    out_dir <- paste0('TreeTraversal_',typeOfTest,'/',comps[1],"_",comps[2],"_",inum,"/")
    dir.create(out_dir)
    } else {
      out_dir <- paste0('TreeTraversal_',typeOfTest,'/',comps[1],"_",comps[2],"_",inum,"/")
      dir.create(out_dir)
    }
}

sink(paste0(out_dir,"REPORT_",threshold,".txt"), append=TRUE, split=TRUE)
Sys.time()
cat("Starting analysis in",comps[1],"vs",comps[2],"with atleast", threshold,"reads per gene \n")


############################
####### Functions ##########

deltaPI <- function(mat,gene){
  pos <- sum(c(rev(mat$delta[mat$delta > 0]),0,0)[1:2])
  neg <- abs(sum(c(mat$delta[mat$delta < 0],0,0)[1:2]))
  change=max(pos,neg)
  if (pos >= neg){
    x = c(rev(mat$delta[mat$delta > 0]),0)[1:2]
    y = mat %>% dplyr::filter(delta %in% x) %>% dplyr::select(IsoID) %>% as.matrix()
    return(list(change,y[1],y[2]))}
  else {
    x = c(mat$delta[mat$delta < 0],0,0)[1:2]
    y = mat %>% dplyr::filter(delta %in% x) %>% dplyr::select(IsoID) %>% as.matrix()
    return(list(-change,y[1],y[2]))}
}

Get_Pval_DeltaPI <- function(mat){
  max_ix1 <- deltaPI(mat,gene)[2]
  max_ix2 <- deltaPI(mat,gene)[3]
  pval <- chisq.test(mat[,3:4])$p.value
  d_pi <- deltaPI(mat,gene)[1]
	return(list(pval,d_pi,max_ix1,max_ix2))
}

############################

coresToUse <- parallel::detectCores()/2


numIsoPerCluster <- numIsoPerCluster[numIsoPerCluster$V2 %in% comparisons]
colnames(numIsoPerCluster) <- c("Isoform","CellType","NumTranscripts")
data_df <- tidyr::separate(numIsoPerCluster, Isoform, into=c("Gene","IsoID"), sep="_")

rm(numIsoPerCluster)

##### Check if hierarchical analysis is true
if(is.hier == TRUE){
  cellGroup <- unlist(strsplit(comps[1],region[1]))[2]
  hierLevel <- FindNode(hier, cellGroup)$level
  if (cellGroup != hier$root$name){
    parentCG <- data.tree::FindNode(hier, cellGroup)$parent$name
    sig <- read.table(paste0("../TreeTraversal_Hier_Iso/",
                             region[1],parentCG,"_",region[2],parentCG,"_",inum,"/",
                             'sigGenes_',region[1],parentCG,"_",region[2],parentCG,threshold),header=TRUE)
    data_df <- data_df[data_df$Gene %in% sig$x, ]
  }
} else {hierLevel <- 1}


data_df$IsoID <- as.integer(data_df$IsoID)
data_df$IsoID[data_df$IsoID >= inum] <- inum


processedDF = data_df %>% dplyr::group_by(Gene, CellType, IsoID) %>%
  dplyr::summarise(Sum = sum(NumTranscripts), .groups='drop' ) %>%
  dplyr::filter(sum(Sum) >= threshold & length(Sum) >= 2) %>%
  dplyr::ungroup() %>% dplyr::group_by(Gene) %>%
  dplyr::filter(length(unique(CellType)) == 2) %>%
	tidyr::spread(CellType,Sum) %>% replace(is.na(.), 0) %>%
  dplyr::mutate(pi1 = Group1/sum(Group1),
	pi2 = Group2/sum(Group2), delta = pi1-pi2) %>%
	dplyr::arrange(delta, .by_group = TRUE) %>% as.data.frame()

PerGene = split(processedDF,processedDF$Gene)

cat("About to start calculating p-values for",length(PerGene),"genes \n")



uncorrected_output <- parallel::mclapply(names(PerGene),
	function(geneName) Get_Pval_DeltaPI(PerGene[[geneName]]),mc.cores=coresToUse)
names(uncorrected_output) <- names(PerGene)


output_DF <- as.data.frame(do.call(rbind, uncorrected_output))
colnames(output_DF) <- c("pvals","dPI","maxDeltaPI_ix1","maxDeltaPI_ix2")
output_DF$FDR <- p.adjust(output_DF$pvals, method = "BH")*hierLevel
output_DF$dPI <- as.numeric(output_DF$dPI)

nums <- processedDF %>% dplyr::select(Gene, IsoID, Group1, Group2) %>% dplyr::group_by(Gene) %>%
	dplyr::arrange(IsoID, .by_group=TRUE) %>% as.data.frame()
colnames(nums)[3:4] <- comps

cat("Total genes tested:",nrow(output_DF),"\n")
cat("Number of significant genes by fdr:",length(which(output_DF$FDR <= 0.05)),"\n")
cat("Number of fdr genes with deltaPI >= 0.1:",length(which(output_DF$FDR <= 0.05 & abs(output_DF$dPI) >= 0.1)),"\n")
cat("Breakdown of isoform of sig genes with biggest change \n")
sigs <- output_DF[which(output_DF$FDR <= 0.05 & abs(output_DF$dPI) >= 0.1),]
print(table(unname(unlist(sigs$maxDeltaPI_ix1))))
cat('\n')
Sys.time()
sink()


print("Writing output files")
output_DF <- output_DF %>% tibble::rownames_to_column("Gene") %>% as.data.frame()
output_DF <- as.matrix(output_DF)

write.table(output_DF,
            file=paste0(out_dir,comps[1],'_',comps[2],'_',threshold,'_results.csv'),
            sep="\t",row.names = FALSE, quote = FALSE)

write.table(nums,file = paste0(out_dir,"NumsForGenesBeingTested_",comps[1],"_",comps[2],"_",threshold),
            quote = FALSE,row.names = FALSE)

colnames(processedDF)[3:4] <- comps
save(processedDF,file=
       paste0(out_dir,comps[1],'_',comps[2],'_',inum,'X2_Atleast',threshold,'reads.Robj'))


sig_Genes <- rownames(sigs)
write.table(sig_Genes,file = paste0(out_dir,'sigGenes_',comps[1],'_',comps[2],threshold),quote=FALSE,row.names = FALSE)
