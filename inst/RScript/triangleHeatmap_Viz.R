#' @import dplyr
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import cowplot

`%>%` <- magrittr::`%>%`

args <- commandArgs(trailingOnly = TRUE)
tree_dir <- args[1]
files = list.files(tree_dir,recursive = TRUE)[grep("results.csv",list.files(tree_dir,recursive = TRUE))]

cT <- unique(as.vector(read.table(args[2])$V1))

teType <- args[3]

summaryMat <- matrix(0,length(cT),length(cT))
rownames(summaryMat) <- colnames(summaryMat) <- cT

for (file in files){
  i = grep(file,files)
  f = read.table(file.path(tree_dir,file),header=TRUE,row.names = 1)
  if(teType == "Exon"){colnames(f)[3] <- "dPI"} 
  perc <- nrow(f %>% dplyr::filter(FDR <= 0.05 & abs(dPI) >= 0.1))/nrow(f)
  s = unlist(strsplit(basename(file),"_"))[1:2]
  if(s[1] %in% cT & s[2] %in% cT){
    summaryMat[s[1],s[2]] = perc
    summaryMat[s[2],s[1]] = perc
  }
}


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(summaryMat)
mS <- reshape2::melt(upper_tri,na.rm = TRUE)


g = ggplot2::ggplot(data = mS, ggplot2::aes(x=Var1, y=Var2, fill=as.numeric(value)*100)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_viridis_c(option = "inferno") +
  cowplot::theme_cowplot() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  ggplot2::coord_fixed()+ggplot2::labs(title=paste0("Percentage significant genes\n",args[4],"_",teType),
                     x ="CellType1", y = "CellType2") +
  ggplot2::geom_text(data = mS,
	ggplot2::aes(Var1, Var2, label = round(as.numeric(value)*100,digits = 1)),
	color = "black", size = 4)


pdf(paste0("Visualizations/",args[4],"_",teType,".pdf"),12,8,useDingbats=FALSE)
print(g)
dev.off()
