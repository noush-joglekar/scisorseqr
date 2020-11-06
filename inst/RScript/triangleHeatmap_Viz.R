#' @import dplyr
#' @importFrom magrittr %>%

`%>%` <- magrittr::`%>%`

args <- commandArgs(trailingOnly = TRUE)
tree_dir <- args[1]

files = list.files(tree_dir,recursive = TRUE)[grep("results.csv",list.files(tree_dir,recursive = TRUE))]

fileNames = unlist(strsplit(files,"_10/"))[c(TRUE,FALSE)]

x <- as.vector(read.table(args[2])$V1)
cT <- unique(unlist(strsplit(x,"_"))[c(FALSE,TRUE)])

summaryMat <- matrix(0,length(cT),length(cT))
rownames(summaryMat) <- colnames(summaryMat) <- cT

for (file in files){
  i = grep(file,files)
  f = read.table(paste0(tree_dir,"/",file),header=TRUE,row.names = 1)
  perc <- nrow(f %>% dplyr::filter(FDR <= 0.05 & abs(dPI) >= 0.1))/nrow(f)
  s = unlist(strsplit(file,"_"))[1:2]
  summaryMat[s[1],s[2]] = perc
  summaryMat[s[2],s[1]] = perc
}


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(summaryMat)
mS <- reshape2::melt(upper_tri,na.rm = TRUE)


g = ggplot2::ggplot(data = mS, aes(x=Var1, y=Var2, fill=as.numeric(value)*100)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "inferno") +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()+labs(title=paste0("Percentage significant genes\n",args[2]),
                     x ="CellType1", y = "CellType2") +
  geom_text(aes(Var1, Var2, label = round(as.numeric(value)*100,digits = 1)), color = "black", size = 4)


pdf(paste0("Visualizations/",args[3],".pdf"),12,8,useDingbats=FALSE)
print(g)
dev.off()
