## By Anoushka Joglekar 2019. Edited 07/2020

devtools::use_package('clusterProfiler','Suggests')
devtools::use_package('enrichplot','Suggests')

args <- commandArgs(trailingOnly=TRUE)
current_dir <- args[1]
setwd(current_dir)

pairwise_dir <- args[2]

files = list.files(pairwise_dir,recursive = TRUE)[grep(
  pattern = "results.csv",list.files(pairwise_dir,recursive = TRUE))]


mf <- list()
bp <- list()
cc <- list()

dir.create('OutputGO/GO_BP/')
dir.create('OutputGO/GO_MF/')
dir.create('OutputGO/GO_CC/')

for (file in files){
  name = strsplit(file,"_10")[[1]][1]
  print(paste("Conducting GO analysis against",name))
  result_file <- read.table(paste0(pairwise_dir,file),sep="\t", header=TRUE)
  sig_ix <- which(result_file$FDR <= 0.05 & abs(result_file$dPI)>=0.1)
  sigGenes <- unlist(strsplit(as.vector(result_file$Gene[sig_ix]),"[.]"))[c(TRUE,FALSE)]
  ego_MF <- clusterProfiler::enrichGO(sigGenes,OrgDb = org.Mm.eg.db,keyType= 'ENSEMBL',ont = "MF")
  ego_BP <- clusterProfiler::enrichGO(sigGenes,OrgDb = org.Mm.eg.db,keyType= 'ENSEMBL',ont = "BP")
  ego_CC <- clusterProfiler::enrichGO(sigGenes,OrgDb = org.Mm.eg.db,keyType= 'ENSEMBL',ont = "CCC")

  mf[name] <- ego_MF
  bp[name] <- ego_BP
  cc[name] <- ego_CC

  rm(ego_MF)
  rm(ego_BP)
  rm(ego_CC)
}

printMF <- function(mf,name,file){
  mf1 <- enrichplot::emapplot(mf[[name]], vertex.label.cex=1.2)

  pdf(file = paste0("OutputGO/GO_MF/",name,".pdf"),12,10)
  print(mf1)
  dev.off()
}

printBP <- function(bp,name,file){
  bp1 <- enrichplot::emapplot(bp[[name]], vertex.label.cex=1.2)

  pdf(file = paste0("OutputGO/GO_BP/",name,".pdf"),12,10)
  print(bp1)
  dev.off()
}

printBP <- function(cc,name,file){
  cc1 <- enrichplot::emapplot(cc[[name]], vertex.label.cex=1.2)

  pdf(file = paste0("OutputGO/GO_CC/",name,".pdf"),12,10)
  print(cc1)
  dev.off()
}

for (file in files){
  name = strsplit(file,"_10")[[1]][1]
  print(name)
  try(printMF(mf,name,file))
  try(printBP(bp,name,file))
  try(printCC(cc,name,file))

}

save(mf,file='OutputGO/GO_MF/MolecularFunction.Robj')
save(bp,file='OutputGO/GO_BP/BiologicalProcess.Robj')
save(cc,file='OutputGO/GO_CC/CellularComponent.Robj')

print("GO plots in OutputGO directory under GO_MF and GO_BP")


