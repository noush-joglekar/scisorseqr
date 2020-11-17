#' @import dplyr
#' @import ggplot2
#' @import scales
#' @importFrom magrittr %>%

`%>%` <- magrittr::`%>%`

args <- commandArgs(trailingOnly=TRUE)
workingDir <- args[1]
print(paste0("Plotting for ",workingDir))

name = rev(unlist(strsplit(workingDir,"/")))[1]
print(name)

files = list.files(workingDir)

numsFile <- read.table(file.path(workingDir,files[grep('Nums',files)]),header = TRUE)
colnames(numsFile)[3:4] <- c("Group1","Group2")
resultsFile <- read.table(file.path(workingDir,files[grep('results',files)]),header=TRUE)

sig_ix <- resultsFile %>% dplyr::filter(FDR<=0.05) %>% dplyr::select(Gene)

########################################################################
######################## Pie chart section #############################
########################################################################

breakup <- c()
for (gene in sig_ix$Gene){
  cpos <- numsFile %>% dplyr::filter(Gene==gene) %>%
    dplyr::mutate(delta=Group1/sum(Group1)-Group2/sum(Group2)) %>%
    dplyr::arrange(-delta) %>% dplyr::filter(delta>0) %>%
    dplyr::mutate(cpos=cumsum(delta)) %>% dplyr::select(cpos)

  cneg <- numsFile %>% dplyr::filter(Gene==gene) %>%
    dplyr::mutate(delta=Group1/sum(Group1)-Group2/sum(Group2)) %>%
    dplyr::arrange(delta) %>% dplyr::filter(delta<=0) %>%
    dplyr::mutate(cneg=cumsum(abs(delta))) %>% dplyr::select(cneg)

  if(length(cpos$cpos) >= 2 && length(cneg$cneg) >= 2){
    if(cpos$cpos[2] >= cneg$cneg[2]){
      breakup <- c(breakup,which(cpos$cpos >= 0.1)[1])
    } else{breakup <- c(breakup,which(cneg$cneg >= 0.1)[1])}
  } else if(length(cpos$cpos) <= 2 || length(cneg$cneg) <= 2){
    breakup <- c(breakup,1)}
}

### Plotting section: Pie chart
total <- nrow(resultsFile)
nonSig <- length(which(resultsFile$FDR>0.05))

df <- as.data.frame(table(breakup))
df$breakup <- as.character(df$breakup)
df <- rbind(c("NonSig",nonSig),df)
df$value <- as.numeric(df$Freq)/total


df <- df %>% dplyr::mutate(ypos = cumsum(value)-0.5*value)
cols <- c('NonSig'="#F9C6BA",'1'="#930077",'2'="#DD6892",'3'="#3c6f9c",'4'="#29cdb5",'5'="#048998",'6'="#553c8b")

final_pie <- ggplot2::ggplot(df, aes(x="", y=value, fill=breakup)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = percent(value)), color = "white", size=6) +
  scale_fill_manual("",values=cols)


pdf(file=paste0("Visualizations/PieCharts_BreakUp_",name,".pdf"),8,6)
print(final_pie)
dev.off()
