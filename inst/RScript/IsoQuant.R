## By Anoushka Joglekar 2019. Edited 03/2020

args <- commandArgs(trailingOnly=TRUE)

print("Reading in file, this might take a while")
all_info <- data.table::fread(args[1])[,c(1:4,6)]

is.iso <- args[2]
is.tss <- args[3]
is.polya <- args[4]

isoDir <- args[5]

colnames(all_info) <- c("Readname","Gene","Celltype","Barcode","Isoform")
if(is.tss != FALSE || is.polya != FALSE ) {
  data_df <- tidyr::separate(all_info, Isoform, into=c("TSS","Intron","PolyA"), sep="&&")}


if(is.iso != FALSE){
	ids1 <- all_info %>% group_by(Gene) %>% count(Isoform) %>% arrange(-n, .by_group=TRUE) %>% mutate(iso.id = row_number()-1)
	df1 <- inner_join(all_info,ids1, by = c("Gene","Isoform")) %>% select(-n, -Isoform)
	write.table(ids1,file=paste0(isoDir,"Iso-IsoID.csv"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
}

###

if(is.tss != FALSE){
	ids2 <- data_df %>% group_by(Gene) %>% count(TSS) %>% arrange(-n, .by_group=TRUE) %>% mutate(tss.id = row_number()-1)
	data_df <- inner_join(data_df,ids2, by = c("Gene","TSS")) %>% select(-n, -TSS)
	write.table(ids2,file=paste0(isoDir,"TSS-ID.csv"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
}

###

if(is.polya != FALSE){
	ids3 <- data_df %>% group_by(Gene) %>% count(PolyA) %>% arrange(-n, .by_group=TRUE) %>% mutate(polya.id = row_number()-1)
	data_df <- inner_join(data_df,ids3, by = c("Gene","PolyA")) %>% select(-n, -PolyA)
	write.table(ids3,file=paste0(isoDir,"PolyA-PolyID.csv"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
}

if(exists("data_df")){
tS <- colnames(data_df)[c(1, grep("id",colnames(data_df)))]
data_df <- data_df %>% select(tS)
outDF <- inner_join(df1,data_df)
} else {outDF <- df1}


write.table(outDF,file=paste0(isoDir,"Gene_Cluster_Cell_IsoIDs_PerRead.csv"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")



