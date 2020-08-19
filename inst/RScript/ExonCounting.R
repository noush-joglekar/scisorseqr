library(data.table)
library(dplyr)
library(tidyr)
library(parallel)
library(tictoc)

args <- commandArgs(trailingOnly=TRUE)
## arguments : AllInfo file, exonInfo.gff, threshold, grouping factor


allInfo <- fread(args[1])[,-c(7:8)]     ## Read in file with read - gene - celltype assignments
colnames(allInfo) <- c("Read","Gene","Celltype","Barcode","UMI","Stretch")
groupingFactor <- args[4]

geneInfo <- allInfo %>% select(Read,Gene,groupingFactor)

print("Reading in GFF file, will take a while")

exonGFF <- read.table(gzfile(args[2]))[,-c(2,3,6,8,9,11)]	## GFF file with exon info
colnames(exonGFF) <- c("chr","start","end","strand","readname")
exonGFF <- exonGFF %>% separate(readname,into=c("Read","path"),sep=".path") %>%
        select(-path) %>% unite(exon, c(chr,start,end,strand),sep="_",remove=FALSE)
exonGFF <- right_join(exonGFF,geneInfo) ## Merge gene information with GFF in case of multi-mapping exons

threshold <- as.integer(args[3])

tic()
uniqExons <- exonGFF %>% select(-c(Read,groupingFactor)) %>% distinct()		## Get unique exons to iterate over
edges <- exonGFF %>% group_by(Read) %>% mutate(s=min(start),e=max(end)) %>%
	select(Read,Gene,s,e,exon,groupingFactor) %>% distinct() %>% as.data.frame()
toc()
print(paste0("About to start processing exons per ",groupingFactor))

calcPSI <- function(x){		## x for line in uniqExons
	geneDF <- edges %>% filter(Gene == uniqExons$Gene[x]) %>% group_by(.dots = groupingFactor) %>%
			filter(s <= uniqExons$start[x] & e >= uniqExons$end[x]) %>% as.data.frame()
	numReads <- geneDF %>% group_by(.dots = groupingFactor, Read) %>% select(Read,groupingFactor) %>%
			ungroup() %>% group_by(.dots = groupingFactor) %>% distinct() %>%
			add_count() %>% filter(n>=threshold) %>% as.data.frame()
	geneDF <- geneDF %>% filter(Read %in% numReads$Read)  %>% as.data.frame()
	u_gf <- geneDF %>% select(groupingFactor) %>% distinct()	## make list of unique grouping factors
	psiGF <- NULL
	for(gf in u_gf[,1]){	## If total reads spanning exon more than sampling rate, extract those reads
		reads <- geneDF %>% filter(get(groupingFactor) == gf)
		psi = NULL
		sr <- reads %>% select(Read) %>% distinct()
		inc <- reads %>% filter(Read %in% sr$Read) %>%
			filter(exon == uniqExons$exon[x]) %>% nrow()
		psi <- inc/nrow(sr)	## Total reads spanning the exon = sr
		psi <- as.data.frame(psi)
		psi$inclusion <- inc
		psi$exclusion <- nrow(sr) - inc
		psi$exon <- uniqExons$exon[x]
		psi$Gene <- uniqExons$Gene[x]
		psi[,as.character(groupingFactor)] <- gf
		psiGF <- rbind(psiGF,psi) }
	return(psiGF)
}


tic()
#sampledPSI <- mclapply(1:nrow(uniqExons), function(x) calcPSI(x), mc.cores=28)
sampledPSI <- mclapply(1:200, function(x) calcPSI(x), mc.cores=28) ## testing purposes
toc()

print("Converting to dataframe")
psiDF <- bind_rows(sampledPSI) %>% select(exon,Gene,groupingFactor,inclusion,exclusion,psi)

write.table(psiDF, file=paste0("ExonInfo/InclusionExclusionCounts",groupingFactor,"_min",threshold,"Reads.csv"),
        sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
