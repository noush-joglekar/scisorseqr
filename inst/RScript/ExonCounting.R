#' @import dplyr
#' @importFrom magrittr %>%

`%>%` <- magrittr::`%>%`

args <- commandArgs(trailingOnly=TRUE)
## arguments : AllInfo file, exonInfo.gff, threshold, grouping factor


allInfo <- data.table::fread(args[1],sep="\t")[,-c(7:9)]     ## Read in file with read - gene - celltype assignments
colnames(allInfo) <- c("Read","Gene","Celltype","Barcode","UMI","Stretch")
groupingFactor <- args[4]

geneInfo <- allInfo %>% dplyr::select(Read,Gene,groupingFactor)

print("Reading in GFF file, will take a while")

exonGFF <- read.table(gzfile(args[2]))[,-c(2,3,6,8,9,11)]	## GFF file with exon info
colnames(exonGFF) <- c("chr","start","end","strand","readname")
exonGFF <- exonGFF %>% tidyr::separate(readname,into=c("Read","path"),sep=".path") %>%
  dplyr::select(-path) %>% tidyr::unite(exon, c(chr,start,end,strand),sep="_",remove=FALSE)
exonGFF <- dplyr::right_join(exonGFF,geneInfo) ## Merge gene information with GFF in case of multi-mapping exons

threshold <- as.integer(args[3])
numThreads <- as.integer(args[5])


uniqExons <- exonGFF %>% dplyr::select(-c(Read,groupingFactor)) %>% dplyr::distinct()		## Get unique exons to iterate over
edges <- exonGFF %>% dplyr::group_by(Read) %>% dplyr::mutate(s=min(start),e=max(end)) %>%
  dplyr::select(Read,Gene,s,e,exon,groupingFactor) %>% dplyr::distinct() %>% dplyr::as_data_frame()
print(paste0("About to start processing exons per ",groupingFactor))

calcPSI <- function(x){		## x for line in uniqExons
	geneDF <- edges %>% dplyr::filter(Gene == uniqExons$Gene[x]) %>% dplyr::group_by(.dots = groupingFactor) %>%
	  dplyr::filter(s <= uniqExons$start[x] & e >= uniqExons$end[x]) %>% as.data.frame()
	numReads <- geneDF %>% dplyr::group_by(.dots = groupingFactor, Read) %>% dplyr::select(Read,groupingFactor) %>%
	  dplyr::ungroup() %>% dplyr::group_by(.dots = groupingFactor) %>% dplyr::distinct() %>%
	  dplyr::add_count() %>% dplyr::filter(n>=threshold) %>% as.data.frame()
	geneDF <- geneDF %>% dplyr::filter(Read %in% numReads$Read)  %>% as.data.frame()
	u_gf <- geneDF %>% dplyr::select(groupingFactor) %>% dplyr::distinct()	## make list of unique grouping factors
	psiGF <- NULL
	for(gf in u_gf[,1]){	## If total reads spanning exon more than sampling rate, extract those reads
		reads <- geneDF %>% dplyr::filter(get(groupingFactor) == gf)
		psi = NULL
		sr <- reads %>% dplyr::select(Read) %>% dplyr::distinct()
		inc <- reads %>% dplyr::filter(Read %in% sr$Read) %>%
		  dplyr::filter(exon == uniqExons$exon[x]) %>% nrow()
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


sampledPSI <- parallel::mclapply(1:nrow(uniqExons), function(x) calcPSI(x), mc.cores=numThreads)
#sampledPSI <- parallel::mclapply(1:200, function(x) calcPSI(x), mc.cores=28) ## testing purposes

print("Converting to dataframe")
psiDF <- dplyr::bind_rows(sampledPSI) %>% dplyr::select(exon,Gene,groupingFactor,inclusion,exclusion,psi)

write.table(psiDF, file="ExonQuantOutput/InclusionExclusionCounts.tsv",
        sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
