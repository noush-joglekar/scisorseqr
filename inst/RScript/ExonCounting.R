#' @importFrom magrittr %>%
#' @import parallel
#' @import dplyr
#' @export

`%>%` <- magrittr::`%>%`

args <- commandArgs(trailingOnly=TRUE)
## arguments : AllInfo file, groupingFactor, nThreads, threshold

allInfo <- data.table::fread(args[1],sep="\t")
if(ncol(allInfo) == 11){
        allInfo <- allInfo[,c(1:4,9)]
} else { allInfo <- allInfo[,c(1:4,7)]}  ## Read in file with read - gene - celltype assignments

colnames(allInfo) <- c("Read","Gene","Celltype","Barcode","Exons")
groupingFactor <- args[2]

nThreads <- as.integer(args[3])
rpg <- as.integer(args[4])

outDir <- 'ExonQuantOutput/'

allInfo_SE <- allInfo %>% dplyr::group_by(Gene) %>% dplyr::add_count() %>%
        dplyr::filter(n >= rpg) %>% dplyr::select(-n) %>% dplyr::rowwise() %>%
        dplyr::mutate(start = unlist(strsplit(Exons,"_"))[2],
                end = rev(unlist(strsplit(Exons,"_")))[2]) %>%
        as.data.frame()

internalExons = allInfo_SE %>% tidyr::separate_rows(Exons,sep = ";%;") %>% dplyr::filter(Exons != "") %>%
        dplyr::group_by(Read) %>% dplyr::add_count() %>% dplyr::filter(n>=3) %>% dplyr::slice(2:(dplyr::n()-1))

uniqExons <- internalExons %>% dplyr::ungroup() %>% dplyr::select(Exons, Gene) %>% dplyr::distinct()

inclusionCounts <- internalExons %>% dplyr::ungroup() %>% dplyr::select(Exons,Gene,all_of(groupingFactor)) %>%
        dplyr::group_by(Exons,Gene,.dots=groupingFactor) %>% dplyr::add_count(name = "inclusion") %>%
        dplyr::distinct() %>% as.data.frame()

readSE <- internalExons %>% dplyr::select(Read,Gene,all_of(groupingFactor),start,end) %>%
	dplyr::distinct() %>% as.data.frame()

checkSpanningReads <- function(gene){
	tid <- as.character(Sys.getpid())
        exons <- uniqExons %>% dplyr::filter(Gene == gene) %>%
                tidyr::separate(Exons,into=c("chr","s","e","strand"),sep = "_",remove=FALSE)
        reads <- readSE %>% dplyr::filter(Gene == gene)
        spanningReads <- dplyr::left_join(reads,exons,by = "Gene") %>%
		mutate_at(c('s','e','start','end'), as.integer) %>%
		dplyr::filter(s >= start & e <= end) %>%
                dplyr::select(Exons,Gene,all_of(groupingFactor)) %>% dplyr::group_by(Exons,Gene,.dots = groupingFactor) %>%
                dplyr::add_count(name = "Total") %>% dplyr::distinct() %>% as.data.frame()
        inclusionReads <- inclusionCounts %>% dplyr::filter(Gene == gene)
        inc_tot <- dplyr::right_join(inclusionReads,spanningReads,by = c('Exons',"Gene",groupingFactor)) %>%
                replace(is.na(.),0) %>%
                dplyr::mutate(PSI = inclusion/Total, exclusion = Total-inclusion) %>% as.data.frame()
	inc_tot <- inc_tot[,c(1:4,7,5,6)]
        fName <- file.path(outDir,paste0("PID_",tid,"_InclusionExclusionCounts.tsv"))
        if(file.exists(fName)){
                write.table(inc_tot,fName,sep ="\t",
                        append= T, quote = F, row.names = F, col.names = FALSE)
        } else{
               	write.table(inc_tot,fName,sep ="\t",
                        quote = F, row.names = F, col.names = TRUE)
        }
	return
}

byGene = parallel::mclapply(unique(uniqExons$Gene), function(g) checkSpanningReads(g), mc.cores = nThreads)
