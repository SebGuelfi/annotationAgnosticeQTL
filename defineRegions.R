### Manuel Sebastian Guelfi
### This script is implemented to use DERFINDER software
### to define the agnostic-regions
## Load libraries
## we select the chromosome
## print( chr <- as.numeric( commandArgs(trailingOnly=T)[1] ) )
chr <- "21"

library('derfinder')
library('GenomicRanges')
library('GenomicFeatures')

setwd("/home/sguelfi/projects/R/hipp/data/derfinder/fullCoverage/derfinder/")

load(paste0("fullCoverageChr",chr,".rda"))

load("~/projects/R/hipp/data/general/STAR.summary.rda")

librarySize <- as.numeric(as.character(unlist(lapply(summary.output.STAR,function(x){x[8,2]}))))
names(librarySize) <- names(summary.output.STAR)


# load("/SAN/neuroscience/WT_BRAINEAC/hipp/derfinder/scripts/STAR.summary.rda")
# listSamples <- names(summary.output.STAR)

txdb_ens87 <- makeTxDbFromGFF("/data/references/GTF/Homo_sapiens.GRCh38.87.gtf",format="gtf", dataSource="ENSEMBL", organism="Homo sapiens", taxonomyId="9606", circ_seqs=DEFAULT_CIRC_SEQS, chrominfo=NULL, miRBaseBuild="GRCh38")
txdb_ens87 <- keepSeqlevels(txdb_ens87, c(1:22, 'X', 'Y', 'MT'))

## the problem of the annotation have been described by Leo here:
## https://support.bioconductor.org/p/93235/
hg38_chrominfo <- fetchExtendedChromInfoFromUCSC("hg38")
new_info <- hg38_chrominfo$UCSC_seqlength[match(seqlevels(txdb_ens87), mapSeqlevels(hg38_chrominfo$UCSC_seqlevel, 'Ensembl'))]
names(new_info) <- names(seqlengths(txdb_ens87))
#seqlengths(txdb_ens87) <- new_info
rm(txdb_ens87)

library('rtracklayer')
gr <- import("/data/references/GTF/Homo_sapiens.GRCh38.87.gtf")
gr <- keepSeqlevels(gr, c(1:22, 'X', 'Y', 'MT'))
seqlengths(gr) <- new_info
genome(gr) <- 'hg38'
txdb_ens87 <- makeTxDbFromGRanges(gr)
rm(gr)

genomicState <- makeGenomicState(txdb_ens87, chrs = chr, style = 'Ensembl', currentStyle = 'Ensembl')

## getting the information about the library size
# libSize <- c()
# for (i in listSamples){
#   libSize[i] <- as.numeric(as.character(summary.output.STAR[[i]][grep("Uniquely mapped reads number |",summary.output.STAR[[i]][,"V1"],fixed=TRUE),"V2"]))
# }

## check that we select the rigth samples
stopifnot(identical(colnames(fullCov[[paste0("chr",chr)]]),names(librarySize)))
##dim(fullCovPUTM$chr21)
## we first filter
filteredCov <- lapply(fullCov, filterData, cutoff = 5, totalMapped=libSize, filter="mean",returnMean=TRUE)
## we delete the fullCov object so to empty the memory.
rm(fullCov)
## identify regions
expressedRegions <-regionMatrix(fullCov = filteredCov, maxRegionGap = 10L, maxClusterGap = 300L, L = 100, totalMapped=libSize)
## annotate regions
annotatedRegions <- annotateRegions(regions=expressedRegions[[paste0("chr",chr)]]$regions,genomicState=genomicState$fullGenome,fullOrCoding='full')

#### now we get the intergenic regions only
intergenicRegions <-  expressedRegions[[paste0("chr",chr)]]$regions[which(annotatedRegions$countTable$intergenic == 1)]
coverageIntergenic <- expressedRegions[[paste0("chr",chr)]]$coverageMatrix[which(annotatedRegions$countTable$intergenic == 1),]

## we save all the objects
save(intergenicRegions,coverageIntergenic,annotatedRegions,expressedRegions,filteredCov,file=paste0("/SAN/neuroscience/ukbec/analysis/hipp/derfinder/results/",chr,".rda"))

rm(fullCov,intergenicRegions,coverageIntergenic,annotatedRegions,expressedRegions,filteredCov)
