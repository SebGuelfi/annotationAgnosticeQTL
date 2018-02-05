### Global variables
chr <- 21 
path <- "/data/hipp/derfinder/results/" 
fastaFile <- "/data/references/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
exprPath <- "/home/sguelfi/projects/R/hipp/data/expression/derfinder/RPKM.cqn/"



load(paste0(path,chr,".rda"))
rm(filteredCov,coverageIntergenic)


## filter those regions shorter than 3bp
gr <- expressedRegions[[paste0("chr",chr)]]$regions[!width(expressedRegions[[paste0("chr",chr)]]$regions)<3]
counts <- expressedRegions[[paste0("chr",chr)]]$coverageMatrix[!width(expressedRegions[[paste0("chr",chr)]]$regions)<3,] 

## assign indexes
rownames(counts)<- c(paste0("ER",gsub("chr","",seqnames(gr)),":",start(gr),"-",end(gr)))
gr$names <- c(paste0("ER",gsub("chr","",seqnames(gr)),":",start(gr),"-",end(gr)))


## get the GC content 

## RPKMs ##################################################################
library(easyRNASeq)
message(Sys.time()," Calculates RPKM")

load("~/projects/R/hipp/data/general/STAR.summary.rda")

librarySize <- as.numeric(as.character(unlist(lapply(summary.output.STAR,function(x){x[8,2]}))))
names(librarySize) <- names(summary.output.STAR)

regLength <- width(gr)
names(regLength) <- gr$names

stopifnot(identical(names(librarySize),colnames(counts)))
stopifnot(identical(names(regLength),rownames(counts)))

#convert in RPKM
RPKM.std <- RPKM(counts, NULL, 
                 lib.size=librarySize[colnames(counts)], 
                 feature.size=regLength)



RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
message(Sys.time()," regions filtered ",round((1-(nrow(RPKM.std)/nrow(counts)))*100,2),"%")

genesList <- rownames(RPKM.std)
rm(RPKM.std,regLength)

counts <- counts[genesList,]
gr <- gr[match(rownames(counts),gr$names)]

stopifnot(identical(gr$names,rownames(counts)))

### GC content per region
message("calculates the GC content per region")

regions.bed <- data.frame(seqnames=gsub("chr","",seqnames(gr)),
                          starts=start(gr)-1,
                          ends=end(gr),
                          names=as.character(gr$names),
                          scores=c(rep(".", length(gr))),
                          strands=strand(gr))


tmpf.bed <- tempfile("gr",fileext = ".bed")
write.table(regions.bed, file=tmpf.bed, quote=F, sep="\t", row.names=F, col.names=F)

tmpf.gc <- tempfile()
system(paste0("/tools/bedtools2/bin/bedtools nuc -fi ",fastaFile," -bed ",tmpf.bed," > ",tmpf.gc))
GCcontentTab <- read.delim(tmpf.gc)
rm(tmpf.bed,tmpf.gc)

library(cqn)
library(scales)

stopifnot(identical(rownames(counts),as.character(GCcontentTab$X4_usercol)))
message("Conditional quantile normalisation")

my.cqn <- cqn(counts=counts, lengths = GCcontentTab$X15_seq_len,x = GCcontentTab$X8_pct_gc, sizeFactors = librarySize , verbose = TRUE)

png(paste0("~/projects/R/hipp/plots/cqn/chr",chr,"cqn.png"), type="cairo")
par(mfrow=c(1,2))
cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
dev.off()
RPKM.cqn <- my.cqn$y + my.cqn$offset

save(RPKM.cqn,file=paste0(exprPath,chr,".rda"))

expr.gr <- gr

## server path
## pathData <- "/home/seb/projectsR/hipp/data/"
# pathImputed data
## pathImputed <- "/home/seb/hipp/imputed_genotype/imputed_genotype_hg38/"
### local path
pathData <- "/home/sguelfi/projects/R/hipp/data/"
## pathImputed data
pathImputed <- "/home/sguelfi/projects/R/hipp/data/imputed_data/imputed_genotype_hg38/"



library(doParallel)
library(foreach)

cl <- makeCluster(20)
clusterExport(cl, c("eQTL_analysis","findOverlaps","seqnames","GRanges","IRanges"))
registerDoParallel(cl)
getDoParWorkers()

final_geneeQTL <- foreach(i=1:23,.combine=rbind,.verbose=F)%dopar%eQTL_analysis(i,expr.gr,expr.qn,pathImputed,covs,outputFolder="/data1/users/seb/hipp/data/results/geneLevel/fullResults/")
save(final_geneeQTL,file="/data1/users/seb/hipp/data/results/geneLeveleQTL.rda")

stopCluster(cl)
rm(cl)


load("/data1/users/seb/hipp/data/results/geneLeveleQTL.rda")
length((unique(final_geneeQTL$gene)))
(unique(final_geneeQTL$gene))
table(final_geneeQTL$rank)
length(unique(paste0(final_geneeQTL$gene,final_geneeQTL$rank)))


