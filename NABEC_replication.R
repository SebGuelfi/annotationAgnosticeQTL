## Replication in the NABEC dataset


## here we will use the definition of the regions that we have identified in the 
## HIPP to then quantify the regions using the NABEC and then performed the 
## eQTL analysis

# load the libraries
library(DeepBlueR)
library(derfinder)
library(GenomicRanges)

## load the regions to validate
load("/home/sguelfi/projects/R/hipp/data/expression/derfinder/mergedDerfinder.rda")
## to pass leo
## save(gr,file="/home/sguelfi/projects/R/hipp/data/expression/derfinder/regions.hg38.rda")
## load the library size
load("~/FCTX-NIH/flagstat/flagstat.summary.rda")

listSamples <- read.table("~/FCTX-NIH/seb/listSamples.txt",header=F,colClasses="character")
libSize <- c()
for (i in listSamples$V1){
    libSize <- c(libSize,flagstatList[[i]][1,"reads"])
}
names(libSize) <- listSamples$V1
rm(listSamples,i)


for(i in 11:22)
{
    ## load the fullCov
    print(i)
    load(paste0("~/FCTX-NIH/derfinder/FCTXFullCoverage/fullCoverageChr",i,".rda"))
    
    ## select only chromosome 21
    gr.tmp <- gr[seqnames(gr)==i]
    ## change the seqnames to have the same format in both datasets
    seqlevels(gr.tmp) <- paste0('chr',seqlevels(gr.tmp))
    
    ## Convert the regions to the hg19 
    gr.tmp.19 <- deepblue_liftover(gr.tmp, source = "hg38",target = "hg19")
    
    gr.tmp.19 <- gr.tmp.19[seqnames(gr.tmp.19)==paste0("chr",i)]
    
    ## get the list of samples analysed
    rm(gr.tmp)
    stopifnot(identical(names(libSize), gsub("fctx.fil.bam","",colnames(fullCov[[1]]))))
    ## TRUE
    
    ## get the region coverage
    countsNABEC <- getRegionCoverage(fullCov = fullCov, gr.tmp.19, totalMapped =libSize)
    ## get the mean
    counts <- lapply(countsNABEC, function(x){
        apply(x,2,mean)
    })
    
    names(counts) <- names(gr.tmp.19)
    counts <- do.call(rbind,counts)
    save(counts,gr.tmp.19,file=paste0("/home/sguelfi/projects/R/hipp/data/validationNABEC/chr",i,".rda"))
    rm(fullCov,gr.tmp.19)
}

## merge all results
## when converting from 38 to 19 some of the regions get splitted and those get quantified
## However, for replication purposes, we will remove those regions.

counts.all <- NULL
gr <- NULL

for(i in 1:22)
{
    print(i)
    load(paste0("/home/sguelfi/projects/R/hipp/data/validationNABEC/chr",i,".rda"))
    if (i ==1)
    {
        toRemove <- rownames(counts)[duplicated(rownames(counts))]
        counts <- counts[-which(rownames(counts) %in% toRemove),]
        gr.tmp.19 <- gr.tmp.19[-which(names(gr.tmp.19) %in% toRemove)]
        counts.all <- counts
        gr <- gr.tmp.19
    }else{
        stopifnot(identical(colnames(counts.all),colnames(counts)))
        toRemove <- rownames(counts)[duplicated(rownames(counts))]
        counts <- counts[-which(rownames(counts) %in% toRemove),]
        gr.tmp.19 <- gr.tmp.19[-which(names(gr.tmp.19) %in% toRemove)]
        counts.all <- rbind(counts.all,counts)
        gr <- makeGRangesFromDataFrame(rbind(as.data.frame(gr),as.data.frame(gr.tmp.19)),keep.extra.columns=T)
    }
}

save(counts.all,gr,file=paste0("/home/sguelfi/projects/R/hipp/data/validationNABEC/chr.all.rda"))

## we validate the expression

library(easyRNASeq)

load(paste0("/home/sguelfi/projects/R/hipp/data/validationNABEC/chr.all.rda"))

colnames(counts.all) <- gsub("fctx.fil.bam","",colnames(counts.all))

load("~/FCTX-NIH/flagstat/flagstat.summary.rda")
listSamples <- read.table("~/FCTX-NIH/seb/listSamples.txt",header=F,colClasses="character")
libSize <- c()
for (i in listSamples$V1){
    libSize <- c(libSize,flagstatList[[i]][1,"reads"])
}
names(libSize) <- listSamples$V1

regLength <- width(gr)
names(regLength) <- gr$names

stopifnot(identical(names(libSize),colnames(counts.all)))
stopifnot(identical(names(regLength),rownames(counts.all)))

RPKM.std <- RPKM(as.matrix(counts.all), NULL,
                 lib.size=libSize[colnames(counts.all)],
                 feature.size=regLength)

RPKM.std.filtered <- RPKM.std[rowSums(RPKM.std>0.1)>=(0.8*ncol(RPKM.std)),]

nrow(RPKM.std.filtered)/nrow(RPKM.std)

load("/home/sguelfi/projects/R/hipp/data/expression/derfinder/mergedDerfinder.rda")

tmp.notFilter <- ann.reg[rownames(RPKM.std),]

tmp.filtered <- ann.reg[rownames(RPKM.std.filtered),]

barplot((apply(tmp.filtered,2,sum)/apply(tmp.notFilter,2,sum))*100,ylim=c(0,100),ylab="percentage of validation", main="Hipp data - Validation NABEC")






## we then validate the eQTLs
fastaFile <- "/data/references/fasta/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
exprPath <- "/home/sguelfi/projects/R/hipp/data/validationNABEC/"

## RPKMs ##################################################################
library(easyRNASeq)
message(Sys.time()," Calculates RPKM")

listSamples <- read.table("~/FCTX-NIH/seb/listSamples.txt",header=F,colClasses="character")
libSize <- c()
for (i in listSamples$V1){
    libSize <- c(libSize,flagstatList[[i]][1,"reads"])
}
names(libSize) <- listSamples$V1

regLength <- width(gr)
names(regLength) <- gr$names

stopifnot(identical(names(librarySize),colnames(counts)))
stopifnot(identical(names(regLength),rownames(counts)))

#convert in RPKM
RPKM.std <- RPKM(as.matrix(counts), NULL,
                 lib.size=librarySize[colnames(counts)],
                 feature.size=regLength)

RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
message(Sys.time()," regions filtered ",round((1-(nrow(RPKM.std)/nrow(counts)))*100,2),"%")

genesList <- rownames(RPKM.std)
rm(RPKM.std,regLength)

counts <- counts[genesList,-match(c("Sample_A653_817","Sample_A653_856","Sample_A653_813",
                                    "Sample_A653_1073","Sample_A653_1001","Sample_A653_1278"),
                                  colnames(counts))]
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
librarySize <- librarySize[colnames(counts)]

my.cqn <- cqn(counts=counts, lengths = GCcontentTab$X15_seq_len,x = GCcontentTab$X8_pct_gc, sizeFactors = librarySize , verbose = TRUE)

png(paste0("~/projects/R/hipp/plots/cqn/all.cqn.png"), type="cairo")
par(mfrow=c(1,2))
cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
dev.off()
RPKM.cqn <- my.cqn$y + my.cqn$offset

## save the RPKM
save(RPKM.cqn,file=paste0(exprPath,"RPKM.cqn.rda"))
rm(ann.reg,counts,GCcontentTab,my.cqn,genesList,librarySize,regions.bed)












