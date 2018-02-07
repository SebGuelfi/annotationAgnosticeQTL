### Global variables

## first collapse all expression data

# path <- "/data/hipp/derfinder/results/" 
# gr <- NULL
# counts <- NULL
# ann.reg <- NULL
# 
# for (chr in c(1:22,"X")){
#     print(chr)
#     load(paste0(path,chr,".rda"))
#     
#     ## regions location
#     gr.tmp <- expressedRegions[[paste0("chr",chr)]]$regions[!width(expressedRegions[[paste0("chr",chr)]]$regions)<3]
#     gr.tmp$names <- c(paste0("ER",gsub("chr","",seqnames(gr.tmp)),"_",start(gr.tmp),"-",end(gr.tmp)))
#     names(gr.tmp) <- gr.tmp$names
#     
#     ## rename the "chr" from the chromosome name
#     gr.tmp <- renameSeqlevels(gr.tmp,unique(gsub("chr","",as.character(seqnames(gr.tmp)))))
# 
#     ## region counts
#     counts.tmp <- as.data.frame(expressedRegions[[paste0("chr",chr)]]$coverageMatrix[!width(expressedRegions[[paste0("chr",chr)]]$regions)<3,]) 
#     rownames(counts.tmp)<- names(gr.tmp)
#     
#     ##region annotation
#     ann.tmp <- annotatedRegions$countTable[!width(expressedRegions[[paste0("chr",chr)]]$regions)<3,]
#     rownames(ann.tmp) <- names(gr.tmp)
#     
#     
#     gr <- rbind(gr,as.data.frame(gr.tmp))
#     counts <- rbind(counts,counts.tmp)
#     ann.reg <- rbind(ann.reg,ann.tmp)
#     rm(gr.tmp,count.tmp,ann.tmp)
# }
# 
# gr <- GRanges(gr)
# stopifnot(identical(nrow(ann.reg),length(gr),nrow(counts)))
# save(gr,counts,ann.reg,file="/home/sguelfi/projects/R/hipp/data/expression/derfinder/mergedDerfinder.rda")



# load("/home/sguelfi/projects/R/hipp/data/expression/derfinder/mergedDerfinder.rda")
# fastaFile <- "/data/references/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# exprPath <- "/home/sguelfi/projects/R/hipp/data/expression/derfinder/RPKM.cqn/"
# 
# ## RPKMs ##################################################################
# library(easyRNASeq)
# message(Sys.time()," Calculates RPKM")
# 
# load("~/projects/R/hipp/data/general/STAR.summary.rda")
# 
# librarySize <- as.numeric(as.character(unlist(lapply(summary.output.STAR,function(x){x[8,2]}))))
# names(librarySize) <- names(summary.output.STAR)
# 
# regLength <- width(gr)
# names(regLength) <- gr$names
# 
# stopifnot(identical(names(librarySize),colnames(counts)))
# stopifnot(identical(names(regLength),rownames(counts)))
# 
# #convert in RPKM
# RPKM.std <- RPKM(as.matrix(counts), NULL, 
#                  lib.size=librarySize[colnames(counts)], 
#                  feature.size=regLength)
# 
# 
# 
# RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
# message(Sys.time()," regions filtered ",round((1-(nrow(RPKM.std)/nrow(counts)))*100,2),"%")
# 
# genesList <- rownames(RPKM.std)
# rm(RPKM.std,regLength)
# 
# counts <- counts[genesList,]
# gr <- gr[match(rownames(counts),gr$names)]
# 
# stopifnot(identical(gr$names,rownames(counts)))
# 
# ### GC content per region
# message("calculates the GC content per region")
# 
# regions.bed <- data.frame(seqnames=gsub("chr","",seqnames(gr)),
#                           starts=start(gr)-1,
#                           ends=end(gr),
#                           names=as.character(gr$names),
#                           scores=c(rep(".", length(gr))),
#                           strands=strand(gr))
# 
# 
# tmpf.bed <- tempfile("gr",fileext = ".bed")
# write.table(regions.bed, file=tmpf.bed, quote=F, sep="\t", row.names=F, col.names=F)
# 
# tmpf.gc <- tempfile()
# system(paste0("/tools/bedtools2/bin/bedtools nuc -fi ",fastaFile," -bed ",tmpf.bed," > ",tmpf.gc))
# GCcontentTab <- read.delim(tmpf.gc)
# rm(tmpf.bed,tmpf.gc)
# 
# library(cqn)
# library(scales)
# 
# stopifnot(identical(rownames(counts),as.character(GCcontentTab$X4_usercol)))
# message("Conditional quantile normalisation")
# 
# my.cqn <- cqn(counts=counts, lengths = GCcontentTab$X15_seq_len,x = GCcontentTab$X8_pct_gc, sizeFactors = librarySize , verbose = TRUE)
# 
# png(paste0("~/projects/R/hipp/plots/cqn/all.cqn.png"), type="cairo")
# par(mfrow=c(1,2))
# cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
# cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
# dev.off()
# RPKM.cqn <- my.cqn$y + my.cqn$offset
# 
# ## save the RPKM
# save(RPKM.cqn,file=paste0(exprPath,"RPKM.cqn.rda"))
# rm(ann.reg,counts,GCcontentTab,my.cqn,genesList,librarySize,regions.bed)


## have to remove some duplicared snps as well as those not present in the hg38
## system(paste0("grep -v rs4517039 ",pathImputed,"7.dosage.traw > ",pathImputed,"7.dosage.traw.tmp ; mv ",pathImputed,"7.dosage.traw.tmp ",pathImputed,"7.dosage.traw"))
## system(paste0("grep -v 144235357 ",pathImputed,"8.dosage.traw > ",pathImputed,"8.dosage.traw.tmp ; mv ",pathImputed,"8.dosage.traw.tmp ",pathImputed,"8.dosage.traw"))
## system(paste0("grep -v rs373336509 ",pathImputed,"10.dosage.traw > ",pathImputed,"10.dosage.traw.tmp ; mv ",pathImputed,"10.dosage.traw.tmp ",pathImputed,"10.dosage.traw"))
## system(paste0("grep -v rs200903336 ",pathImputed,"10.dosage.traw > ",pathImputed,"10.dosage.traw.tmp ; mv ",pathImputed,"10.dosage.traw.tmp ",pathImputed,"10.dosage.traw"))



library(GenomicRanges)
exprPath <- "/home/sguelfi/projects/R/hipp/data/expression/derfinder/RPKM.cqn/"
load(paste0(exprPath,"RPKM.cqn.rda"))
load("/home/sguelfi/projects/R/hipp/data/expression/derfinder/mergedDerfinder.rda")
rm(counts)
### data path
pathData <- "/home/sguelfi/projects/R/hipp/data/"
## pathImputed data
pathImputed <- "/home/sguelfi/projects/R/hipp/data/imputed_data/imputed_genotype_hg38/"


## add covariates 

load(file=paste0(pathData,"general/HIPPInfo.rda"))
rownames(hipp.info) <- paste0("Sample_",hipp.info$Sample_ID)
head(hipp.info)
genetic.pca <- read.delim("/home/sguelfi/projects/R/hipp/data/genotype/hipp.genetic.pca.eigenvec",sep =" ",header = F,
                          stringsAsFactors = F,colClasses = "character",as.is = T)
rownames(genetic.pca) <- paste(genetic.pca$V1,genetic.pca$V2,sep="_")
genetic.pca <- genetic.pca[,3:5]


load(paste0(pathData,"general/peer.axes.rda"))


nPeer <- 22
covs <- cbind(genetic.pca,as.numeric(as.factor(hipp.info[which(as.character(hipp.info$SD.No) %in% rownames(genetic.pca)),"Gender"]))
              ,hipp.info[which(as.character(hipp.info$SD.No) %in% rownames(genetic.pca)),"Age"],peer.axes[rownames(genetic.pca),1:(nPeer)])

expr.cqn <- RPKM.cqn

expr.cqn <- expr.cqn[,rownames(hipp.info)]

stopifnot(identical(colnames(expr.cqn),rownames(hipp.info)))
colnames(expr.cqn) <- hipp.info$SD.No

expr.gr <- gr[match(rownames(expr.cqn),names(gr))]
ann.reg <- ann.reg[match(rownames(expr.cqn),rownames(ann.reg)),]

source("/home/sguelfi/projects/R/annotationAgnosticeQTL/eQTLAnalysis_function.R")

library(doParallel)
library(foreach)

cl <- makeCluster(20)
clusterExport(cl, c("eQTL_analysis_derfinder","findOverlaps","seqnames","GRanges","IRanges"))
registerDoParallel(cl)
getDoParWorkers()

final_derfindereQTL <- foreach(i=1:23,.combine=rbind,.verbose=F)%dopar%eQTL_analysis_derfinder(i,expr.gr,expr.cqn,pathImputed,covs,outputFolder="/home/sguelfi/projects/R/hipp/data/results/derfinder/fullResults/")
save(final_derfindereQTL,file="/home/sguelfi/projects/R/hipp/data/results/derfinder/final_derfinder.rda")

stopCluster(cl)
rm(cl)


# seqnames(expr.gr)
# 
# load("/data1/users/seb/hipp/data/results/geneLeveleQTL.rda")
# length((unique(final_geneeQTL$gene)))
# (unique(final_geneeQTL$gene))
# table(final_geneeQTL$rank)
# length(unique(paste0(final_geneeQTL$gene,final_geneeQTL$rank)))
# 
# 
# rm(my.cov,my.expr,my.markers,outputFolder,results,j,chr,store,covs.tmp,dosage,dosage.gr,expr.tmp,dosage.tmp,expr.tmp.gr)
