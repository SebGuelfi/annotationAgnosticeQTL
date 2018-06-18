## Replication in the NABEC dataset


## here we will use the definition of the regions that we have identified in the 
## HIPP to then quantify the regions using the NABEC and then performed the 
## eQTL analysis

# load the libraries
library(DeepBlueR)
library(derfinder)
library(GenomicRanges)

pathImputed <- "/home/sguelfi/FCTX-NIH/imputed/"
pathDosage <- "/home/sguelfi/FCTX-NIH/imputed/dosage/"


## to do once
chr <- c(1:22,"X")
for(i in chr)
{
  print(i)
  system(paste0("vcftools --gzvcf ",pathImputed,"chr",i,".dose.vcf.gz --freq2 --out ",pathDosage,i))
  
  system(paste0("gunzip ",pathImputed,"chr",i,".info.gz"))
  tmp <- read.delim(paste0(pathImputed,"chr",i,".info"))
  head(tmp)
  tmp <- tmp[which(as.numeric(as.character(tmp$Rsq))>0.4 & tmp$MAF>=0.05),]
  
  tmp <- do.call(rbind,strsplit(as.character(tmp$SNP),":"))
  
  tmp <- tmp[,c(1,2,2)]
  colnames(tmp) <- c("CHR","start","end")
  write.table(tmp,file=paste0(pathImputed,i,".toFilter.BED"),row.names = F,quote = F)
  system(paste0("vcftools --gzvcf ",pathImputed,"chr",i,".dose.vcf.gz --recode --recode-INFO-all --exclude-bed ",pathImputed,i,".toFilter.BED --stdout | gzip -c > ",pathDosage,"chr",i,".tmp.vcf.gz"))
  system(paste0("mv ",pathDosage,"chr",i,".tmp.vcf.gz ",pathDosage,"chr",i,".vcf.gz"))
  #system(paste0("~/programs/plink/plink --vcf ",pathImputed,i,".vcf.gz --recode A-transpose --out ",pathImputed,i,".dosage --allow-extra-chr"))
  system(paste0("/tools/plink/plink --vcf ",pathDosage,"chr",i,".vcf.gz --recode A-transpose --out ",pathDosage,"chr",i,".dosage --allow-extra-chr"))
}

## load the regions to validate
load("/home/sguelfi/projects/R/hipp/data/expression/derfinder/mergedDerfinder.rda")
## to pass leo
## save(gr,file="/home/sguelfi/projects/R/hipp/data/expression/derfinder/regions.hg38.rda")


getCoverage_par <- function(i,L,gr){
    
    #gr.tmp <- regs[seqnames(regs)==paste0("chr",i)]
    ## Convert the regions to the hg19 
    gr.tmp.19 <- deepblue_liftover(gr, source = "hg38",target = "hg19")
    
    gr.tmp.19 <- gr.tmp.19[seqnames(gr.tmp.19)==paste0("chr",i)]
    ## get the list of samples analysed
    rm(gr.tmp)
    load(paste0("~/FCTX-NIH/derfinder/FCTXFullCoverage/fullCoverageChr",i,".rda"))
    #gr <- regs[seqnames(regs)==paste0("chr",i)]
    #gr <- gr[1:20]
    coverage <- getRegionCoverage(fullCov = fullCov, gr, totalMapped=NULL)
    covMat <- lapply(coverage, colSums)
    covMat <- do.call(rbind, covMat)
    covMat <- covMat / L
    rownames(covMat) <- paste0(seqnames(gr),"_",start(gr),"_",end(gr))
    
    save(counts,gr.tmp.19,file=paste0("/home/sguelfi/projects/R/hipp/data/validationNABEC/chr",i,".rda"))
    rm(fullCov,gr.tmp.19)
    return(TRUE)
}

no_cores <- 3 
library(foreach)
library(doParallel)

cl<-makeCluster(no_cores)
clusterExport(cl,c("deepblue_liftover", "rbind", "colSums", "lapply","do.call",
                   "seqnames","start","end"))
registerDoParallel(cl)

out <- foreach(i=c(1:22),.combine = c)  %dopar%  getCoverage_par(i,100,regs[seqnames(regs)==paste0("chr",i)])




# for(i in 1:22)
# {
#     ## load the fullCov
#     print(i)
#     load(paste0("~/FCTX-NIH/derfinder/FCTXFullCoverage/fullCoverageChr",i,".rda"))
#     
#     ## select only chromosome 21
#     gr.tmp <- regs[seqnames(regs)==paste0("chr",i)]
#     ## change the seqnames to have the same format in both datasets
#     #seqlevels(gr.tmp) <- paste0('chr',seqlevels(gr.tmp))
#     
#     ## Convert the regions to the hg19 
#     gr.tmp.19 <- deepblue_liftover(gr.tmp, source = "hg38",target = "hg19")
#     
#     gr.tmp.19 <- gr.tmp.19[seqnames(gr.tmp.19)==paste0("chr",i)]
#     
#     ## get the list of samples analysed
#     rm(gr.tmp)
#     #stopifnot(identical(names(libSize), gsub("fctx.fil.bam","",colnames(fullCov[[1]]))))
#     ## TRUE
#     
#     ## get the region coverage
#     coverage <- getRegionCoverage(fullCov = fullCov, gr.tmp.19, totalMapped=NULL)
#     covMat <- lapply(coverage, colSums)
#     covMat <- do.call(rbind, covMat)
#     covMat <- covMat / 10
#     rownames(covMat) <- paste0(seqnames(gr),"_",start(gr),"_",end(gr))
#     
#     
#     save(counts,gr.tmp.19,file=paste0("/home/sguelfi/projects/R/hipp/data/validationNABEC/chr",i,".rda"))
#     return(TRUE)
#     rm(fullCov,gr.tmp.19)
# }

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

load("/home/sguelfi/projects/R/hipp/data/expression/derfinder/RPKM.cqn/RPKM.cqn.rda")

ann.reg <- ann.reg[rownames(RPKM.cqn),]
tmp.notFilter <- ann.reg[intersect(rownames(ann.reg),rownames(RPKM.std)),]

tmp.filtered <- ann.reg[intersect(rownames(ann.reg),rownames(RPKM.std.filtered)),]

barplot((apply(tmp.filtered,2,sum)/apply(tmp.notFilter,2,sum))*100,ylim=c(0,100),ylab="percentage of validation", main="Hipp data - Validation NABEC")

load(file="~/projects/R/hipp/data/general/annotated.regions.rda")

intronic <- ann.reg[which(ann.reg$exon==0 & ann.reg$intron==1 & ann.reg$intergenic==0),]
intergenic <- ann.reg[which(ann.reg$exon==0 & ann.reg$intron==0 & ann.reg$intergenic==1),]

tmp.notFilter.intronic <- annotated.regions[which(as.character(annotated.regions$names) %in% rownames(intronic)),]

tmp.filtered <- annotated.regions[which(as.character(annotated.regions$names) %in% rownames(intronic)),]
tmp.filtered <- tmp.filtered[which(rownames(RPKM.std.filtered) %in% as.character(tmp.filtered$names)),]

table(is.na(tmp.notFilter.intronic$junID))["FALSE"]/nrow(tmp.notFilter.intronic)
## 2.5%

head(tmp.filtered)

table(is.na(tmp.filtered$junID))["FALSE"]/table(is.na(tmp.notFilter.intronic$junID))["FALSE"]
table(is.na(tmp.filtered$junID))["TRUE"]/table(is.na(tmp.notFilter.intronic$junID))["TRUE"]
## there is no difference between split annotated and split unannoted



## rm(list=ls())

## we then validate the eQTLs
fastaFile <- "/data/references/fasta/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
exprPath <- "/home/sguelfi/projects/R/hipp/data/validationNABEC/"

## RPKMs ##################################################################
library(easyRNASeq)
message(Sys.time()," Calculates RPKM")

load("~/FCTX-NIH/flagstat/flagstat.summary.rda")
listSamples <- read.table("~/FCTX-NIH/seb/listSamples.txt",header=F,colClasses="character")
libSize <- c()
for (i in listSamples$V1){
    libSize <- c(libSize,flagstatList[[i]][1,"reads"])
}
names(libSize) <- listSamples$V1
rm(i)

load(paste0("/home/sguelfi/projects/R/hipp/data/validationNABEC/chr.all.rda"))
colnames(counts.all) <- gsub("fctx.fil.bam","",colnames(counts.all))

load("/home/sguelfi/projects/R/hipp/data/results/final_derfinder.rda")

##load("/home/sguelfi/projects/R/hipp/data/expression/derfinder/mergedDerfinder.rda")

counts.all <- counts.all[as.character(intersect(rownames(counts.all),as.character(final_derfindereQTL$gene))),]
gr <- gr[as.character(intersect(rownames(counts.all),as.character(final_derfindereQTL$gene)))]

regLength <- width(gr)
names(regLength) <- gr$names

stopifnot(identical(names(libSize),colnames(counts.all)))
stopifnot(identical(names(regLength),rownames(counts.all)))

#convert in RPKM
RPKM.std <- RPKM(as.matrix(counts.all), NULL,
                 lib.size=libSize[colnames(counts.all)],
                 feature.size=regLength)

RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
#message(Sys.time()," regions filtered ",round((1-(nrow(RPKM.std)/nrow(counts)))*100,2),"%")

tmp.notFilter <- ann.reg[rownames(counts.all),]

tmp.filtered <- ann.reg[rownames(RPKM.std),]

barplot((apply(tmp.filtered,2,sum)/apply(tmp.notFilter,2,sum))*100,ylim=c(0,100),ylab="percentage of validation", main="Hipp eQTL data - Validation NABEC")

genesList <- rownames(RPKM.std)
rm(RPKM.std,regLength)

## select those genes that are expressed
counts <- counts.all[genesList,]

stopifnot(identical(gr$names,rownames(counts.all)))

### GC content per region
message("calculates the GC content per region")
library(GenomicFeatures)
library(GenomicRanges)
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

stopifnot(identical(rownames(counts.all),as.character(GCcontentTab$X4_usercol)))
message("Conditional quantile normalisation")
libSize <- libSize[colnames(counts.all)]

my.cqn <- cqn(counts=counts.all, lengths = GCcontentTab$X15_seq_len,x = GCcontentTab$X8_pct_gc, sizeFactors = libSize , verbose = TRUE)

head(my.cqn)
#png(paste0("~/projects/R/hipp/plots/cqn/all.cqn.png"), type="cairo")
#par(mfrow=c(1,2))
cqnplot(my.cqn, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(my.cqn, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
#dev.off()
RPKM.cqn <- my.cqn$y + my.cqn$offset

## save the RPKM
save(RPKM.cqn,file=paste0(exprPath,"RPKM.cqn.rda"))
rm(ann.reg,counts,GCcontentTab,my.cqn,genesList,libSize,regions.bed,tmp.notFilter,tmp.filtered,flagstatList,counts.all)

library(GenomicRanges)
exprPath <- "/home/sguelfi/projects/R/hipp/data/validationNABEC/"
load(paste0(exprPath,"RPKM.cqn.rda"))
#load("/home/sguelfi/projects/R/hipp/data/expression/derfinder/mergedDerfinder.rda")
### data path
pathData <- "/home/sguelfi/FCTX-NIH/"
## pathImputed data
pathImputed <- "/home/sguelfi/FCTX-NIH/imputed/dosage/"


## add covariates 
peer <- read.csv(paste0(pathData,"peer/X.csv"),header=F)
covars <- read.table(paste0(pathData,"sample_info/nabec.postQC.core.cohort.peer.covars.named.txt"),header=T,row.names=1)

## Check the covariates and get the appropriate headings
for(i in 1:19){
    print(paste(i,identical(as.numeric(covars[,i]),as.numeric(peer[,i]))))
}
# [1] "1 TRUE"
# [1] "2 TRUE"
# [1] "3 TRUE"
# [1] "4 TRUE"
# [1] "5 TRUE"
# [1] "6 TRUE"
# [1] "7 TRUE"
# [1] "8 TRUE"
# [1] "9 TRUE"
# [1] "10 TRUE"
# [1] "11 TRUE"
# [1] "12 TRUE"
# [1] "13 TRUE"
# [1] "14 TRUE"
# [1] "15 TRUE"
# [1] "16 TRUE"
# [1] "17 TRUE"
# [1] "18 TRUE"
# [1] "19 TRUE"


rownames(peer) <- rownames(covars)
# check the samples name match
identical(rownames(peer),
          gsub("_","-",
               gsub("KEN","UKY-",	
                    gsub("UM","UMARY-",
                         gsub("SH","SH-",colnames(RPKM.cqn))))))
colnames(RPKM.cqn) <- gsub("_","-",gsub("KEN","UKY-",gsub("UM","UMARY-",gsub("SH","SH-",colnames(RPKM.cqn)))))

identical(rownames(peer),colnames(RPKM.cqn))
## TRUE

## I live the first three pca axes and remove the rest

peer <- peer[,-c(11:20)]
expr.cqn <- RPKM.cqn

source("/home/sguelfi/projects/R/annotationAgnosticeQTL/eQTLAnalysis_function.R")

load("/home/sguelfi/projects/R/hipp/data/results/final_derfinder.rda")

## select only those regions that pass the expression cut-off, also we select those with the first rank.
toValidate <- final_derfindereQTL[match(intersect(as.character(final_derfindereQTL$gene),rownames(expr.cqn)),
                                                  as.character(final_derfindereQTL$gene)),]

## cannot find in the dosage all individuals so I will run this with those that are share
#expr.cqn <- expr.cqn[,intersect(colnames(expr.cqn),colnames(dosage))]
expr.gr <- expr.gr[match(rownames(expr.cqn),names(expr.gr))]

covs <- peer

## here I also remove the column for the UKY brain bank
# covs$V5 <- NULL
# covs$V6 <- NULL
# 
# covs$V4[covs$V4==0] <- 2


length(expr.gr)
nrow(expr.cqn)

library(doParallel)
library(foreach)


cl <- makeCluster(20)
clusterExport(cl, c("eQTL_analysis_derfinder_NABEC","findOverlaps","seqnames","GRanges","IRanges"))
registerDoParallel(cl)
getDoParWorkers()

final_derfindereQTL_NABEC <- foreach(i=1:22,.combine=rbind,.verbose=F)%dopar%eQTL_analysis_derfinder_NABEC(i,expr.gr,expr.cqn,pathImputed,covs,outputFolder="/home/sguelfi/projects/R/hipp/data/validationNABEC/results/")
save(final_derfindereQTL_NABEC,expr.gr,ann.reg,file="/home/sguelfi/projects/R/hipp/data/validationNABEC/results/final_derfinder.rda")

stopCluster(cl)
rm(cl)



# gr.tmp <- GRanges(seqnames = unlist(lapply(strsplit(as.character(final_derfindereQTL$snps),split = ":",fixed=T),function(x){return(x[1])})),
#         ranges = IRanges(as.numeric(unlist(lapply(strsplit(as.character(final_derfindereQTL$snps),split = ":",fixed=T),function(x){return(x[2])}))),
#                          as.numeric(unlist(lapply(strsplit(as.character(final_derfindereQTL$snps),split = ":",fixed=T),function(x){return(x[2])})))))
# 
# ## change the seqnames to have the same format in both datasets
# seqlevels(gr.tmp) <- paste0('chr',seqlevels(gr.tmp))
# 
# ## Convert the regions to the hg19 
# gr.tmp.19 <- deepblue_liftover(gr.tmp, source = "hg38",target = "hg19")
# 
# gr.tmp.19 <- gr.tmp.19[seqnames(gr.tmp.19)!=paste0("chr",i)]
# 
# ## we only test, 
# 
# head(final_derfindereQTL)
# 

load("~/FCTX-NIH/sample.totRNA.rda")


library(GenomicRanges)
library(DeepBlueR)

load("/home/sguelfi/projects/R/hipp/data/validationNABEC/results/final_derfinder.rda")
load("/home/sguelfi/projects/R/hipp/data/results/final_derfinder.rda")

toValidate <- final_derfindereQTL[-grep("X:",as.character(final_derfindereQTL$snps)),]
toValidate <- toValidate[-which(toValidate$rank >0),]


gr.tmp <- GRanges(seqnames = unlist(lapply(strsplit(as.character(toValidate$snps),split = ":",fixed=T),function(x){return(x[1])})),
                    ranges = IRanges(as.numeric(unlist(lapply(strsplit(as.character(toValidate$snps),split = ":",fixed=T),function(x){return(x[2])}))),
                    as.numeric(unlist(lapply(strsplit(as.character(toValidate$snps),split = ":",fixed=T),function(x){return(x[2])})))))

seqlevels(gr.tmp) <- paste0('chr',seqlevels(gr.tmp))

names(gr.tmp) <- as.character(toValidate$snps)

## Convert the regions to the hg19
gr.tmp.19 <- deepblue_liftover(gr.tmp, source = "hg38",target = "hg19")
gr.tmp.19 <- gr.tmp.19[!duplicated(names(gr.tmp.19))]
gr.tmp.19 <- as.data.frame(gr.tmp.19)

gr.tmp.19$hg19 <- paste0(gsub("chr","",gr.tmp.19$seqnames),":",gr.tmp.19$start)

toValidate$hg19 <- gr.tmp.19[as.character(toValidate$snps),"hg19"]

toValidate <- toValidate[!duplicated(as.character(toValidate$gene)),]
toValidate$validation <- "2"

for(i in 1:nrow(toValidate))
{
    #print(i)
    tryCatch({
            tmp <- read.delim(pipe(paste0("grep -w chr", as.character(toValidate$hg19[i]),
                                      " /home/sguelfi/projects/R/hipp/data/validationNABEC/results/chr",
                                      unlist(lapply(strsplit(as.character(toValidate$snps[i]),split = ":",fixed=T),function(x){return(x[1])})),"/",
                                      as.character(toValidate$gene[i]))),header = F)
        print(i)
        if(tmp[,6]<0.05)
        {
            toValidate$validation[i] <- "1"
        }else{
            toValidate$validation[i] <- "0"
        }
        
    }, error = function(e) {})

}


#save(toValidate,file="/home/sguelfi/projects/R/hipp/data/validationNABEC/toValidate.rda")



load(file="/home/sguelfi/projects/R/hipp/data/validationNABEC/toValidate.rda")




load("/home/sguelfi/projects/R/hipp/data/results/final_derfinder.rda")

library(tidyverse)
## first check how many regions in each category were used as input
tmp.all <- NULL
tmp.all["intron"] <- nrow(ann.reg %>% filter(intron>=1 & intergenic==0 & exon==0))
tmp.all["intergenic"] <- nrow(ann.reg %>% filter(intron==0 & intergenic>=1 & exon==0))
tmp.all["exonic"] <- nrow(ann.reg%>% filter(intron==0 & intergenic==0 & exon>=1))
tmp.all["exon-intron"] <- nrow(ann.reg %>% filter(intron>=1 & intergenic==0 & exon>=1))
tmp.all["exon-intergenic"] <- nrow(ann.reg %>% filter(intron==0 & intergenic>=1 & exon>=1))
tmp.all["exon-intron-intergenic"] <- nrow(ann.reg %>% filter(intron>=1 & intergenic>=1 & exon>=1))

barplot(sort(tmp.all,decreasing = T),ylim=c(0,300000),las=2,main="Regions tested")
barplot(sort(tmp.all,decreasing = T)/sum(tmp.all)*100,ylim=c(0,60),las=2,main="Regions tested", 
        ylab="Percentage of eQT")

tmp <- NULL
tmp["intron"] <- nrow(ann.reg[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron>=1 & intergenic==0 & exon==0))
tmp["intergenic"] <- nrow(ann.reg[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron==0 & intergenic>=1 & exon==0))
tmp["exonic"] <- nrow(ann.reg[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron==0 & intergenic==0 & exon>=1))
tmp["exon-intron"] <- nrow(ann.reg[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron>=1 & intergenic==0 & exon>=1))
tmp["exon-intergenic"] <- nrow(ann.reg[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron==0 & intergenic>=1 & exon>=1))
tmp["exon-intron-intergenic"] <- nrow(ann.reg[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron>=1 & intergenic>=1 & exon>=1))

length(unique(final_derfindereQTL$gene))

par(mar=c(11,4.1,3.1,2.1))
barplot(sort(tmp,decreasing = T),las=2,ylim=c(0,20000))

barplot(t(cbind(region_tested=sort(tmp.all,decreasing = T)/sum(tmp.all)*100,
                eQTL_target_region=sort((tmp/sum(tmp))*100,decreasing = T))),beside = T,las=2,ylim=c(0,60),
        legend.text = c("Tested regions","eQTL target regions"),ylab="Percentage of regions")




head(tmp[1:3])

ann.reg[as.character(toValidate$gene),]

validated <- toValidate[toValidate$validation=="1",]


nrow(ann.reg[unique(as.character(validated$gene)),] %>% filter(intron>=1 & intergenic==0 & exon==0))

nrow(ann.reg[unique(as.character(validated$gene)),] %>% filter(intron==0 & intergenic>=1 & exon==0))

nrow(ann.reg[unique(as.character(validated$gene)),] %>% filter(intron==0 & intergenic>=0 & exon>=1))

nrow(ann.reg[unique(as.character(toValidate$gene)),] %>% filter(intron>=1 & intergenic==0 & exon==0))

nrow(ann.reg[unique(as.character(toValidate$gene)),] %>% filter(intron==0 & intergenic>=1 & exon==0))

nrow(ann.reg[unique(as.character(toValidate$gene)),] %>% filter(intron==0 & intergenic>=0 & exon>=1))

rbind(c(intronic=15999,intergenic=2411,exonic=8147),
      c(intronic=305,intergenic=108,exonic=675))






