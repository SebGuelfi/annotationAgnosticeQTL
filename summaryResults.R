load("/home/sguelfi/projects/R/hipp/data/results/final_derfinder.rda")

library(tidyverse)
## first check how many regions in each category were used as input
tmp.all <- NULL
tmp.all["intron"] <- nrow(regs.ann$countTable %>% filter(intron>=1 & intergenic==0 & exon==0))
tmp.all["intergenic"] <- nrow(regs.ann$countTable %>% filter(intron==0 & intergenic>=1 & exon==0))
tmp.all["exonic"] <- nrow(regs.ann$countTable%>% filter(intron==0 & intergenic==0 & exon>=1))
tmp.all["exon-intron"] <- nrow(regs.ann$countTable %>% filter(intron>=1 & intergenic==0 & exon>=1))
tmp.all["exon-intergenic"] <- nrow(regs.ann$countTable %>% filter(intron==0 & intergenic>=1 & exon>=1))
tmp.all["exon-intron-intergenic"] <- nrow(regs.ann$countTable %>% filter(intron>=1 & intergenic>=1 & exon>=1))

barplot(sort(tmp.all,decreasing = T),ylim=c(0,300000),las=2,main="Regions tested")
barplot(sort(tmp.all,decreasing = T)/sum(tmp.all)*100,ylim=c(0,60),las=2,main="Regions tested", 
        ylab="Percentage of eQT")

tmp <- NULL
tmp["intron"] <- nrow(regs.ann$countTable[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron>=1 & intergenic==0 & exon==0))
tmp["intergenic"] <- nrow(regs.ann$countTable[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron==0 & intergenic>=1 & exon==0))
tmp["exonic"] <- nrow(regs.ann$countTable[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron==0 & intergenic==0 & exon>=1))
tmp["exon-intron"] <- nrow(regs.ann$countTable[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron>=1 & intergenic==0 & exon>=1))
tmp["exon-intergenic"] <- nrow(regs.ann$countTable[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron==0 & intergenic>=1 & exon>=1))
tmp["exon-intron-intergenic"] <- nrow(regs.ann$countTable[unique(as.character(final_derfindereQTL$gene)),] %>% filter(intron>=1 & intergenic>=1 & exon>=1))

length(unique(final_derfindereQTL$gene))

par(mar=c(11,4.1,3.1,2.1))
barplot(sort(tmp,decreasing = T),las=2,ylim=c(0,20000))

barplot(t(cbind(region_tested=sort(tmp.all,decreasing = T)/sum(tmp.all)*100,
    eQTL_target_region=sort((tmp/sum(tmp))*100,decreasing = T))),beside = T,las=2,ylim=c(0,60),
    legend.text = c("Tested regions","eQTL target regions"),ylab="Percentage of regions")
      

barplot(sort(tmp,decreasing = T),las=2,ylim=c(0,20000))

## proportions
barplot(sort((tmp/sum(tmp))*100,decreasing = T),las=2,ylim=c(0,60), main="eQTL target regions")

##(17389+1309)/8909
#2.1%
### Annotate with junctions

setwd("/home/sguelfi/projects/R/annotatER/")
library(devtools)
load_all()

## set the path and take the names of the samples
## to create the proper path in ubuntu I have to do mount the neuroscience path like this.
## sshfs skgtmgu@wise.cs.ucl.ac.uk:/SAN/neuroscience /home/sguelfi/neuroscience/

path <- "/data/hipp/2ndPassSTAR/"

list.samples = list.dirs(path = path,recursive = F)
names(list.samples) <- basename(list.samples)
list.samples[1:length(list.samples)] <- paste0(list.samples,"/align_star.bamSJ.out.tab")
## list.samples = the list of samples to collect
## minSamples = filter split reads that do not have this min of number of samples

tmp.table <- NULL
minSamples = 5  ## in percentage

STARSplitRead <- loadSplitReads(list.samples = list.samples,minSamples = 5)

## Separate counts and create the junction ID so that can be mapped back
junctions <- list()
junctions[["annotation"]] <- cbind(junID=1:nrow(STARSplitRead),STARSplitRead[,c(1:4,7)])
junctions[["counts"]] <- cbind(junID=1:nrow(STARSplitRead),STARSplitRead[,8:ncol(STARSplitRead)])

#load("~/projects/R/hipp/data/expression/splitReads/splitReads.filtered.rda")

GTFPath <- "/data/references/GTF/Homo_sapiens.GRCh38.87.gtf"

## annotate the junctions

system.time(splitReadTable <- annotateSplitReads(GTFPath,junctions[["annotation"]]))

library(SummarizedExperiment)
regions <- expr.gr[match(unique(as.character(final_derfindereQTL$gene)),names(expr.gr))]

annotated.regions <- annotatERJunction(regions,splitReadTable)



tmp <- NULL
ann.reg.tmp <- ann.reg[unique(as.character(annotated.regions[!is.na(annotated.regions$junID),"names"])),]
tmp["intron"] <- nrow(ann.reg.tmp %>% filter(intron>=1 & intergenic==0 & exon==0))
tmp["intergenic"] <- nrow(ann.reg.tmp %>% filter(intron==0 & intergenic>=1 & exon==0))
tmp["exonic"] <- nrow(ann.reg.tmp %>% filter(intron==0 & intergenic==0 & exon>=1))
tmp["exon-intron"] <- nrow(ann.reg.tmp %>% filter(intron>=1 & intergenic==0 & exon>=1))
tmp["exon-intergenic"] <- nrow(ann.reg.tmp %>% filter(intron==0 & intergenic>=1 & exon>=1))
tmp["exon-intron-intergenic"] <- nrow(ann.reg.tmp %>% filter(intron>=1 & intergenic>=1 & exon>=1))

barplot(tmp)

tmp["intron"]
tmp["intergenic"]
578+86

578/16324

    
