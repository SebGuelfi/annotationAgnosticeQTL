### Global variables

## first collapse all expression data

path <- "/data/hipp/derfinder/results/"
gr <- NULL
counts <- NULL
ann.reg <- NULL

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

# gr <- GRanges(gr)
# stopifnot(identical(nrow(ann.reg),length(gr),nrow(counts)))
# save(gr,counts,ann.reg,file="/home/sguelfi/projects/R/hipp/data/expression/derfinder/mergedDerfinder.rda")

load("/home/sguelfi/projects/R/hipp/data/expression/derfinder/mergedDerfinder.rda")
nrow(counts)
## 385121 regions analysed

fastaFile <- "/data/references/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
exprPath <- "/home/sguelfi/projects/R/hipp/data/expression/derfinder/RPKM.cqn/"

## RPKMs ##################################################################
library(easyRNASeq)
message(Sys.time()," Calculates RPKM")

load("~/projects/R/hipp/data/general/STAR.summary.rda")

librarySize <- as.numeric(as.character(unlist(lapply(summary.output.STAR,function(x){x[8,2]}))))
names(librarySize) <- names(summary.output.STAR)

regLength <- width(regs)
names(regLength) <- names(regs)

stopifnot(identical(names(librarySize),colnames(counts)))
counts <- counts[names(regs),]
stopifnot(identical(names(regLength),rownames(counts)))

#convert in RPKM
RPKM.std <- RPKM(as.matrix(counts), NULL,
                 lib.size=librarySize[colnames(counts)],
                 feature.size=regLength)

RPKM.std <- RPKM.std[,-match(c("Sample_A653_817","Sample_A653_856","Sample_A653_813",
                               "Sample_A653_1073","Sample_A653_1001","Sample_A653_1278"),
                             colnames(RPKM.std))]

RPKM.std=RPKM.std[rowSums(RPKM.std>=0.1)>(ncol(RPKM.std)-((ncol(RPKM.std)*20)/100)),]
message(Sys.time()," regions filtered ",round((1-(nrow(RPKM.std)/nrow(counts)))*100,2),"%")

genesList <- rownames(RPKM.std)
rm(RPKM.std,regLength)

counts <- counts[genesList,-match(c("Sample_A653_817","Sample_A653_856","Sample_A653_813",
                                    "Sample_A653_1073","Sample_A653_1001","Sample_A653_1278"),
                                  colnames(counts))]
regs <- regs[match(rownames(counts),names(regs))]

stopifnot(identical(names(regs),rownames(counts)))

### GC content per region
message("calculates the GC content per region")

regions.bed <- data.frame(seqnames=gsub("chr","",seqnames(regs)),
                          starts=start(regs)-1,
                          ends=end(regs),
                          names=as.character(names(regs)),
                          scores=c(rep(".", length(regs))),
                          strands=strand(regs))


tmpf.bed <- tempfile("regs",fileext = ".bed")
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
rm(regs.ann,regs,counts,GCcontentTab,my.cqn,genesList,librarySize,
   regions.bed,RPKM.cqn,summary.output.STAR,exprPath,fastaFile)


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
pathImputed <- "/home/sguelfi/projects/R/hipp/data/imputed_genotype_v2/imputed_genotype_hg38_95/"


## add covariates 

load(file=paste0(pathData,"general/HIPPInfo.rda"))
rownames(hipp.info) <- paste0("Sample_",hipp.info$Sample_ID)
head(hipp.info)
genetic.pca <- read.delim("/home/sguelfi/projects/R/hipp/data/genotype/hipp.genetic.pca.eigenvec",sep =" ",header = F,
                          stringsAsFactors = F,colClasses = "character",as.is = T)
rownames(genetic.pca) <- paste(genetic.pca$V1,genetic.pca$V2,sep="_")
genetic.pca <- genetic.pca[as.character(hipp.info$SD.No),]
genetic.pca <- genetic.pca[,3:5]


load(paste0(pathData,"general/peer.axes.rda"))


nPeer <- 19
covs <- cbind(genetic.pca,as.numeric(as.factor(hipp.info[which(as.character(hipp.info$SD.No) %in% rownames(genetic.pca)),"Gender"]))
              ,hipp.info[which(as.character(hipp.info$SD.No) %in% rownames(genetic.pca)),"Age"],
              as.numeric(as.factor(hipp.info[which(as.character(hipp.info$SD.No) %in% rownames(genetic.pca)),"Brain.Bank"])),
              peer.axes[rownames(genetic.pca),1:(nPeer)])

expr.cqn <- RPKM.cqn

expr.cqn <- expr.cqn[,rownames(hipp.info)]

stopifnot(identical(colnames(expr.cqn),rownames(hipp.info)))
colnames(expr.cqn) <- hipp.info$SD.No

expr.gr <- regs[match(rownames(expr.cqn),names(regs))]

regs.ann <- regs.ann[match(rownames(expr.cqn),rownames(regs.ann))]

source("/home/sguelfi/projects/R/annotationAgnosticeQTL/eQTLAnalysis_function.R")

library(doParallel)
library(foreach)

## run in parallel

cl <- makeCluster(20)
clusterExport(cl, c("eQTL_analysis_derfinder","findOverlaps","seqnames","GRanges","IRanges"))
registerDoParallel(cl)
getDoParWorkers()

final_derfindereQTL <- foreach(i=1:23,.combine=rbind,.verbose=F)%dopar%eQTL_analysis_derfinder(i,expr.gr,expr.cqn,pathImputed,covs,outputFolder="/home/sguelfi/projects/R/hipp/data/results/derfinder/fullResults/")
save(final_derfindereQTL,expr.gr,ann.reg,file="/home/sguelfi/projects/R/hipp/data/results/final_derfinder.rda")

stopCluster(cl)
rm(cl)


### collect unsentinalised

load("/home/sguelfi/projects/R/hipp/data/results/final_derfinder.rda")
load(file="/home/sguelfi/projects/R/hipp/data/general/snp.map.rda")

## we first collect the transcript eQTL
#head(final_derfindereQTL)

regionID <- unique(final_derfindereQTL[,"gene"]) 


library(stringr)
library(tidyverse)
eQTL.unsentinalised <- list()
for (i in 1:length(regionID))
{
    print(i)
    
    chr <- gsub("ER","",unlist(str_split(regionID[i],"_"))[1])    
    
    eQTLs <- read_delim(paste0("/home/sguelfi/projects/R/hipp/data/results/derfinder/fullResults/chr",gsub("X","23",chr),"/",as.character(regionID[i])),delim = "\t")
    
    eQTLs <- eQTLs[eQTLs$FDR<=0.05,] 
    
    if(nrow(eQTLs)==0)
    {
        next()   
    }
    
    eQTL.unsentinalised[[unique(as.character(regionID[i]))]] <- as.data.frame(eQTLs)
    
}

#save(eQTL.unsentinalised,file="/home/sguelfi/projects/R/hipp/data/results/final_unsentinalised_derfinder_eQTL.rda")
load("/home/sguelfi/projects/R/hipp/data/results/final_unsentinalised_derfinder_eQTL.rda")


eQTL.unsentinalised.tab <- do.call(rbind,eQTL.unsentinalised)
load("/data/references/STOPGAP/gwas.RData")

table(duplicated(paste0(gwas.data$rsid,gwas.data$PUBMEDID)))
head(sort(table(paste0(gwas.data$rsid,gwas.data$PUBMEDID)),decreasing = T))

tmp.duplicated <- gwas.data[duplicated(paste0(gwas.data$rsid,gwas.data$PUBMEDID)),]
tmp.duplicated <- gwas.data[match(paste0(tmp.duplicated$rsid,tmp.duplicated$PUBMEDID),table = paste0(gwas.data$rsid,gwas.data$PUBMEDID)),]

ann.tmp <- cbind(eQTL.unsentinalised.tab,ann.reg[as.character(eQTL.unsentinalised.tab$gene),])
#ann.tmp <- cbind(ann.tmp,rsid=snp.map[ann.tmp$SNP,"rsid"])

intron <- ann.tmp[which(ann.tmp$exon==0 & ann.tmp$intergenic==0 & ann.tmp$intron==1),]
exon <- ann.tmp[which(ann.tmp$exon==1 & ann.tmp$intergenic==0 & ann.tmp$intron==0),]
intergenic <- ann.tmp[which(ann.tmp$exon==0 & ann.tmp$intergenic==1 & ann.tmp$intron==0),]
novel <- rbind(intron,intergenic)

novel.rs <- as.character(snp.map[match(unique(as.character(novel$SNP)),snp.map$position),"rsid"])
intron.rs <- as.character(snp.map[match(unique(as.character(intron$SNP)),snp.map$position),"rsid"])
exon.rs <- as.character(snp.map[match(unique(as.character(exon$SNP)),snp.map$position),"rsid"])

gwas.data <- gwas.data[gwas.data$pvalue<=5e-8,]


table(gwas.data$msh.cat)
# Neurological/behavioral 4864
nrow(gwas.data)
#84905-4864 = 80041

gwas.neuro <- gwas.data[which(gwas.data$msh.cat=="Neurological/behavioral"),]
nrow(gwas.neuro)
gwas.all <- gwas.data[-which(gwas.data$msh.cat=="Neurological/behavioral"),]
nrow(gwas.all)

novel.neuro <- intersect(novel.rs,as.character(unique(gwas.neuro$rsid)))
length(novel.neuro)
intron.neuro <- intersect(intron.rs,as.character(unique(gwas.neuro$rsid)))
length(intron.neuro)
exon.neuro <- intersect(exon.rs,as.character(unique(gwas.neuro$rsid)))
length(exon.neuro)

novel.all <- intersect(novel.rs,as.character(unique(gwas.all$rsid)))
length(novel.all)
intron.all <- intersect(intron.rs,as.character(unique(gwas.all$rsid)))
length(intron.all)
exon.all <- intersect(exon.rs,as.character(unique(gwas.all$rsid)))
length(exon.all)

## novel
539/4864
7860/80041

## intronic
498/4864
7732/80041

##
325/4864
5519/80041

fisher.test(matrix(c(539,7860,4864,80041),ncol=2))
fisher.test(matrix(c(325,5519,4864,80041),ncol=2))

tmp <- as.data.frame(cbind(type=c("Exons","Novel\ntranscription"),enrichment=c((325/(325+5519))/(4864/(4864+80041)),(539/(539+7860))/(4864/(4864+80041)))))
tmp$enrichment <- as.numeric(as.character((tmp$enrichment)))

define_ggplot_theme <- function(vertical_x_lab, minimal){
    
    if(minimal){
        
        theme_DZ <- 
            theme_minimal()
        
    }else{
        
        theme_DZ <- 
            theme()
        
    }
    
    theme_DZ <- 
        theme_DZ +
        theme(text = element_text(color = "black", family = "sans"),
              plot.margin = margin(1, 1, 1, 1, "cm"),
              axis.text.x = element_text(face = "italic", 
                                         size = rel(2)),
              axis.text.y = element_text(face = "italic", 
                                         size = rel(1.5), 
                                         hjust = 0.5),
              axis.ticks.x = element_line(colour = "black"), 
              axis.ticks.y = element_line(colour = "black"), 
              axis.title.x = element_text(size = rel(2.2),vjust = -2), 
              axis.title.y = element_text(size = rel(2.2), 
                                          margin = margin(t = 0, r = 10, b = 0, l = 0)),
              axis.line = element_line(color = "black", 
                                       size = 0.5),
              legend.text = element_text(size = rel(1)),
              legend.title = element_text(size = rel(1), 
                                          face = "bold"), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank()
        )
    
    if(vertical_x_lab){
        
        theme_DZ <- 
            theme_DZ +
            theme(axis.text.x = 
                      element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1))
        
    }
    
    return(theme_DZ)
    
}

generate_figure_5A <- function(figure_5A_data, theme_5A){
    
    
    col.bar <- c(pal_npg("nrc")(4)[3],  
                 pal_npg("nrc")(4)[2],
                 pal_npg("nrc")(4)[3],
                 pal_npg("nrc")(4)[1])
    
    names(col.bar) <- figure_5A_data$type
    
    figure_5A <- ggplot(figure_5A_data, aes(y=enrichment, x=type, fill=paste(type))) +
        geom_bar(position="dodge", stat="identity",color="white")+
        scale_x_discrete(name = "eQTL type")+
        scale_y_continuous(name = "Enrichment ratio")+
        scale_fill_manual(values = col.bar)+
        labs(fill="")+
        ylim(c(0,1.2))+
        geom_hline(yintercept=1, linetype="dashed", size=1)#+
        #theme_5A
    
    return(figure_5A)
}

# Main ------------------------------------------------------------------------------------------------

theme_5A <- define_ggplot_theme(vertical_x_lab = F, minimal = F)

figure_5A <- generate_figure_5A(tmp, theme_5A)


