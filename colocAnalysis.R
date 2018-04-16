
## Coloc analysis

## Pass the eQTL regions to David

## need to get the rsIDs first


# snp.map <- NULL
# for( i in c(1:22,"X"))
# {
#     print(i)
#     tmp <- read.table(pipe(paste0("cut -f1,2,4,5,6 /home/sguelfi/projects/R/hipp/data/imputed_data/imputed_genotype_hg38/",i,".dosage.traw")), header=T)    
#     position <- paste0(tmp$CHR,":", tmp$POS,":",tmp$COUNTED,"_",tmp$ALT)
#     tmp <- cbind(position,rsid=as.character(tmp$SNP))
#     snp.map <- rbind(snp.map,tmp)
#     rm(tmp,position)
# }
# 
# snp.map <- as.data.frame(snp.map)
## remove duplicates
# snp.map <- snp.map[-which(duplicated(snp.map$position)),]
# rownames(snp.map) <- snp.map$position
# save(snp.map,file="/home/sguelfi/projects/R/hipp/data/general/snp.map.rda")

load("/home/sguelfi/projects/R/hipp/data/results/final_derfinder.rda")
load(file="/home/sguelfi/projects/R/hipp/data/general/snp.map.rda")

## we first collect the transcript eQTL
#head(final_derfindereQTL)

regionID <- unique(final_derfindereQTL[,"gene"]) 


library(stringr)
library(tidyverse)
for (i in 1:length(regionID))
{
    print(i)
    
    chr <- gsub("ER","",unlist(str_split(regionID[i],"_"))[1])    
    
    eQTLs <- read_delim(paste0("/home/sguelfi/projects/R/hipp/data/results/derfinder/fullResults/chr",chr,"/",as.character(regionID[i])),delim = "\t")
    
    tmpf <- tempfile()
    keepRows <- c(unique(as.character(c(unlist(lapply(str_split(eQTLs[,"SNP"],":"),function(x){return(x[2])}))))))
    write.table(data.frame(sort(keepRows)), file=tmpf, row.names=F, col.names=F, quote=F )
    filename <- paste0("/home/sguelfi/projects/R/hipp/data/imputed_data/imputed_genotype_hg38/",chr,".frq")
    
    
    maf <- read.table( pipe( paste("fgrep -w -f", tmpf, filename) ), header=F)
    ## remove duplicates
    maf <- maf[!duplicated(maf$V2),]
    rownames(maf) <- maf$V2
    
    maf.tmp <- maf[as.character(unlist(lapply(str_split(eQTLs[,"SNP"],":"),function(x){return(x[2])}))),"V6"]
    eQTLs[,"maf"] <-  ifelse(maf.tmp>0.5,1-maf.tmp,maf.tmp)
    eQTLs[,"rsid"] <- as.character(snp.map[as.character(eQTLs$SNP),"rsid"])
    rm(maf.tmp,maf)
    
    save(eQTLs,file=paste0("/home/sguelfi/projects/R/hipp/data/results/derfinder/coloc/",regionID[i],".rda"))
}

tail(grep("ER1_",regionID))


## Coloc analysis

## Pass the eQTL regions to David

## need to get the rsIDs first


snp.map <- NULL
for( i in c(1:22,"X"))
{
    print(i)
    tmp <- read.table(pipe(paste0("cut -f1,2,4,5,6 /home/sguelfi/projects/R/hipp/data/imputed_genotype_v2/imputed_genotype_hg38_95/",i,".dosage.traw")), header=T)
    position <- paste0(tmp$CHR,":", tmp$POS,":",tmp$COUNTED,"_",tmp$ALT)
    tmp <- cbind(position,rsid=as.character(tmp$SNP))
    snp.map <- rbind(snp.map,tmp)
    rm(tmp,position)
}
rm(i)
snp.map <- as.data.frame(snp.map)
# remove duplicates there are no duplicates
# snp.map <- snp.map[-which(duplicated(snp.map$position)),]
rownames(snp.map) <- snp.map$position
save(snp.map,file="/home/sguelfi/projects/R/hipp/data/general/snp.map.rda")

load("/home/sguelfi/projects/R/hipp/data/results/final_derfinder.rda")
load(file="/home/sguelfi/projects/R/hipp/data/general/snp.map.rda")

## we first collect the transcript eQTL
#head(final_derfindereQTL)

regionID <- unique(final_derfindereQTL[,"gene"]) 


library(stringr)
library(tidyverse)
for (i in 1:length(regionID))
{
    print(i)
    
    chr <- gsub("ER","",unlist(str_split(regionID[i],"_"))[1])    
    
    eQTLs <- read_delim(paste0("/home/sguelfi/projects/R/hipp/data/results/derfinder/fullResults/chr",gsub("X","23",chr), "/",as.character(regionID[i])),delim = "\t")
    
    tmpf <- tempfile()
    keepRows <- c(unique(as.character(c(unlist(lapply(str_split(eQTLs$SNP,":"),function(x){return(x[2])}))))))
    write.table(data.frame(sort(keepRows)), file=tmpf, row.names=F, col.names=F, quote=F )
    filename <- paste0("/home/sguelfi/projects/R/hipp/data/imputed_genotype_v2/imputed_genotype_hg38_95/",chr,".frq")
        
    maf <- read.table(pipe( paste("fgrep -w -f", tmpf, filename) ), header=F)
    ## remove duplicates
    maf <- maf[!duplicated(maf$V2),]
    rownames(maf) <- maf$V2
    
    maf.tmp <- maf[as.character(unlist(lapply(str_split(eQTLs$SNP,":"),function(x){return(x[2])}))),"V6"]
    eQTLs$maf <-  ifelse(maf.tmp>0.5,1-maf.tmp,maf.tmp)
    eQTLs$rsid <- as.character(snp.map[gsub("X","23",as.character(eQTLs$SNP)),"rsid"])
    rm(maf.tmp,maf)
    
    save(eQTLs,file=paste0("/home/sguelfi/projects/R/hipp/data/results/derfinder/coloc_v2/",regionID[i],".rda"))
}



