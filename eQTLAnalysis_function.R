

eQTL_analysis_derfinder <- function(chr,expr.gr,expr.cqn,pathImputed,covs,outputFolder)
{"/home/seb/projectsR/hipp/data/"
    ##i <- 21  test ## loop for chromosome
    
    system(paste0("mkdir ",outputFolder,"chr",chr))
    if (chr==23)
    {
        expr.tmp.gr <- expr.gr[seqnames(expr.gr) == "chrX"]
        dosage <- read.delim(paste0(pathImputed,"X.dosage.traw"),check.names = F)
        rownames(dosage) <- paste("chrX",dosage$POS,paste0(dosage$COUNTED,"_",dosage$ALT),sep  = ":")
        dosage.gr <- GRanges("chrX",IRanges(dosage$POS,dosage$POS))
    }else{
        expr.tmp.gr <- expr.gr[gsub("chr","",seqnames(expr.gr)) == chr]
        dosage <- read.delim(paste0(pathImputed,chr,".dosage.traw"),check.names = F)
        rownames(dosage) <- paste(paste0("chr",dosage$CHR),dosage$POS,paste0(dosage$COUNTED,"_",dosage$ALT),sep  = ":")
        dosage.gr <- GRanges(paste0("chr",dosage$CHR),IRanges(dosage$POS,dosage$POS))
    }
    
    ## j <-1 ## loop for the gene
    results <- NULL
    for(j in 1:length(expr.tmp.gr))
    {
        print(j)
        expr.tmp <- t(as.matrix(expr.cqn[names(expr.tmp.gr)[j],]))
        rownames(expr.tmp) <- names(expr.tmp.gr)[j]
        dosage.tmp <- as.matrix(dosage[as.matrix(findOverlaps(expr.tmp.gr[j],subject = dosage.gr,maxgap = 1000000))[,2],colnames(expr.tmp)])
        if(nrow(dosage.tmp)>0)
        {
            stopifnot(identical(colnames(dosage.tmp),colnames(expr.tmp)))
            covs.tmp <- as.matrix(t(covs[colnames(expr.tmp),]))
            stopifnot(identical(colnames(covs.tmp),colnames(expr.tmp)))
            
            library(MatrixEQTL)
            
            my.expr <-SlicedData$new()
            my.expr$CreateFromMatrix(expr.tmp)
            my.markers <- SlicedData$new()
            my.markers$CreateFromMatrix(dosage.tmp)
            my.cov <- SlicedData$new()
            my.cov$CreateFromMatrix(covs.tmp)
            
            ## load genetic pca
            
            #my.cov <- SlicedData$new()
            #my.cov$CreateFromMatrix(cov)
            ##rm(expr.tmp, dosage.tmp)
            store <- Matrix_eQTL_main( my.markers, my.expr,my.cov, output_file_name = paste0(outputFolder,"/","chr",chr,"/",names(expr.tmp.gr)[j]),pvOutputThreshold=1, useModel=modelLINEAR, errorCovariance=numeric(0), verbose=T )
            rank <- 0
            if(nrow(store$all$eqtls[store$all$eqtls$FDR <=0.05,])>0){
                print(j)
                results.tmp <- store$all$eqtls[which(store$all$eqtls$pvalue == min(store$all$eqtls$pvalue)),]
                results.tmp$rank <- rank
                results <- rbind(results,results.tmp)
                rm(results.tmp)
                cond_loop <- T
                while(cond_loop){
                    lead.SNP <- store$all$eqtls[1,1]
                    covs.tmp <- rbind(covs.tmp,dosage.tmp[match(as.character(lead.SNP),rownames(dosage.tmp)),colnames(covs.tmp)])
                    my.cov <- SlicedData$new()
                    my.cov$CreateFromMatrix(covs.tmp)
                    my.markers <- SlicedData$new()
                    dosage.tmp <- dosage.tmp[-match(as.character(lead.SNP),rownames(dosage.tmp)),]
                    if(nrow(dosage.tmp)>0)
                    {
                        my.markers$CreateFromMatrix(dosage.tmp)
                        store <- Matrix_eQTL_main( my.markers, my.expr,my.cov, output_file_name = NULL,pvOutputThreshold=0.05, useModel=modelLINEAR, errorCovariance=numeric(0), verbose=T )
                        if(nrow(store$all$eqtls[store$all$eqtls$FDR <=0.05,])>0){
                            rank=rank+1
                            results.tmp <- store$all$eqtls[which(store$all$eqtls$pvalue == min(store$all$eqtls$pvalue)),]
                            results.tmp$rank <- rank
                            results <- rbind(results,results.tmp)
                            rm(results.tmp)
                        }
                        else{
                            cond_loop <- F
                        }
                    }else{
                        cond_loop <- F
                    }
                }
                
            }
        }else{
            write.table(names(expr.tmp.gr)[j],file=paste0(outputFolder,"/noVariantInRegion.txt"),append = T,quote = F,row.names = F,col.names = F)
            
        }
    }
    return(results)
}



## function to validate the eQTL from NABEC
eQTL_analysis_derfinder_NABEC <- function(chr,expr.gr,expr.cqn,pathImputed,covs,outputFolder)
{"/home/seb/projectsR/hipp/data/"
    ##i <- 21  test ## loop for chromosome
    system(paste0("mkdir ",outputFolder,"/","chr",chr))
    expr.tmp.gr <- expr.gr[gsub("chr","",seqnames(expr.gr)) == chr]
    dosage <- read.delim(paste0(pathImputed,"chr",chr,".dosage.traw"),check.names = F)
    dosage <- dosage[!duplicated(dosage$SNP),]
    rownames(dosage) <- paste(dosage$SNP)
    dosage.gr <- GRanges(as.character(unlist(lapply(strsplit(as.character(dosage$SNP),split = ":",fixed=T),function(x){return(x[1])}))),
                         IRanges(as.numeric(unlist(lapply(strsplit(as.character(dosage$SNP),split = ":",fixed=T),function(x){return(x[2])}))),
                                 as.numeric(unlist(lapply(strsplit(as.character(dosage$SNP),split = ":",fixed=T),function(x){return(x[2])})))))
    dosage <- dosage[,-c(1:6)]
    
    ## j <-1 ## loop for the gene
    results <- NULL
    for(j in 1:length(expr.tmp.gr))
    {
        print(j)
        expr.tmp <- t(as.matrix(expr.cqn[names(expr.tmp.gr)[j],]))
        rownames(expr.tmp) <- names(expr.tmp.gr)[j]
        
        colnames(dosage) <- gsub("_","-",gsub("UKY_","UKY-",gsub("UMARY_","UMARY-",gsub("SH_","SH-",colnames(dosage)))))
        dosage.tmp <- as.matrix(dosage[as.matrix(findOverlaps(expr.tmp.gr[j],subject = dosage.gr,maxgap = 1000000))[,2],colnames(expr.tmp)])
        if(nrow(dosage.tmp)>0)
        {
            stopifnot(identical(colnames(dosage.tmp),colnames(expr.tmp)))
            covs.tmp <- as.matrix(t(covs[colnames(expr.tmp),]))
            stopifnot(identical(colnames(covs.tmp),colnames(expr.tmp)))
            
            library(MatrixEQTL)
            
            my.expr <-SlicedData$new()
            my.expr$CreateFromMatrix(expr.tmp)
            my.markers <- SlicedData$new()
            my.markers$CreateFromMatrix(dosage.tmp)
            my.cov <- SlicedData$new()
            my.cov$CreateFromMatrix(covs.tmp)
            
            ## load genetic pca
            
            #my.cov <- SlicedData$new()
            #my.cov$CreateFromMatrix(cov)
            ##rm(expr.tmp, dosage.tmp)
            store <- Matrix_eQTL_main( my.markers, my.expr,my.cov, output_file_name = paste0(outputFolder,"/","chr",chr,"/",names(expr.tmp.gr)[j]),pvOutputThreshold=1, useModel=modelLINEAR, errorCovariance=numeric(0), verbose=T )
            rank <- 0
            if(nrow(store$all$eqtls[store$all$eqtls$FDR <=0.05,])>0){
                print(j)
                results.tmp <- store$all$eqtls[which(store$all$eqtls$pvalue == min(store$all$eqtls$pvalue)),]
                results.tmp$rank <- rank
                results <- rbind(results,results.tmp)
                rm(results.tmp)
                cond_loop <- T
                while(cond_loop){
                    lead.SNP <- store$all$eqtls[1,1]
                    covs.tmp <- rbind(covs.tmp,dosage.tmp[match(as.character(lead.SNP),rownames(dosage.tmp)),colnames(covs.tmp)])
                    my.cov <- SlicedData$new()
                    my.cov$CreateFromMatrix(covs.tmp)
                    my.markers <- SlicedData$new()
                    dosage.tmp <- dosage.tmp[-match(as.character(lead.SNP),rownames(dosage.tmp)),]
                    if(nrow(dosage.tmp)>0)
                    {
                        my.markers$CreateFromMatrix(dosage.tmp)
                        store <- Matrix_eQTL_main( my.markers, my.expr,my.cov, output_file_name = NULL,pvOutputThreshold=0.05, useModel=modelLINEAR, errorCovariance=numeric(0), verbose=T )
                        if(nrow(store$all$eqtls[store$all$eqtls$FDR <=0.05,])>0){
                            rank=rank+1
                            results.tmp <- store$all$eqtls[which(store$all$eqtls$pvalue == min(store$all$eqtls$pvalue)),]
                            results.tmp$rank <- rank
                            results <- rbind(results,results.tmp)
                            rm(results.tmp)
                        }
                        else{
                            cond_loop <- F
                        }
                    }else{
                        cond_loop <- F
                    }
                }
                
            }
        }else{
            write.table(names(expr.tmp.gr)[j],file=paste0(outputFolder,"/noVariantInRegion.txt"),append = T,quote = F,row.names = F,col.names = F)
            
        }
    }
    return(results)
}

