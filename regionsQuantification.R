## author: Sebastian Guelfi
## This script quantify the regions defined using the optimisation step
## 1. Select the cutoff and maxgap based on the optimisation 
## 2. Get the coverage for UKBEC hippocampus
## 3. 

library(tidyverse)
library(derfinder)
library(GenomicRanges)

load("/home/sguelfi/projects/R/OMIM/outfile/ER_annotation_exon_delta_details_hipp_cutoff_0.2_to_15_by_0.1_40M_maxgap.rda")

get_min_exon_delta_max_exon_delta_eq_0 <- function(exon_details){
  
  # get maxgap for the highest num exon delta = 0
  return(exon_details %>% 
    filter(exon_delta_median == min(exon_delta_median)) %>% 
    filter(num_exon_delta_eq_0 == max(num_exon_delta_eq_0)))
}

exon_details <- get_min_exon_delta_max_exon_delta_eq_0(ER_annotation_exon_delta_details_all_maxgap)

# extract the regions 
load(paste0("/home/sguelfi/projects/R/OMIM/outfile//maxgap_normalised_40M_reads/region_cuts_maxgap_",
            exon_details$maxgap,".rda"))

regs <- region_cuts_maxGap[[as.character(exon_details$mean_coverage_cut_off)]]

regs <- regs[width(regs)>2]

no_cores <- 3 
library(foreach)
library(doParallel)

cl<-makeCluster(no_cores)
clusterExport(cl,c("getRegionCoverage", "rbind", "colSums", "lapply","do.call",
                      "seqnames","start","end"))
registerDoParallel(cl)

getCoverage_par <- function(i,L,gr){
    load(paste0("/home/sguelfi/projects/R/hipp/data/derfinder/fullCoverage/derfinder/fullCoverageChr",i,".rda"))
    #gr <- regs[seqnames(regs)==paste0("chr",i)]
    #gr <- gr[1:20]
    coverage <- getRegionCoverage(fullCov = fullCov, gr, totalMapped=NULL)
    covMat <- lapply(coverage, colSums)
    covMat <- do.call(rbind, covMat)
    covMat <- covMat / L
    rownames(covMat) <- paste0(seqnames(gr),"_",start(gr),"_",end(gr))
    
    save(covMat,file=paste0("/home/sguelfi/projects/R/hipp/data/expression/derfinder/coverage/chr",i,".rda"))
    return(TRUE)
}
    
out <- foreach(i=c(1:22,"X"),.combine = c)  %dopar%  getCoverage_par(i,100,regs[seqnames(regs)==paste0("chr",i)])

stopCluster(cl)
rm(cl)


counts <-NULL
for(i in list.files("/home/sguelfi/projects/R/hipp/data/expression/derfinder/coverage/"))
{
    load(paste0("/home/sguelfi/projects/R/hipp/data/expression/derfinder/coverage/",i))
    counts <- rbind(counts,covMat)
}

names(regs)<- paste0(seqnames(regs),"_",start(regs),"_",end(regs))
regs <- regs[rownames(counts)]


source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_txDb_from_gtf.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_genomic_state.R")

ensembl_grch38_v87_TxDb <- generate_txDb_from_gtf(gtf_gff3_path = "/data/references/ensembl/gtf_gff3/v87/Homo_sapiens.GRCh38.87.gtf",
                           output_path = "/data/references/ensembl/txdb_sqlite/v87/ensembl_grch38_v87_txdb.sqlite",
                           seq_levels_to_keep = c(1:22, "X", "Y", "MT"), 
                           genome_build = "hg38")

genomic_state_ensembl_grch38_v87 <- 
    generate_genomic_state(ensembl_grch38_v87_TxDb, 
                           output_path = "/data/references/ensembl/genomic_state/v87/ensembl_grch38_v87_genomic_state.rda", 
                           chrs_to_keep = str_c("chr", c(1:22, "X", "Y", "M")))


regs.ann <- annotateRegions(regions = regs, 
                    genomicState = genomic_state_ensembl_grch38_v87$fullGenome, 
                    annotate = T)



save(regs,counts,regs.ann,file="/home/sguelfi/projects/R/hipp/data/expression/derfinder/mergedDerfinder.rda")

