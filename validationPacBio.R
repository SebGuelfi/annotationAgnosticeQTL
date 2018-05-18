## annotation.

library(rtracklayer)

brain.pacBio <- import.gff("/home/sguelfi/pacBio/brain.all_size.5merge.collapsed.gff")

load("/home/sguelfi/projects/R/hipp/data/results/final_derfinder.rda")

brain.pacBio <- brain.pacBio[as.character(brain.pacBio$type)=="exon",]

ann.reg <- cbind(ann.reg,regNames=rownames(ann.reg))
intron <- ann.reg %>% filter(intron>=1 & intergenic==0 & exon==0)
intergenic <- ann.reg %>% filter(intron==0 & intergenic>=1 & exon==0)
exonic <- ann.reg%>% filter(intron==0 & intergenic==0 & exon>=1)

expr.gr.exonic <- expr.gr[as.character(exonic$regNames)]
expr.gr.intronic <- expr.gr[as.character(intron$regNames)]
expr.gr.intergenic <- expr.gr[as.character(intergenic$regNames)]

seqlevels(expr.gr.exonic) <- paste0('chr',seqlevels(expr.gr.exonic))
seqlevels(expr.gr.intronic) <- paste0('chr',seqlevels(expr.gr.intronic))
seqlevels(expr.gr.intergenic) <- paste0('chr',seqlevels(expr.gr.intergenic))

(table(countOverlaps(expr.gr.exonic,brain.pacBio,type="within")>0)/length(expr.gr.exonic))*100
(table(countOverlaps(expr.gr.intronic,brain.pacBio,type="within")>0)/length(expr.gr.intronic))*100
(table(countOverlaps(expr.gr.intergenic,brain.pacBio,type="within")>0)/length(expr.gr.intergenic))*100


brain.pacBio <- import.gff("/home/sguelfi/pacBio/IsoSeq_Alzheimer_2016edition_polished.confident.fusion.hg38.gff")

load("/home/sguelfi/projects/R/hipp/data/results/final_derfinder.rda")

brain.pacBio <- brain.pacBio[as.character(brain.pacBio$type)=="exon",]

ann.reg <- cbind(ann.reg,regNames=rownames(ann.reg))
intron <- ann.reg %>% filter(intron>=1 & intergenic==0 & exon==0)
intergenic <- ann.reg %>% filter(intron==0 & intergenic>=1 & exon==0)
exonic <- ann.reg%>% filter(intron==0 & intergenic==0 & exon>=1)

expr.gr.exonic <- expr.gr[as.character(exonic$regNames)]
expr.gr.intronic <- expr.gr[as.character(intron$regNames)]
expr.gr.intergenic <- expr.gr[as.character(intergenic$regNames)]

seqlevels(expr.gr.exonic) <- paste0('chr',seqlevels(expr.gr.exonic))
seqlevels(expr.gr.intronic) <- paste0('chr',seqlevels(expr.gr.intronic))
seqlevels(expr.gr.intergenic) <- paste0('chr',seqlevels(expr.gr.intergenic))

(table(countOverlaps(expr.gr.exonic,brain.pacBio,type="within")>0)/length(expr.gr.exonic))*100
(table(countOverlaps(expr.gr.intronic,brain.pacBio,type="within")>0)/length(expr.gr.intronic))*100
(table(countOverlaps(expr.gr.intergenic,brain.pacBio,type="within")>0)/length(expr.gr.intergenic))*100






