#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Generate gene lists associated to DMP from rTg4510 and J20 RRBS and genotype
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------
## 
## 1. EWalker ran array analysis for rTg4510 and J20 (note missing interaction for J20 array)
## 2. SLeung ran DMP analysis for rTg4510 and J20
## 2. SLeung obtained gene list from rTg4510 and J20 RRBS analysis for genes with differentially methylated probes (using chipseeker)


## ------------ input ------------

# directory input
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "0_PaperFigs/paper_import.config.R"))


#-------------- significant DMP -------------

# list of genes with significant DMP in rTg4510 and J20 
# note not split by platform
# significance defined by gwas threshold 
mouseDMP <- list(
  rTg4510_genotype = unique(sigRes$rTg4510$Genotype %>% filter(P.value_Genotype < 9e-8) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  rTg4510_pathology = unique(sigRes$rTg4510$Pathology %>% filter(P.value_Pathology < 9e-8) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  J20_genotype = unique(sigRes$J20$Genotype %>% filter(P.value_Genotype < 9e-8) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  J20_pathology = unique(sigRes$J20$Pathology %>% filter(P.value_Pathology < 0.05) %>% dplyr::select(ChIPseeker_GeneSymbol))
)


#-------------- write out gene list -------------

for(i in 1:length(mouseDMP)){
  name = names(mouseDMP [i])
  write.table(mouseDMP[[i]], paste0(dirnames$humanAnnot,name,"_geneList.txt"),row.names=F,col.names=F,sep="\t",quote=F)
}