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
  # entorhinal cortex
  rTg4510_genotype_ECX = unique(sigRes$rTg4510$Genotype %>% filter(FDR_adj_Genotype < 0.05) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  rTg4510_interaction_ECX = unique(sigRes$rTg4510$GenotypeAge %>% filter(FDR_adj_GenotypeAge < 0.05) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  rTg4510_pathology_ECX = unique(sigRes$rTg4510$PathologyCommonInteraction %>% filter(FDR_adj_Pathology < 0.05) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  J20_genotype_ECX = unique(sigRes$J20$Genotype %>% filter(FDR_adj_Genotype < 9e-8) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  J20_interaction_ECX = unique(sigRes$J20$GenotypeAge %>% filter(FDR_adj_GenotypeAge < 0.05) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  J20_pathology_ECX = unique(sigRes$J20$PathologyCommonInteraction %>% filter(FDR_adj_Pathology < 0.05) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  # hippocampus
  rTg4510_genotype_HIP = unique(sigResArrayHIP$rTg4510$Genotype %>% filter(FDR_adj_Genotype < 0.05) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  J20_genotype_HIP = unique(sigResArrayHIP$J20$Genotype %>% filter(FDR_adj_Genotype < 0.05) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  rTg4510_pathology_HIP = unique(sigResArrayHIP$rTg4510$PathologyCommonInteraction %>% filter(FDR_adj_Pathology < 0.05) %>% dplyr::select(ChIPseeker_GeneSymbol)),
  J20_pathology_HIP = unique(sigResArrayHIP$J20$PathologyCommonInteraction %>% filter(FDR_adj_Pathology < 0.05) %>% dplyr::select(ChIPseeker_GeneSymbol))
)


#-------------- write out gene list -------------

for(i in 1:length(mouseDMP)){
  name = names(mouseDMP [i])
  write.table(mouseDMP[[i]], paste0(dirnames$humanAnnot,name,"_geneList.txt"),row.names=F,col.names=F,sep="\t",quote=F)
}
