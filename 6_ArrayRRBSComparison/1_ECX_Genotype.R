#!/usr/bin/env Rscript
## ----------Script-----------------  
##
## Purpose: Comparison of Array and RRBS results - ECX at Genotype level
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------
## 
## Entorinal cortex: RRBS and Array, rTg4510 and J20 
## Model: 
##      Array: E.Walker - MixedModelResultsECX
##      RRBS: I.Castanho - DMPsInteractionModel_rTg4510_sig_genotype
## Annotations with ChipSeeker (keeping original Symbol column):
##      Array: E.Walker using 1500bp 
##      RRBS: S.Leung using same function as E.Walker
## Comparison at significant sites (FDR-adjusted < 0.05) 
##      Position, Gene (Chipseeker annotation)
##

#-------------- input -------------

# source functions to run script
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/4_DMP_annotation/0_source_functions.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/5_Array/5a_MEM_analysisFunctions.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/methylation.config.R")

# results after applying BiSeq and smoothing
# load from two mouse models and save into list
RRBS_complete <- list()
load(file = paste0(dirnames$biseq, "/rTg4510_RRBSbetasComplete.RData"))
RRBS_complete$rTg4510 <- RRBS_completebetas
load(file = paste0(dirnames$biseq, "/J20_RRBSbetasComplete.RData"))
RRBS_complete$J20 <- RRBS_completebetas
RRBS_complete <- lapply(RRBS_complete, function(x) x %>% tibble::rownames_to_column(., var = "position"))

# metadata of the combined phenotype and beta values
Array_complete <- list(
  rTg4510 = PrepBetasPhenoFiles(raw$array$rTg4510Array, "rTg4510"),
  J20 = PrepBetasPhenoFiles(raw$array$J20Array, "J20")
)

#-------------- overlap of significant sites ------------
sig <- list(
  array = list(
    rTg4510_ECX_geno = anno$array$rTg4510_ECX_geno %>% mutate(FDR_adj_genotype = p.adjust(PrZ.GenotypeTG , method = "fdr"), Position = position, ChIPseeker_GeneSymbol = SYMBOL),
    J20_ECX_geno = anno$array$J20_ECX_geno %>% mutate(FDR_adj_genotype = p.adjust(PrZ.GenotypeTG , method = "fdr"), Position = position, ChIPseeker_GeneSymbol = SYMBOL)
  ),
  rrbs = list(
    rTg4510_ECX_geno = anno$rrbs$rTg4510_ECX_geno,
    J20_ECX_geno = anno$rrbs$J20_ECX_geno
  )
)
sig$array <- lapply(sig$array, function(x) x %>% filter(FDR_adj_genotype < 0.05))
sig$rrbs <- lapply(sig$rrbs, function(x) x %>% filter(FDR_adj_genotype < 0.05))

# number of significant sites
# rrbs
length(sig$rrbs$rTg4510_ECX_geno$Position)
length(sig$rrbs$J20_ECX_geno$Position)
# array
length(sig$array$rTg4510_ECX_geno$Position)
length(sig$array$J20_ECX_geno$Position)

# number of significant genes
length(unique(sig$rrbs$rTg4510_ECX_geno$ChIPseeker_GeneSymbol))
length(unique(sig$rrbs$J20_ECX_geno$ChIPseeker_GeneSymbol))
length(unique(sig$array$rTg4510_ECX_geno$ChIPseeker_GeneSymbol))
length(unique(sig$array$J20_ECX_geno$ChIPseeker_GeneSymbol))

# position level
length(intersect(sig$array$rTg4510_ECX_geno$Position, sig$rrbs$rTg4510_ECX_geno$Position))
length(intersect(sig$array$J20_ECX_geno$Position, sig$rrbs$J20_ECX_geno$Position))

# gene level
length(intersect(sig$array$rTg4510_ECX_geno$ChIPseeker_GeneSymbol, sig$rrbs$rTg4510_ECX_geno$ChIPseeker_GeneSymbol))
length(intersect(sig$array$J20_ECX_geno$ChIPseeker_GeneSymbol, sig$rrbs$J20_ECX_geno$ChIPseeker_GeneSymbol))

arrayCols = c("cpg","Position","PrZ.GenotypeTG","FDR_adj_genotype","ChIPseeker_GeneSymbol")
rrbsCols = c("Position","p.val.Genotype","FDR_adj_genotype","ChIPseeker_GeneSymbol")
sigGene <- list(
  J20 = merge(sig$array$J20_ECX_geno[arrayCols],sig$rrbs$J20_ECX_geno[rrbsCols], by = "ChIPseeker_GeneSymbol", suffixes = c("Arr","RRBS")),
  rTg4510 = merge(sig$array$rTg4510_ECX_geno[arrayCols],sig$rrbs$rTg4510_ECX_geno[rrbsCols], by = "ChIPseeker_GeneSymbol", suffixes = c("Arr","RRBS"))
)

length(unique(sigGene$rTg4510$PositionRRBS))
length(unique(sigGene$rTg4510$PositionArr))


#-------------- output --------------
for(i in 1:2){
  write.csv(sigGene[[i]], paste0(dirnames$output, names(sigGene)[[i]],"_ArrRRBS_sigGene.csv"))
}

plot_array_rrbs_gene <- function(gene,model,pheno){
  
  filtered <- sigGene[[model]] %>% filter(ChIPseeker_GeneSymbol == gene) 
  arrFiltered <- Array_complete[[model]]$betas %>% reshape2::melt() %>% filter(Var1 %in% filtered$cpg) %>% `colnames<-`(c("Cpg","Sample_ID", "Meth"))
  arrFiltered <- merge(arrFiltered,pheno[[model]][,c("Sample_ID","Genotype","Age_months")], by = "Sample_ID") 
  arrFiltered <- merge(arrFiltered,filtered[,c("cpg","PositionArr")]) %>% dplyr::rename("position"="PositionArr") %>% dplyr::select(Sample_ID,position,Meth,Genotype,Age_months) %>% 
    `colnames<-`(c("Sample","position", "Meth","Genotype","Age")) %>%
    mutate(Dataset = "array")
  
  rrbsFiltered <- RRBS_complete[[model]] %>% filter(position %in% filtered$PositionRRBS) %>% reshape2::melt() 
  rrbsFiltered <- merge(rrbsFiltered, pheno[[model]][,c("Sample_ID","Genotype","Age_months")], by.x = "variable", by.y = "Sample_ID") %>% `colnames<-`(c("Sample", "position", "Meth","Genotype","Age")) %>% mutate(Dataset = "rrbs")
  
  merged <- rbind(arrFiltered,rrbsFiltered) %>% mutate(coordinate = word(position,c(2),sep=fixed(":"))) %>% mutate(coordinate = as.numeric(coordinate))
  
  p <- ggplot(merged, aes(x = Genotype, y = Meth, colour = Dataset)) + geom_boxplot() + 
    theme_classic() +
    labs(x = "Genotype", y = "Methylation", title = paste0(gene,"\n")) +
    scale_colour_manual(values = c("black",label_colour(model))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  if(length(merged$coordinate) > 1){
    p <- p + facet_grid(~as.factor(coordinate)) 
  }
  
  return(p)
}

plot_array_rrbs_gene("Kirrel3","rTg4510",pheno)
plot_array_rrbs_gene("Plxnc1","rTg4510",pheno)
plot_array_rrbs_gene("Osbpl3","rTg4510",pheno)
plot_array_rrbs_gene("Satb1","rTg4510",pheno)
plot_array_rrbs_gene("Arx","J20",pheno)
plot_array_rrbs_gene("Sorcs1","J20",pheno)
plot_array_rrbs_gene("Bnc2","J20",pheno)

plot_RRBS <- function(gene, model, dmp){
  
  rrbsFiltered <- RRBS_complete[[model]] %>% filter(position %in% dmp) %>% reshape2::melt() 
  rrbsFiltered <- merge(rrbsFiltered, coldata[,c("Sample_ID","Genotype","Age_months")], by.x = "variable", by.y = "Sample_ID") %>% `colnames<-`(c("Sample", "position", "Meth","Genotype","Age")) %>% mutate(Dataset = "rrbs")
  
  p <- ggplot(rrbsFiltered, aes(x = Genotype, y = Meth)) + geom_boxplot() + facet_grid(~as.factor(position)) +
    #geom_jitter(aes(colour = Age))
    theme_classic() +
    labs(x = "Genotype", y = "Methylation", title = paste0(gene,"\n")) +
    scale_colour_manual(values = c("black",label_colour(model))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  return(p)
  
}

plot_RRBS("Anp32b","rTg4510",as.vector(sigResults$rTg4510.genotype[sigResults$rTg4510.genotype$ChIPseeker_GeneSymbol == "Anp32b","Position"]))
