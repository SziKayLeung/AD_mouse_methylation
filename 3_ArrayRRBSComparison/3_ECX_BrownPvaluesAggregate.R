#!/usr/bin/env Rscript
## ----------Script-----------------  
##
## Purpose: Apply brown's method to aggregate P-values from differentially expressed results
##         
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------


#-------------- packages -------------

suppressMessages(library("EmpiricalBrownsMethod"))
suppressMessages(library("ggbio"))

#-------------- input -------------

rootDir = "/lustre/projects/Research_Project-191406"
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/"
source(paste0(scriptDir, "import.config"))

sigBeta <- list(
  rTg4510 = get(load(paste0(dirnames$annotated,"/final/rTg4510_ECX_sigBeta.RData"))),
  J20 = get(load(paste0(dirnames$annotated,"/final/J20_ECX_sigBeta.RData")))
);  rm(rTg4510_sig_betaFull); rm(J20_sig_betaFull)

sigRes <- list(
  rTg4510 = get(load(file = paste0(dirnames$annotated,"/final/rTg4510_ECX_sigResultsDMPs.RData"))),
  J20 = get(load(paste0(dirnames$annotated,"/final/J20_ECX_sigResultsDMPs.RData")))
);  rm(rTg4510_ECX_sigResultsFull); rm(J20_ECX_sigResultsFull)



#-------------- function -------------

brownAggregationPvalues <- function(betaMatrix, res, gene, testCol){
  
  #print(gene)
  # extract the significant probes annotated to gene
  resGene <- res %>% filter(ChIPseeker_GeneSymbol == gene)
  #print(resGene)
  #print(resGene$Position)
  
  # extract only the relevant positions
  betaGene <- betaMatrix[resGene$Position,]
  
  # replace NAN with NA to impute mean
  betaGene[betaGene == "NaN"] <- NA
  
  # impute mean for NA values
  impute_mean <- function(x) {
    na_index <- is.na(x)
    x[na_index] <- mean(x, na.rm = TRUE)
    return(x)
  }
  
  betaGene_imputed <- as.data.frame(t(apply(betaGene, 1, impute_mean)))
  brownPvalue <- empiricalBrownsMethod(data_matrix=betaGene_imputed, p_values=resGene[[testCol]], extra_info=TRUE)$P_test
  return(c(gene, brownPvalue))
}

apply_brownAggregation <- function(betaMatrix, res, testCol){
  
  geneList <- unique(na.omit(res$ChIPseeker_GeneSymbol))
  resBrownResList <- lapply(geneList, function(x) 
    brownAggregationPvalues(betaMatrix, res, x, testCol)
  )
  
  resBrownRes <- do.call(rbind.data.frame, resBrownResList)
  colnames(resBrownRes) <- c("Gene","Brown_Pvalue")
  
  resBrownRes <- resBrownRes %>% arrange(as.numeric(as.character(Brown_Pvalue)))
  return(resBrownRes)
}

merge_brownAggregation <- function(sigBeta, sigRes){
  
  dat_sig_brown <- list(
    Genotype = apply_brownAggregation(sigBeta$Genotype, sigRes$Genotype, "FDR_adj_Genotype"),
    Pathology = apply_brownAggregation(sigBeta$Pathology, sigRes$Pathology, "FDR_adj_Pathology")
  )
  
  # merge genotype and pathology output into one table
  dat_sig_brown_merged <- merge(dat_sig_brown$Genotype, dat_sig_brown$Pathology, by = "Gene", all = T)
  colnames(dat_sig_brown_merged) <- c("Gene", "BrownPvalue_Genotype", "BrownPvalue_Pathology")
  dat_sig_brown_merged  <- dat_sig_brown_merged %>% mutate(BrownPvalue_Genotype = as.numeric(as.character(BrownPvalue_Genotype)), 
                                                           BrownPvalue_Pathology = as.numeric(as.character(BrownPvalue_Pathology))) %>% 
    arrange(BrownPvalue_Genotype)
  
  return(dat_sig_brown_merged )
}

#-------------- apply function -------------

rTg4510_sig_brown <- merge_brownAggregation(sigBeta$rTg4510, sigRes$rTg4510)
J20_sig_brown <- merge_brownAggregation(sigBeta$J20, sigRes$J20)

#-------------- output -------------

write.csv(rTg4510_sig_brown, file = paste0(dirnames$annotated,"/final/rTg4510_sigBrown.csv"), quote = F, row.names = F)
write.csv(J20_sig_brown, file = paste0(dirnames$annotated,"/final/J20_sigBrown.csv"), quote = F, row.names = F)
