### For the mixed models to run, I had to remove chip as that correlated with the tissues - so does removing chip influence our previous model?

## libraries
library(ggplot2)
library(tidyverse)
library(tidyr)
library(qqman)
library(cgwtools)
library(BiSeq)
library(reshape2)
library(gridExtra)
library(kableExtra)

# Load the unique dataset
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat")
source("/gpfs/ts0/projects/Research_Project-191406/Aisha/script/rrbs-ad-mice/PlotFunctions.R")
source("/gpfs/ts0/projects/Research_Project-191406/Aisha/script/rrbs-ad-mice/Array/BetaRegressionArray.R")
# Filter for rTg samples
pheno <- QCmetrics[which(QCmetrics$Tissue == "CortexEntorhinalis" & QCmetrics$AD_model == "rTg4510"),]
betas <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno$Basename ]
pheno$AgeMonths = ifelse(pheno$AgeDays < 100, 2, ifelse(pheno$AgeDays < 150, 4, ifelse(pheno$AgeDays < 200, 6, 8)))
pheno$Genotype <- as.factor(pheno$Genotype)
pheno$Sample_ID_ECX <- gsub(".*_", "", pheno$ExternalSampleID)

pathology <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Tg4510_coldata_VertebrateArray.csv",
                      header = T, stringsAsFactors = F)
phenowPath <- merge(pheno,pathology, by = "Sample_ID_ECX")
colnames(phenowPath)[colnames(phenowPath) == 'Genotype.x'] <- 'Genotype'
pheno <- phenowPath

betas <- betas[,order(colnames(betas))]
pheno <- pheno[order(pheno$Basename),]
identical(colnames(betas), pheno$Basename)
pheno$Genotype <- factor(pheno$Genotype, levels = c("WT","TG"))


#### parallel processors
# betas <- as.matrix(betas)
# library(doParallel)
# cl<-makeCluster(16)
# registerDoParallel(cl)
# clusterEvalQ(cl,library(glmmTMB))
# 
# #betas.org <- betas
# print("---------------- Starting Running Parallel --------------")
# ###Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel
# 
# min.meth <- min(betas[betas > 0], na.rm=TRUE)
# max.meth <- max(betas[betas < 1], na.rm=TRUE)
# 
# res<-foreach(j=1:3, .combine=rbind) %dopar%{
#   testCpG(betas[j,], pheno,
#           max.meth = max.meth, min.meth = min.meth,
#           link = "probit",
#           formula = ~Genotype + Age_months+ Genotype*Age_months,
#           formulaNull = ~Genotype + Age_months)
# }


betaResultsECXrTg <- BetaRegGenotypeNoChip(betas = betas,
                           pheno = pheno,
                           formula = ~Genotype + Age_months+ Genotype*Age_months,
                           formulaNull = ~Genotype + Age_months,
                           link = "probit")




### Hip rTg4510
pheno <- QCmetrics[which(QCmetrics$Tissue == "Hippocampus" & QCmetrics$AD_model == "rTg4510"),]
betas <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno$Basename ]
pheno$AgeMonths = ifelse(pheno$AgeDays < 100, 2, ifelse(pheno$AgeDays < 150, 4, ifelse(pheno$AgeDays < 200, 6, 8)))
pheno$Genotype <- as.factor(pheno$Genotype)
pathology_pheno <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Tg4510_coldata_VertebrateArray.csv", header = T, stringsAsFactors = F)
pheno$Sample_ID_HIP <- gsub(".*_", "", pheno$ExternalSampleID)
phenowPath <- merge(pheno,pathology_pheno, by = "Sample_ID_HIP")
colnames(phenowPath)[colnames(phenowPath) == 'Genotype.x'] <- 'Genotype'
pheno <- phenowPath
betas <- betas[,order(colnames(betas))]
pheno <- pheno[order(pheno$Basename),]
identical(colnames(betas), pheno$Basename)
pheno$Genotype <- factor(pheno$Genotype, levels = c("WT","TG"))

betaResultsHIPrTg <- BetaRegGenotypeNoChip(betas = betas,
                                          pheno = pheno,
                                          formula = ~Genotype + Age_months+ Genotype*Age_months,
                                          formulaNull = ~Genotype + Age_months,
                                          link = "probit")


### ECX J20
pheno <- QCmetrics[which(QCmetrics$Tissue == "CortexEntorhinalis" & QCmetrics$AD_model == "J20"),]
betas <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno$Basename ]
pheno$Genotype <- as.factor(pheno$Genotype)
pathology_pheno <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_ECX_Pathology.csv", header = T, stringsAsFactors = F)
colnames(pathology_pheno)[colnames(pathology_pheno) == 'X...Sample_ID_ECX'] <- 'Sample_ID_ECX'
pheno$Sample_ID_ECX <- gsub(".*_", "", pheno$ExternalSampleID)
phenowPath <- merge(pheno,pathology_pheno, by = "Sample_ID_ECX")
colnames(phenowPath)[colnames(phenowPath) == 'Genotype.x'] <- 'Genotype'
pheno <- phenowPath
betas <- betas[,order(colnames(betas))]
pheno <- pheno[order(pheno$Basename),]
identical(colnames(betas), pheno$Basename)
pheno$Genotype <- factor(pheno$Genotype, levels = c("WT","TG"))

betaResultsECXJ20 <- BetaRegGenotypeNoChip(betas = betas,
                                         pheno = pheno,
                                         formula = ~Genotype + Age_months+ Genotype*Age_months,
                                         formulaNull = ~Genotype + Age_months,
                                         link = "probit")

### Hip J20
pheno <- QCmetrics[which(QCmetrics$Tissue == "Hippocampus" & QCmetrics$AD_model == "J20"),]
betas <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno$Basename ]
pheno$Genotype <- as.factor(pheno$Genotype)
pathology_pheno <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_coldata_VertebrateArray.csv", header = T, stringsAsFactors = F)
pheno$Sample_ID_HIP <- gsub(".*_", "", pheno$ExternalSampleID) 
phenowPath <- merge(pheno,pathology_pheno, by = "Sample_ID_HIP")
colnames(phenowPath)[colnames(phenowPath) == 'Genotype.x'] <- 'Genotype'
pheno <- phenowPath
betas <- betas[,order(colnames(betas))]
pheno <- pheno[order(pheno$Basename),]
identical(colnames(betas), pheno$Basename)
pheno$Genotype <- factor(pheno$Genotype, levels = c("WT","TG"))

betaResultsHIPJ20 <- BetaRegGenotypeNoChip(betas = betas,
                                         pheno = pheno,
                                         formula = ~Genotype + Age_months+ Genotype*Age_months,
                                         formulaNull = ~Genotype + Age_months,
                                         link = "probit")




save(betaResultsECXrTg, betaResultsHIPrTg, betaResultsECXJ20, betaResultsHIPJ20,
     file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegwithoutChip.RData")



