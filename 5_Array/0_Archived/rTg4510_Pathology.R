## Need to run three pathology models. Can submit as job on isca
## rTg ECX & HIP Pathology


####################################  		 Entorhinal Cortex 			##############################################
library(ggplot2)
library(tidyverse)
library(tidyr)
library(qqman)
library(cgwtools)
library(BiSeq)
library(reshape2)
library(gridExtra)
library(kableExtra)
## Load data
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat")
source("/gpfs/ts0/projects/Research_Project-191406/Aisha/script/rrbs-ad-mice/PlotFunctions.R")
source("/gpfs/ts0/projects/Research_Project-191406/Aisha/script/rrbs-ad-mice/BetaRegressionArray.R")
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

########## Pathology + Chip
ResultsPathology <- BetaRegressionArrayPathology(betas = betas,
                                                     pheno = pheno,
                                                     formula = ~Pathology_ECX + Chip_ID,
                                                     link = "probit")
save(ResultsPathology, file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_ECX_rTg4510.RData")

########## Pathology + Age + Chip
ResultsPathologyAge <- BetaRegressionArrayPathology(betas = betas,
                                                     pheno = pheno,
                                                     formula = ~Pathology_ECX + Age_months + Chip_ID,
                                                     link = "probit")
resave(ResultsPathologyAge, file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_ECX_rTg4510.RData")

########## Pathology + Genotype + Chip
ResultsPathologyGenotype <- BetaRegressionArrayPathology(betas = betas,
                                                     pheno = pheno,
                                                     formula = ~Pathology_ECX + Genotype + Chip_ID,
                                                     link = "probit")
resave(ResultsPathologyGenotype, file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_ECX_rTg4510.RData")


########## Pathology + Age + Genotype + Chip
ResultsPathologyAgeGenotype <- BetaRegressionArrayPathology(betas = betas,
                                                     pheno = pheno,
                                                     formula = ~Pathology_ECX + Age_months + Genotype + Chip_ID,
                                                     link = "probit")
resave(ResultsPathologyAgeGenotype, file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_ECX_rTg4510.RData")




####################################  		 Hippocampus 			##############################################
rm(list=ls())
library(ggplot2)
library(tidyverse)
library(tidyr)
library(qqman)
library(cgwtools)
library(BiSeq)
library(reshape2)
library(gridExtra)
library(kableExtra)

load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat")
source("/gpfs/ts0/projects/Research_Project-191406/Aisha/script/rrbs-ad-mice/PlotFunctions.R")
source("/gpfs/ts0/projects/Research_Project-191406/Aisha/script/rrbs-ad-mice/BetaRegressionArray.R")
# Filter for rTg samples
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

########## Pathology + Chip
ResultsPathology <- BetaRegressionArrayPathology(betas = betas,
                                                     pheno = pheno,
                                                     formula = ~Pathology_HIP + Chip_ID,
                                                     link = "probit")
save(ResultsPathology, file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_HIP_rTg4510.RData")

########## Pathology + Age + Chip
ResultsPathologyAge <- BetaRegressionArrayPathology(betas = betas,
                                                     pheno = pheno,
                                                     formula = ~Pathology_HIP + Age_months + Chip_ID,
                                                     link = "probit")
resave(ResultsPathologyAge, file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_HIP_rTg4510.RData")


########## Pathology + Genotype + Chip
ResultsPathologyGenotype <- BetaRegressionArrayPathology(betas = betas,
                                                     pheno = pheno,
                                                     formula = ~Pathology_HIP + Genotype + Chip_ID,
                                                     link = "probit")
resave(ResultsPathologyGenotype, file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_HIP_rTg4510.RData")


########## Pathology + Age + Genotype + Chip
ResultsPathologyAgeGenotype <- BetaRegressionArrayPathology(betas = betas,
                                                     pheno = pheno,
                                                     formula = ~Pathology_HIP + Age_months + Genotype + Chip_ID,
                                                     link = "probit")
resave(ResultsPathologyAgeGenotype, file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_HIP_rTg4510.RData")
