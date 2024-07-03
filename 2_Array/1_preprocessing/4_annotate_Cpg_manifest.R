#!/usr/bin/env Rscript
## ----------Script-----------------  
##
## Purpose: Annotate cpG probe ID from array beta with horvath manifest file
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------
## 
## Using A.Dahir normalised file for array data (rTg4510 J20 ECX)
## annotate CpG with manifest file for co-ordinate position 
## subset rTg4510 and J20 ECX beta files
## ---------------------------------

library("dplyr")
library("stringr")

rootDir = "/lustre/projects/Research_Project-191406"
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config"))

# load Aisha's normalized dataset
load(paste0(dirnames$processed, "/array/Normalised_Data_Sesame.rdat"))

# datawrangle for ECX J20 and rTg4510
pheno <- QCmetrics
betas <- Normalised_Sesame_Betas

# ECX
pheno_ecx <- pheno[which(pheno$Tissue %in% "CortexEntorhinalis"),] %>% mutate(sample_ID = word(ExternalSampleID, c(2),sep=fixed("_")))
betas_ecx <- betas[, colnames(betas) %in% pheno_ecx$Basename] 

# rTg4510
pheno_ecx_Tg4510 <- pheno_ecx[which(pheno_ecx$AD_model %in% "rTg4510"),] 
betas_ecx_rTg4510  <- betas_ecx[, colnames(betas_ecx) %in% pheno_ecx_Tg4510$Basename] 
colnames(betas_ecx_rTg4510) <- unlist(lapply(colnames(betas_ecx_rTg4510), function(x) pheno_ecx[pheno_ecx$Basename == x, "sample_ID"]))
betas_ecx_rTg4510 <- merge(betas_ecx_rTg4510, mm10_Manifest, by = 0, all.x = T)
colnames(betas_ecx_rTg4510)[1] <- "cpg"
save(betas_ecx_rTg4510, file = paste0(dirnames$processed, "/array/Normalised_rTg4510_array_ECX.RData"))

# J20 
pheno_ecx_J20 <- pheno_ecx[which(pheno_ecx$AD_model %in% "J20"),] 
betas_ecx_J20 <- betas_ecx[, colnames(betas_ecx) %in% pheno_ecx_J20$Basename] 
colnames(betas_ecx_J20) <- unlist(lapply(colnames(betas_ecx_J20), function(x) pheno_ecx[pheno_ecx$Basename == x, "sample_ID"]))
betas_ecx_J20 <- merge(betas_ecx_J20, mm10_Manifest, by = 0, all.x = T)
colnames(betas_ecx_J20)[1] <- "cpg"
save(betas_ecx_J20, file = paste0(dirnames$processed, "/array/Normalised_J20_array_ECX.RData"))
