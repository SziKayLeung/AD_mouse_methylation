#Emma Walker
#17/05/21
#E.M.Walker@exeter.ac.uk
# adapted from Aisha's script
# using the beta regression as mixed models effect
# PathologyModel <- methylation ~ Pathology + Chip_ID
# Entorhinal Cortex only

library(lmerTest)
library(lmtest)
library(ggfortify)
library(glmmTMB)

setwd("/gpfs/ts0/projects/Research_Project-191406/EmmaW/")

######################### Set up pheno and betas files ready for analysis =============
#######################################################################################
#######################################################################################

load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat")
QCmetrics$SampleID <- gsub(".*_", "", QCmetrics$ExternalSampleID)

rTg4510_path <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Tg4510_coldata_VertebrateArray.csv",
                         header = T, stringsAsFactors = F)

## ! The same samples do not have the same sample ID, these are matched on the pathology file. If we get the 
## same order of the sample ID as the pathology file then they are matched correctly
### ECX
pheno_rTg4510_ECX <- QCmetrics[which(QCmetrics$Tissue == "CortexEntorhinalis" & QCmetrics$AD_model == "rTg4510"),]
betas_rTg4510_ECX <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno_rTg4510_ECX$Basename ]
identical(pheno_rTg4510_ECX$Basename, colnames(betas_rTg4510_ECX))
colnames(betas_rTg4510_ECX) <- pheno_rTg4510_ECX$SampleID
rownames(pheno_rTg4510_ECX) <- pheno_rTg4510_ECX$SampleID 
pheno_rTg4510_ECX <- pheno_rTg4510_ECX[rTg4510_path$Sample_ID_ECX, ] #order based on pathology samples id (NA rows will be generated)
pheno_rTg4510_ECX <- cbind(pheno_rTg4510_ECX, rTg4510_path)
### Hip
pheno_rTg4510_HIP <- QCmetrics[which(QCmetrics$Tissue == "Hippocampus" & QCmetrics$AD_model == "rTg4510"),]
betas_rTg4510_HIP <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno_rTg4510_HIP$Basename ]
identical(pheno_rTg4510_HIP$Basename, colnames(betas_rTg4510_HIP))
colnames(betas_rTg4510_HIP) <- pheno_rTg4510_HIP$SampleID
rownames(pheno_rTg4510_HIP) <- pheno_rTg4510_HIP$SampleID 
pheno_rTg4510_HIP <- pheno_rTg4510_HIP[rTg4510_path$Sample_ID_HIP, ] #order based on pathology samples id (NA rows will be generated)
pheno_rTg4510_HIP <- cbind(pheno_rTg4510_HIP, rTg4510_path)

# Match the same mouse with 2 brain regions
samples <- as.data.frame(cbind(pheno_rTg4510_ECX$SampleID, pheno_rTg4510_HIP$SampleID))
colnames(samples) <- c("ECX", "HIP")
samples$ECX <- as.character(samples$ECX)
samples$HIP <- as.character(samples$HIP)
for(i in 1:nrow(samples)){
  samples[i, "MouseID"] <- paste("M", i, sep="")
}
## 36 matched samples!!

# Now filter these matching samples on both tissues
### ecx
identical(pheno_rTg4510_ECX$SampleID, samples$ECX)
pheno_rTg4510_ECX$MouseID <- samples$MouseID[match(samples$ECX, pheno_rTg4510_ECX$SampleID)] #Now the orders are the same - add in the mouseIDs
identical(colnames(betas_rTg4510_ECX) , pheno_rTg4510_ECX$SampleID)
pheno_rTg4510_ECX <- pheno_rTg4510_ECX[colnames(betas_rTg4510_ECX), ]
identical(colnames(betas_rTg4510_ECX) , pheno_rTg4510_ECX$SampleID)


### hip
identical(pheno_rTg4510_HIP$SampleID, samples$HIP)
pheno_rTg4510_HIP$MouseID <- samples$MouseID[match(samples$HIP, pheno_rTg4510_HIP$SampleID)] #Now the orders are the same - add in the mouseIDs
identical(colnames(betas_rTg4510_HIP) , pheno_rTg4510_HIP$SampleID)
pheno_rTg4510_HIP <- pheno_rTg4510_HIP[colnames(betas_rTg4510_HIP), ]
identical(colnames(betas_rTg4510_HIP) , pheno_rTg4510_HIP$SampleID)


#combine brain regions into one file
pheno <- rbind(pheno_rTg4510_ECX, pheno_rTg4510_HIP)
betas <- cbind(betas_rTg4510_ECX, betas_rTg4510_HIP)
identical(pheno$SampleID, colnames(betas))

# make sure model terms are the right data type
pheno$Genotype <- factor(pheno$Genotype, levels = c("WT","TG"))
pheno$Tissue <- factor(pheno$Tissue, levels = c("CortexEntorhinalis","Hippocampus"))
pheno$Chip_ID <- as.factor(pheno$Chip_ID)
pheno$MouseID <- as.factor(pheno$MouseID)

#subset pheno file to only contain model terms
pheno <- pheno[,c("Tissue","Genotype", "Age_months", "Chip_ID", "PlateNumber", "MouseID", "Pathology_ECX", "Pathology_HIP")]
betas <- as.matrix(betas)

#combine pathology into one column
pheno$Pathology <- rep(NA, nrow(pheno))
for(i in 1:nrow(pheno)){
  pheno[i,"Pathology"] <- ifelse(pheno[i,"Tissue"] == 'CortexEntorhinalis', pheno[i,"Pathology_ECX"], 
                         pheno[i,"Pathology_HIP"])
}
pheno <- pheno[,c("Tissue","Genotype", "Age_months", "Chip_ID", "PlateNumber", "MouseID", "Pathology")]


#### parallel processors
betas <- as.matrix(betas)
library(doParallel)
cl<-makeCluster(16)
registerDoParallel(cl)
clusterEvalQ(cl,library(glmmTMB))

############################### Run ECX Pathology model ============================
####################################################################################
####################################################################################

#subset to just entorhinal cortex samples
phenoTg_ECX <- pheno[which(pheno$Tissue == "CortexEntorhinalis"),]
phenoTg_ECX$Chip_ID <- factor(phenoTg_ECX$Chip_ID)
betasTg_ECX <- betas[,colnames(betas) %in% rownames(phenoTg_ECX)]
identical(rownames(phenoTg_ECX), colnames(betasTg_ECX))

print("---------------- Starting Running Parallel --------------")

#Create function which performs analysis for a probes
testCpG<-function(row, phenoTg_ECX){
  
  library(lmtest)
  
  modelglmm<-glmmTMB(betasTg_ECX[i,] ~ Pathology + Chip_ID, data=phenoTg_ECX, REML = FALSE, family=beta_family(link = "probit"))
  
  # extract info for model terms of interest
  Betas.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",1]
  SE.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",2]
  Z.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",3]
  PrZ.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",4]
  
  Betas.Pathology <-summary(modelglmm)$coefficients$cond["Pathology",1]
  SE.Pathology <-summary(modelglmm)$coefficients$cond["Pathology",2]
  Z.Pathology <-summary(modelglmm)$coefficients$cond["Pathology",3]
  PrZ.Pathology <- summary(modelglmm)$coefficients$cond["Pathology",4]
  
  ## extract info for chip in case useful later...
  Betas.Chip_ID204027420007 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420007",1]
  SE.Chip_ID204027420007 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420007",2]
  Z.Chip_ID204027420007 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420007",3]
  PrZ.Chip_ID204027420007 <- summary(modelglmm)$coefficients$cond["Chip_ID204027420007",4]
  
  Betas.Chip_ID204027420008 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420008",1]
  SE.Chip_ID204027420008 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420008",2]
  Z.Chip_ID204027420008 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420008",3]
  PrZ.Chip_ID204027420008 <- summary(modelglmm)$coefficients$cond["Chip_ID204027420008",4]
  
  Betas.Chip_ID204027420014 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420014",1]
  SE.Chip_ID204027420014 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420014",2]
  Z.Chip_ID204027420014 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420014",3]
  PrZ.Chip_ID204027420014 <- summary(modelglmm)$coefficients$cond["Chip_ID204027420014",4]  
  
  return(c(Betas.Baseline, SE.Baseline, Z.Baseline, PrZ.Baseline,
           Betas.Pathology, SE.Pathology, Z.Pathology, PrZ.Pathology,
           Betas.Chip_ID204027420007, SE.Chip_ID204027420007, Z.Chip_ID204027420007, PrZ.Chip_ID204027420007,
           Betas.Chip_ID204027420008, SE.Chip_ID204027420008, Z.Chip_ID204027420008, PrZ.Chip_ID204027420008,
           Betas.Chip_ID204027420014, SE.Chip_ID204027420014, Z.Chip_ID204027420014, PrZ.Chip_ID204027420014))
}


###Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel
res<-foreach(i=1:nrow(betasTg_ECX), .combine=rbind) %dopar%{
  testCpG(betasTg_ECX[i,], phenoTg_ECX)
}

#write results to csv
write.csv(res,"/gpfs/ts0/projects/Research_Project-191406/EmmaW/Array/Results/rTg4510/MixedModelResultsECX_Patholgy_rTg4510.csv")

print("---------------- Finished Running Parallel --------------")
## add colnames

columnnames <- c("Betas.Baseline", "SE.Baseline", "Z.Baseline", "PrZ.Baseline",
                 "Betas.Pathology", "SE.Pathology", "Z.Pathology", "PrZ.Pathology",
                 "Betas.Chip_ID204027420007", "SE.Chip_ID204027420007", "Z.Chip_ID204027420007", "PrZ.Chip_ID204027420007",
                 "Betas.Chip_ID204027420008", "SE.Chip_ID204027420008", "Z.Chip_ID204027420008", "PrZ.Chip_ID204027420008",
                 "Betas.Chip_ID204027420014", "SE.Chip_ID204027420014", "Z.Chip_ID204027420014", "PrZ.Chip_ID204027420014")

colnames(res)<-columnnames
rownames(res)<-rownames(betasTg_ECX)
res<-as.data.frame(res)
head(res)

#write results to csv
write.csv(res,"/gpfs/ts0/projects/Research_Project-191406/EmmaW/Array/Results/rTg4510/MixedModelResultsECX_Patholgy_rTg4510.csv")


