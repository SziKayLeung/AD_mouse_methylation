#Emma Walker
#17/05/21
#E.M.Walker@exeter.ac.uk
# adapted from Aisha's script
# using the beta regression as mixed models effect
# interactionModel <- methylation ~ Genotype + Age_months + Genotype*Age_months + Chip_ID
# Entorhinal Cortex only

library(lmerTest)
library(lmtest)
library(ggfortify)
library(glmmTMB)

setwd("/gpfs/ts0/projects/Research_Project-191406/EmmaW/")


######################### Set up pheno and betas files ready for analysis =============
#######################################################################################
#######################################################################################

# load in normalised QC'd data
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat")
QCmetrics$SampleID <- gsub(".*_", "", QCmetrics$ExternalSampleID)

J20_path <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_coldata_VertebrateArray.csv",
                     header = T, stringsAsFactors = F)
## ! The same samples do not have the same sample ID, these are matched on the pathology file. If we get the 
## same order of the sample ID as the pathology file then they are matched correctly
### ECX
pheno_J20_ECX <- QCmetrics[which(QCmetrics$Tissue == "CortexEntorhinalis" & QCmetrics$AD_model == "J20"),]
betas_J20_ECX <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno_J20_ECX$Basename ]
identical(pheno_J20_ECX$Basename, colnames(betas_J20_ECX))
colnames(betas_J20_ECX) <- pheno_J20_ECX$SampleID
rownames(pheno_J20_ECX) <- pheno_J20_ECX$SampleID 
pheno_J20_ECX <- pheno_J20_ECX[J20_path$Sample_ID_ECX, ] #order based on pathology samples id (NA rows will be generated)
pheno_J20_ECX <- cbind(pheno_J20_ECX, J20_path)
### Hip
pheno_J20_HIP <- QCmetrics[which(QCmetrics$Tissue == "Hippocampus" & QCmetrics$AD_model == "J20"),]
betas_J20_HIP <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno_J20_HIP$Basename ]
identical(pheno_J20_HIP$Basename, colnames(betas_J20_HIP))
colnames(betas_J20_HIP) <- pheno_J20_HIP$SampleID
rownames(pheno_J20_HIP) <- pheno_J20_HIP$SampleID 
pheno_J20_HIP <- pheno_J20_HIP[J20_path$Sample_ID_HIP, ] #order based on pathology samples id (NA rows will be generated)
pheno_J20_HIP <- cbind(pheno_J20_HIP, J20_path)

# Match the same mouse with 2 brain regions
samples <- as.data.frame(cbind(pheno_J20_ECX$SampleID, pheno_J20_HIP$SampleID))
colnames(samples) <- c("ECX", "HIP")
samples$ECX <- as.character(samples$ECX)
samples$HIP <- as.character(samples$HIP)
for(i in 1:nrow(samples)){
  samples[i, "MouseID"] <- paste("M", i, sep="")
}
## 36 matched samples!!

# Now filter these matching samples on both tissues
### ecx
identical(pheno_J20_ECX$SampleID, samples$ECX)
pheno_J20_ECX$MouseID <- samples$MouseID[match(samples$ECX, pheno_J20_ECX$SampleID)] #Now the orders are the same - add in the mouseIDs
identical(colnames(betas_J20_ECX) , pheno_J20_ECX$SampleID)
pheno_J20_ECX <- pheno_J20_ECX[colnames(betas_J20_ECX), ]
identical(colnames(betas_J20_ECX) , pheno_J20_ECX$SampleID)


### hip
identical(pheno_J20_HIP$SampleID, samples$HIP)
pheno_J20_HIP$MouseID <- samples$MouseID[match(samples$HIP, pheno_J20_HIP$SampleID)] #Now the orders are the same - add in the mouseIDs
identical(colnames(betas_J20_HIP) , pheno_J20_HIP$SampleID)
pheno_J20_HIP <- pheno_J20_HIP[colnames(betas_J20_HIP), ]
identical(colnames(betas_J20_HIP) , pheno_J20_HIP$SampleID)

#combine brain regions into one file
pheno <- rbind(pheno_J20_ECX, pheno_J20_HIP)
betas <- cbind(betas_J20_ECX, betas_J20_HIP)
identical(pheno$SampleID, colnames(betas))

# make sure model terms are the right data type
pheno$Genotype <- factor(pheno$Genotype, levels = c("WT","TG"))
pheno$Tissue <- factor(pheno$Tissue, levels = c("CortexEntorhinalis","Hippocampus"))
pheno$Chip_ID <- as.factor(pheno$Chip_ID)
pheno$MouseID <- as.factor(pheno$MouseID)

#subset pheno file to only contain model terms
pheno <- pheno[,c("Tissue","Genotype", "Age_months", "Chip_ID", "PlateNumber", "MouseID","Pathology_ECX", "Pathology_HIP")]
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


############################### Run ECX model =======================================
####################################################################################
####################################################################################

#subset to just entorhinal cortex samples
phenoJ20_ECX <- pheno[which(pheno$Tissue == "CortexEntorhinalis"),]
phenoJ20_ECX$Chip_ID <- factor(phenoJ20_ECX$Chip_ID)
betasJ20_ECX <- betas[,colnames(betas) %in% rownames(phenoJ20_ECX)]
identical(rownames(phenoJ20_ECX), colnames(betasJ20_ECX))

print("---------------- Starting Running Parallel --------------")

#define function containing model to be used in parallel mode
testCpG<-function(row, phenoJ20_ECX){
  
  library(lmtest)
  
  modelglmm<-glmmTMB(betasJ20_ECX[i,] ~ Genotype + Age_months + Genotype*Age_months + Chip_ID, data=phenoJ20_ECX, REML = FALSE, family=beta_family(link = "probit"))

  # extract info for model terms of interest
  Betas.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",1]
  SE.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",2]
  Z.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",3]
  PrZ.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",4]
  
  Betas.GenotypeTG <-summary(modelglmm)$coefficients$cond["GenotypeTG",1]
  SE.GenotypeTG <-summary(modelglmm)$coefficients$cond["GenotypeTG",2]
  Z.GenotypeTG <-summary(modelglmm)$coefficients$cond["GenotypeTG",3]
  PrZ.GenotypeTG <- summary(modelglmm)$coefficients$cond["GenotypeTG",4]
  
  Betas.Age_months <-summary(modelglmm)$coefficients$cond["Age_months",1]
  SE.Age_months <-summary(modelglmm)$coefficients$cond["Age_months",2]
  Z.Age_months <-summary(modelglmm)$coefficients$cond["Age_months",3]
  PrZ.Age_months <- summary(modelglmm)$coefficients$cond["Age_months",4]
  
  Betas.GenotypeTG.Age_months <-summary(modelglmm)$coefficients$cond["GenotypeTG:Age_months",1]
  SE.GenotypeTG.Age_months <-summary(modelglmm)$coefficients$cond["GenotypeTG:Age_months",2]
  Z.GenotypeTG.Age_months <-summary(modelglmm)$coefficients$cond["GenotypeTG:Age_months",3]
  PrZ.GenotypeTG.Age_months <- summary(modelglmm)$coefficients$cond["GenotypeTG:Age_months",4]
  
  ## extract info for chip in case useful later...
  Betas.Chip_ID204027420024 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420024",1]
  SE.Chip_ID204027420024 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420024",2]
  Z.Chip_ID204027420024 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420024",3]
  PrZ.Chip_ID204027420024 <- summary(modelglmm)$coefficients$cond["Chip_ID204027420024",4]
  
  Betas.Chip_ID204027420027 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420027",1]
  SE.Chip_ID204027420027 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420027",2]
  Z.Chip_ID204027420027 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420027",3]
  PrZ.Chip_ID204027420027 <- summary(modelglmm)$coefficients$cond["Chip_ID204027420027",4]
  
  Betas.Chip_ID204027420037 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420037",1]
  SE.Chip_ID204027420037 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420037",2]
  Z.Chip_ID204027420037 <-summary(modelglmm)$coefficients$cond["Chip_ID204027420037",3]
  PrZ.Chip_ID204027420037 <- summary(modelglmm)$coefficients$cond["Chip_ID204027420037",4] 
 
  return(c(Betas.Baseline, SE.Baseline, Z.Baseline, PrZ.Baseline,
           Betas.GenotypeTG, SE.GenotypeTG, Z.GenotypeTG, PrZ.GenotypeTG,
           Betas.Age_months, SE.Age_months, Z.Age_months, PrZ.Age_months,
           Betas.GenotypeTG.Age_months, SE.GenotypeTG.Age_months, Z.GenotypeTG.Age_months, PrZ.GenotypeTG.Age_months,
           Betas.Chip_ID204027420024, SE.Chip_ID204027420024, Z.Chip_ID204027420024, PrZ.Chip_ID204027420024,
           Betas.Chip_ID204027420027, SE.Chip_ID204027420027, Z.Chip_ID204027420027, PrZ.Chip_ID204027420027,
           Betas.Chip_ID204027420037, SE.Chip_ID204027420037, Z.Chip_ID204027420037, PrZ.Chip_ID204027420037))
}

print("---------------- Start Running Parallel --------------")

###Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel
res<-foreach(i=1:nrow(betasJ20_ECX), .combine=rbind) %dopar%{
  testCpG(betasJ20_ECX[i,], phenoJ20_ECX)
}

#write results to csv
write.csv(res,"/gpfs/ts0/projects/Research_Project-191406/EmmaW/Array/Results/rTg4510/MixedModelResultsECX_J20.csv")

print("---------------- Finished Running Parallel --------------")

## add colnames
columnnames <- c("Betas.Baseline", "SE.Baseline", "Z.Baseline", "PrZ.Baseline",
                 "Betas.GenotypeTG", "SE.GenotypeTG", "Z.GenotypeTG", "PrZ.GenotypeTG",
                 "Betas.Age_months", "SE.Age_months", "Z.Age_months", "PrZ.Age_months",
                 "Betas.GenotypeTG.Age_months", "SE.GenotypeTG.Age_months", "Z.GenotypeTG.Age_months", "PrZ.GenotypeTG.Age_months",
                 "Betas.Chip_ID204027420024", "SE.Chip_ID204027420024", "Z.Chip_ID204027420024", "PrZ.Chip_ID204027420024",
                 "Betas.Chip_ID204027420027", "SE.Chip_ID204027420027", "Z.Chip_ID204027420027", "PrZ.Chip_ID204027420027",
                 "Betas.Chip_ID204027420037", "SE.Chip_ID204027420037", "Z.Chip_ID204027420037", "PrZ.Chip_ID204027420037")

colnames(res)<-columnnames
rownames(res)<-rownames(betasJ20_ECX)
res<-as.data.frame(res)
head(res)

write.csv(res,"/gpfs/ts0/projects/Research_Project-191406/EmmaW/Array/Results/rTg4510/MixedModelResultsECX_J20.csv")
