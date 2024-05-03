## 29.07.20
## Clean the two datasets of brain regions. Find the matching samples then run both regions together 
## using the beta regression as mixed models effect


# Run two brain regions together 
library(lmerTest)
library(lmtest)
library(ggfortify)
library(glmmTMB)

load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat")
QCmetrics$SampleID <- gsub(".*_", "", QCmetrics$ExternalSampleID)

###################################### rTg4510  #######################################
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
# samples <- samples[!is.na(samples$ECX),]
# samples <- samples[!is.na(samples$HIP),]
samples$ECX <- as.character(samples$ECX)
samples$HIP <- as.character(samples$HIP)
for(i in 1:nrow(samples)){
  samples[i, "MouseID"] <- paste("M", i, sep="")
}
## 36 matched samples!!

# Now filter these matching samples on both tissues
### ecx
#pheno_rTg4510_ECX <- pheno_rTg4510_ECX[samples$ECX,]
identical(pheno_rTg4510_ECX$SampleID, samples$ECX)
pheno_rTg4510_ECX$MouseID <- samples$MouseID[match(samples$ECX, pheno_rTg4510_ECX$SampleID)] #Now the orders are the same - add in the mouseIDs
#betas_rTg4510_ECX <- betas_rTg4510_ECX[,pheno_rTg4510_ECX$SampleID]
identical(colnames(betas_rTg4510_ECX) , pheno_rTg4510_ECX$SampleID)
pheno_rTg4510_ECX <- pheno_rTg4510_ECX[colnames(betas_rTg4510_ECX), ]
identical(colnames(betas_rTg4510_ECX) , pheno_rTg4510_ECX$SampleID)


### hip
#pheno_rTg4510_HIP <- pheno_rTg4510_HIP[samples$HIP,]
identical(pheno_rTg4510_HIP$SampleID, samples$HIP)
pheno_rTg4510_HIP$MouseID <- samples$MouseID[match(samples$HIP, pheno_rTg4510_HIP$SampleID)] #Now the orders are the same - add in the mouseIDs
#betas_rTg4510_HIP <- betas_rTg4510_HIP[,pheno_rTg4510_HIP$SampleID]
identical(colnames(betas_rTg4510_HIP) , pheno_rTg4510_HIP$SampleID)
pheno_rTg4510_HIP <- pheno_rTg4510_HIP[colnames(betas_rTg4510_HIP), ]
identical(colnames(betas_rTg4510_HIP) , pheno_rTg4510_HIP$SampleID)


pheno <- rbind(pheno_rTg4510_ECX, pheno_rTg4510_HIP)
betas <- cbind(betas_rTg4510_ECX, betas_rTg4510_HIP)
identical(pheno$SampleID, colnames(betas))


pheno$Genotype <- factor(pheno$Genotype, levels = c("WT","TG"))
pheno$Tissue <- factor(pheno$Tissue, levels = c("CortexEntorhinalis","Hippocampus"))
pheno$Chip_ID <- as.factor(pheno$Chip_ID)
pheno$MouseID <- as.factor(pheno$MouseID)

pheno <- pheno[,c("Tissue","Genotype", "Age_months", "Chip_ID", "PlateNumber", "MouseID", "Pathology_ECX", "Pathology_HIP")]
betas <- as.matrix(betas)

#combine pathology into one column
pheno$Pathology <- rep(NA, nrow(pheno))
for(i in 1:nrow(pheno)){
  pheno[i,"Pathology"] <- ifelse(pheno[i,"Tissue"] == 'CortexEntorhinalis', pheno[i,"Pathology_ECX"], 
                         pheno[i,"Pathology_HIP"])
}
pheno <- pheno[,c("Tissue","Genotype", "Age_months", "Chip_ID", "PlateNumber", "MouseID", "Pathology")]

## Plot PCA's for the  two brain regions 
# PCAs
# betas <- as.data.frame(betas)
# pheno <- as.data.frame(pheno)
# pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg4510_PCA.pdf")
# autoplot(prcomp(t(betas)), data = pheno, colour = 'Genotype', shape = "Tissue")
# dev.off()

#### parallel processors
betas <- as.matrix(betas)
library(doParallel)
cl<-makeCluster(16)
registerDoParallel(cl)
clusterEvalQ(cl,library(glmmTMB))

#betas.org <- betas
print("---------------- Starting Running Parallel --------------")
# #Create function which performs analysis for a probes
# testCpG<-function(row, pheno){
#   library(lmtest)
#   modelglmm<-glmmTMB(betas[i,] ~ Tissue + Tissue*Genotype + Genotype + Age_months + Genotype*Age_months + (1|MouseID),data=pheno, REML = FALSE,
#                      family=beta_family(link = "probit"))
#   Betas.Tissue <-summary(modelglmm)$coefficients$cond["TissueHippocampus",1]
#   SE.Tissue <-summary(modelglmm)$coefficients$cond["TissueHippocampus",2]
#   Z.Tissue <-summary(modelglmm)$coefficients$cond["TissueHippocampus",3]
#   PrZ.Tissue <- summary(modelglmm)$coefficients$cond["TissueHippocampus",4]
# 
#   Betas.TissueHippocampus.GenotypeTG <-summary(modelglmm)$coefficients$cond["TissueHippocampus:GenotypeTG",1]
#   SE.TissueHippocampus.GenotypeTG <-summary(modelglmm)$coefficients$cond["TissueHippocampus:GenotypeTG",2]
#   Z.TissueHippocampus.GenotypeTG <-summary(modelglmm)$coefficients$cond["TissueHippocampus:GenotypeTG",3]
#   PrZ.TissueHippocampus.GenotypeTG <- summary(modelglmm)$coefficients$cond["TissueHippocampus:GenotypeTG",4]
# 
#   Betas.GenotypeTG <-summary(modelglmm)$coefficients$cond["GenotypeTG",1]
#   SE.GenotypeTG <-summary(modelglmm)$coefficients$cond["GenotypeTG",2]
#   Z.GenotypeTG <-summary(modelglmm)$coefficients$cond["GenotypeTG",3]
#   PrZ.GenotypeTG <- summary(modelglmm)$coefficients$cond["GenotypeTG",4]
# 
#   Betas.GenotypeTG.Age_months <-summary(modelglmm)$coefficients$cond["GenotypeTG:Age_months",1]
#   SE.GenotypeTG.Age_months <-summary(modelglmm)$coefficients$cond["GenotypeTG:Age_months",2]
#   Z.GenotypeTG.Age_months <-summary(modelglmm)$coefficients$cond["GenotypeTG:Age_months",3]
#   PrZ.GenotypeTG.Age_months <- summary(modelglmm)$coefficients$cond["GenotypeTG:Age_months",4]
# 
# 
#   Betas.Age_months <-summary(modelglmm)$coefficients$cond["Age_months",1]
#   SE.Age_months <-summary(modelglmm)$coefficients$cond["Age_months",2]
#   Z.Age_months <-summary(modelglmm)$coefficients$cond["Age_months",3]
#   PrZ.Age_months <- summary(modelglmm)$coefficients$cond["Age_months",4]
# 
#   model.null<-glmmTMB(betas[i,] ~ Genotype + Age_months + Genotype*Age_months + (1|MouseID), data=pheno, REML = FALSE,
#                       family=beta_family(link = "probit"))
#   P<-lrtest(modelglmm,model.null)
#   P <-P$`Pr(>Chisq)`[2]
#   return(c(Betas.Tissue, SE.Tissue, Z.Tissue, PrZ.Tissue,
#            Betas.TissueHippocampus.GenotypeTG, SE.TissueHippocampus.GenotypeTG, Z.TissueHippocampus.GenotypeTG, PrZ.TissueHippocampus.GenotypeTG,
#            Betas.GenotypeTG, SE.GenotypeTG, Z.GenotypeTG, PrZ.GenotypeTG,
#            Betas.GenotypeTG.Age_months, SE.GenotypeTG.Age_months, Z.GenotypeTG.Age_months, PrZ.GenotypeTG.Age_months,
#            Betas.Age_months, SE.Age_months, Z.Age_months, PrZ.Age_months ,
#            P))
# }
# ###Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel
# res<-foreach(i=1:nrow(betas), .combine=rbind) %dopar%{
#   testCpG(betas[i,], pheno)
# }
# 
# print("---------------- Finished Running Parallel --------------")
# ## add colnames
# 
# columnnames <- c("Betas.Tissue", "SE.Tissue", "Z.Tissue", "PrZ.Tissue",
#                  "Betas.TissueHippocampus.GenotypeTG", "SE.TissueHippocampus.GenotypeTG", "Z.TissueHippocampus.GenotypeTG", "PrZ.TissueHippocampus.GenotypeTG",
#                  "Betas.GenotypeTG", "SE.GenotypeTG", "Z.GenotypeTG", "PrZ.GenotypeTG",
#                  "Betas.GenotypeTG.Age_months", "SE.GenotypeTG.Age_months", "Z.GenotypeTG.Age_months", "PrZ.GenotypeTG.Age_months",
#                  "Betas.Age_months", "SE.Age_months", "Z.Age_months", "PrZ.Age_months",
#                  "lrt.Pval")
# 
# colnames(res)<-columnnames
# rownames(res)<-rownames(betas)
# res<-as.data.frame(res)
# head(res)
# 
# 
# write.csv(res,"/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/MixedModelResults_rTg4510.csv")


############################### Another model with pathology

#Create function which performs analysis for a probes
testCpG<-function(row, pheno){
  library(lmtest)
  modelglmm<-glmmTMB(betas[i,] ~ Tissue + Tissue*Genotype + Genotype + Age_months + Pathology + Tissue*Pathology + (1|MouseID),data=pheno, REML = FALSE,
                     family=beta_family(link = "probit"))
  
  Betas.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",1]
  SE.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",2]
  Z.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",3]
  PrZ.Baseline <- summary(modelglmm)$coefficients$cond["(Intercept)",4]
  
  Betas.Tissue <-summary(modelglmm)$coefficients$cond["TissueHippocampus",1]
  SE.Tissue <-summary(modelglmm)$coefficients$cond["TissueHippocampus",2]
  Z.Tissue <-summary(modelglmm)$coefficients$cond["TissueHippocampus",3]
  PrZ.Tissue <- summary(modelglmm)$coefficients$cond["TissueHippocampus",4]

  Betas.TissueHippocampus.GenotypeTG <-summary(modelglmm)$coefficients$cond["TissueHippocampus:GenotypeTG",1]
  SE.TissueHippocampus.GenotypeTG <-summary(modelglmm)$coefficients$cond["TissueHippocampus:GenotypeTG",2]
  Z.TissueHippocampus.GenotypeTG <-summary(modelglmm)$coefficients$cond["TissueHippocampus:GenotypeTG",3]
  PrZ.TissueHippocampus.GenotypeTG <- summary(modelglmm)$coefficients$cond["TissueHippocampus:GenotypeTG",4]

  Betas.GenotypeTG <-summary(modelglmm)$coefficients$cond["GenotypeTG",1]
  SE.GenotypeTG <-summary(modelglmm)$coefficients$cond["GenotypeTG",2]
  Z.GenotypeTG <-summary(modelglmm)$coefficients$cond["GenotypeTG",3]
  PrZ.GenotypeTG <- summary(modelglmm)$coefficients$cond["GenotypeTG",4]

  Betas.Age_months <-summary(modelglmm)$coefficients$cond["Age_months",1]
  SE.Age_months <-summary(modelglmm)$coefficients$cond["Age_months",2]
  Z.Age_months <-summary(modelglmm)$coefficients$cond["Age_months",3]
  PrZ.Age_months <- summary(modelglmm)$coefficients$cond["Age_months",4]

  Betas.Pathology <-summary(modelglmm)$coefficients$cond["Pathology",1]
  SE.Pathology <-summary(modelglmm)$coefficients$cond["Pathology",2]
  Z.Pathology <-summary(modelglmm)$coefficients$cond["Pathology",3]
  PrZ.Pathology <- summary(modelglmm)$coefficients$cond["Pathology",4]
  
  Betas.TissueHippocampus.Pathology <-summary(modelglmm)$coefficients$cond["TissueHippocampus:Pathology",1]
  SE.TissueHippocampus.Pathology <-summary(modelglmm)$coefficients$cond["TissueHippocampus:Pathology",2]
  Z.TissueHippocampus.Pathology <-summary(modelglmm)$coefficients$cond["TissueHippocampus:Pathology",3]
  PrZ.TissueHippocampus.Pathology <- summary(modelglmm)$coefficients$cond["TissueHippocampus:Pathology",4]

  model.null<-glmmTMB(betas[i,] ~ Genotype + Age_months + Pathology + (1|MouseID), data=pheno, REML = FALSE,
                      family=beta_family(link = "probit"))
  P<-lrtest(modelglmm,model.null)
  P <-P$`Pr(>Chisq)`[2]
  return(c(Betas.Baseline, SE.Baseline, Z.Baseline, PrZ.Baseline,
           Betas.Tissue, SE.Tissue, Z.Tissue, PrZ.Tissue,
           Betas.TissueHippocampus.GenotypeTG, SE.TissueHippocampus.GenotypeTG, Z.TissueHippocampus.GenotypeTG, PrZ.TissueHippocampus.GenotypeTG,
           Betas.GenotypeTG, SE.GenotypeTG, Z.GenotypeTG, PrZ.GenotypeTG,
           Betas.Age_months, SE.Age_months, Z.Age_months, PrZ.Age_months ,
           Betas.Pathology, SE.Pathology, Z.Pathology, PrZ.Pathology,
           Betas.TissueHippocampus.Pathology, SE.TissueHippocampus.Pathology, Z.TissueHippocampus.Pathology, PrZ.TissueHippocampus.Pathology,
           P))
}
###Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel
res<-foreach(i=1:nrow(betas), .combine=rbind) %dopar%{
  testCpG(betas[i,], pheno)
}

print("---------------- Finished Running Parallel --------------")
## add colnames

columnnames <- c("Betas.Baseline", "SE.Baseline", "Z.Baseline", "PrZ.Baseline",
                 "Betas.Tissue", "SE.Tissue", "Z.Tissue", "PrZ.Tissue",
                 "Betas.TissueHippocampus.GenotypeTG", "SE.TissueHippocampus.GenotypeTG", "Z.TissueHippocampus.GenotypeTG", "PrZ.TissueHippocampus.GenotypeTG",
                 "Betas.GenotypeTG", "SE.GenotypeTG", "Z.GenotypeTG", "PrZ.GenotypeTG",
                 "Betas.Age_months", "SE.Age_months", "Z.Age_months", "PrZ.Age_months",
                 "Betas.Pathology", "SE.Pathology", "Z.Pathology", "PrZ.Pathology",
                 "Betas.TissueHippocampus.Pathology", "SE.TissueHippocampus.Pathology", "Z.TissueHippocampus.Pathology", "PrZ.TissueHippocampus.Pathology",
                 "lrt.Pval")

colnames(res)<-columnnames
rownames(res)<-rownames(betas)
res<-as.data.frame(res)
head(res)


write.csv(res,"/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/MixedModelResultswPathology_rTg4510.csv")


