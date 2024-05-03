#### Does Aging in TG differ to WT

library(tidyr)
library(dplyr)
library(lme4)
#################### ECX rTg4510
## load data

## read data in
setwd("C:/Users/user/Desktop/Exe/vertarray/AD Project")
load("data/Array/Normalised_Data_Sesame.rdat") #load this in to have passed QC samples only as samplesheet has failed samples
pheno <- QCmetrics[which(QCmetrics$Tissue %in% c("CortexEntorhinalis", "Hippocampus")),] #take tissue samples further


samplesheet <- read.csv("OutputMouseClocks.csv", header = T,stringsAsFactors = F)
samplesheet <- samplesheet[,c(1:50)] #filter unnecessary values
samplesheet <- samplesheet[which(samplesheet$Basename %in% pheno$Basename),] #filter to just tissue samples to join with pathology data
samplesheet$SampleID <- gsub(".*_", "",samplesheet$ExternalSampleID)
## Add pathology data in
J20_pathdata <- read.csv("data/Array/J20_coldata_VertebrateArray.csv", header = T, stringsAsFactors = F)
rTg4510_pathdata <- read.csv("data/Array/Tg4510_coldata_VertebrateArray.csv", header = T, stringsAsFactors = F)

#Pathology knows which samples are paired so add the column for it now
for(i in 1:nrow(J20_pathdata)){
  J20_pathdata[i, "MouseID"] <- paste("J20_M", i, sep="")
}

for(i in 1:nrow(rTg4510_pathdata)){
  rTg4510_pathdata[i, "MouseID"] <- paste("rTg4510_M", i, sep="")
}

#split the tissue and model into sep data
J20_ecx <- J20_pathdata[,c("Sample_ID_ECX","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology_ECX", "MouseID" )]
colnames(J20_ecx) <- c("SampleID","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology", "MouseID")
J20_hip <- J20_pathdata[,c("Sample_ID_HIP","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology_HIP", "MouseID")]
colnames(J20_hip) <- c("SampleID","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology", "MouseID")

rTg4510_ecx <- rTg4510_pathdata[,c("Sample_ID_ECX","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology_ECX", "MouseID")]
colnames(rTg4510_ecx) <- c("SampleID","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology", "MouseID")
rTg4510_hip <- rTg4510_pathdata[,c("Sample_ID_HIP","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology_HIP", "MouseID")]
colnames(rTg4510_hip) <- c("SampleID","Mouse_ID","Genotype","Age_months","Group_ID","Histology_no","Pathology", "MouseID")

pathology <- rbind(J20_ecx, J20_hip, rTg4510_ecx, rTg4510_hip)
samplesheet <- merge(samplesheet, pathology, by = "SampleID")
samplesheet <- samplesheet[,-53]
colnames(samplesheet)[colnames(samplesheet) == 'Genotype.x'] <- 'Genotype'

samplesheet <- separate(data = samplesheet, col = Basename, 
                        into = c("Chip", "Chip_Position"), sep="_")
samplesheet$Genotype <- factor(samplesheet$Genotype, levels = c("WT","TG"))





################ DNAmAge ~ Genotype + Age + Genotype*Age,

#create a results table
res <- matrix(NA, nrow = 10, ncol= 4)
colnames(res) <- c("rTg4510_ECX", "rTg4510_HIP", "J20_ECX", "J20_HIP")
rownames(res) <- c("Intercept", "GenotypeTG.Beta", "GenotypeTG.SE", "GenotypeTG.P",
                   "Age.Beta","Age.SE", "Age.P",
                   "GenotypeTG.Age.Beta","GenotypeTG.Age.SE", "GenotypeTG.Age.P")

#rTg4510 ECX 
Tissue <- "CortexEntorhinalis" 
Model  <- "rTg4510"
dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]

lmodel <- lm(DNAmAgeClockCortex ~ Genotype + Age + Genotype*Age, data=dat)
summary(lmodel)$coef

res[1,1] <- summary(lmodel)$coef["(Intercept)",1]
res[2,1] <- summary(lmodel)$coef["GenotypeTG",1]
res[3,1] <- summary(lmodel)$coef["GenotypeTG",2]
res[4,1] <- summary(lmodel)$coef["GenotypeTG",4]


res[5,1] <- summary(lmodel)$coef["Age",1]
res[6,1] <- summary(lmodel)$coef["Age",2]
res[7,1] <- summary(lmodel)$coef["Age",4]


res[8,1] <- summary(lmodel)$coef["GenotypeTG:Age",1]
res[9,1] <- summary(lmodel)$coef["GenotypeTG:Age",2]
res[10,1] <- summary(lmodel)$coef["GenotypeTG:Age",4]

#rTg4510 HIP
Tissue <- "Hippocampus" 
Model  <- "rTg4510"
dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]

lmodel <- lm(DNAmAgeClockCortex ~ Genotype + Age + Genotype*Age, data=dat)
summary(lmodel)$coef

res[1,2] <- summary(lmodel)$coef["(Intercept)",1]
res[2,2] <- summary(lmodel)$coef["GenotypeTG",1]
res[3,2] <- summary(lmodel)$coef["GenotypeTG",2]
res[4,2] <- summary(lmodel)$coef["GenotypeTG",4]

res[5,2] <- summary(lmodel)$coef["Age",1]
res[6,2] <- summary(lmodel)$coef["Age",2]
res[7,2] <- summary(lmodel)$coef["Age",4]

res[8,2] <- summary(lmodel)$coef["GenotypeTG:Age",1]
res[9,2] <- summary(lmodel)$coef["GenotypeTG:Age",2]
res[10,2] <- summary(lmodel)$coef["GenotypeTG:Age",4]

#J20 ECX
Tissue <- "CortexEntorhinalis" 
Model  <- "J20"
dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]

lmodel <- lm(DNAmAgeClockCortex ~ Genotype + Age + Genotype*Age, data=dat)
summary(lmodel)$coef

res[1,3] <- summary(lmodel)$coef["(Intercept)",1]
res[2,3] <- summary(lmodel)$coef["GenotypeTG",1]
res[3,3] <- summary(lmodel)$coef["GenotypeTG",2]
res[4,3] <- summary(lmodel)$coef["GenotypeTG",4]

res[5,3] <- summary(lmodel)$coef["Age",1]
res[6,3] <- summary(lmodel)$coef["Age",2]
res[7,3] <- summary(lmodel)$coef["Age",4]

res[8,3] <- summary(lmodel)$coef["GenotypeTG:Age",1]
res[9,3] <- summary(lmodel)$coef["GenotypeTG:Age",2]
res[10,3] <- summary(lmodel)$coef["GenotypeTG:Age",4]


#J20 HIP
Tissue <- "Hippocampus" 
Model  <- "J20"
dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]

lmodel <- lm(DNAmAgeClockCortex ~ Genotype + Age + Genotype*Age, data=dat)
summary(lmodel)$coef

res[1,4] <- summary(lmodel)$coef["(Intercept)",1]
res[2,4] <- summary(lmodel)$coef["GenotypeTG",1]
res[3,4] <- summary(lmodel)$coef["GenotypeTG",2]
res[4,4] <- summary(lmodel)$coef["GenotypeTG",4]

res[5,4] <- summary(lmodel)$coef["Age",1]
res[6,4] <- summary(lmodel)$coef["Age",2]
res[7,4] <- summary(lmodel)$coef["Age",4]

res[8,4] <- summary(lmodel)$coef["GenotypeTG:Age",1]
res[9,4] <- summary(lmodel)$coef["GenotypeTG:Age",2]
res[10,4] <- summary(lmodel)$coef["GenotypeTG:Age",4]


write.csv(res, "DNAmAgeClockCortex_GenotypeAgeInteraction_results.csv")











################ DNAmAge ~ Genotype + Age + Genotype*Age,
#create a results table
res <- matrix(NA, nrow = 10, ncol= 4)
colnames(res) <- c("rTg4510_ECX", "rTg4510_HIP", "J20_ECX", "J20_HIP")
rownames(res) <- c("Intercept", "GenotypeTG.Beta", "GenotypeTG.SE", "GenotypeTG.P",
                   "Age.Beta","Age.SE", "Age.P",
                   "Pathology.Beta","Pathology.SE", "Pathology.P")

#rTg4510 ECX 
Tissue <- "CortexEntorhinalis" 
Model  <- "rTg4510"
dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]

lmodel <- lm(DNAmAgeClockCortex ~ Genotype + Age + Pathology, data=dat)
summary(lmodel)$coef

res[1,1] <- summary(lmodel)$coef["(Intercept)",1]
res[2,1] <- summary(lmodel)$coef["GenotypeTG",1]
res[3,1] <- summary(lmodel)$coef["GenotypeTG",2]
res[4,1] <- summary(lmodel)$coef["GenotypeTG",4]


res[5,1] <- summary(lmodel)$coef["Age",1]
res[6,1] <- summary(lmodel)$coef["Age",2]
res[7,1] <- summary(lmodel)$coef["Age",4]


res[8,1] <- summary(lmodel)$coef["Pathology",1]
res[9,1] <- summary(lmodel)$coef["Pathology",2]
res[10,1] <- summary(lmodel)$coef["Pathology",4]

#rTg4510 HIP
Tissue <- "Hippocampus" 
Model  <- "rTg4510"
dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]

lmodel <- lm(DNAmAgeClockCortex ~ Genotype + Age + Pathology, data=dat)
summary(lmodel)$coef

res[1,2] <- summary(lmodel)$coef["(Intercept)",1]
res[2,2] <- summary(lmodel)$coef["GenotypeTG",1]
res[3,2] <- summary(lmodel)$coef["GenotypeTG",2]
res[4,2] <- summary(lmodel)$coef["GenotypeTG",4]

res[5,2] <- summary(lmodel)$coef["Age",1]
res[6,2] <- summary(lmodel)$coef["Age",2]
res[7,2] <- summary(lmodel)$coef["Age",4]

res[8,2] <- summary(lmodel)$coef["Pathology",1]
res[9,2] <- summary(lmodel)$coef["Pathology",2]
res[10,2] <- summary(lmodel)$coef["Pathology",4]

#J20 ECX
Tissue <- "CortexEntorhinalis" 
Model  <- "J20"
dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]

lmodel <- lm(DNAmAgeClockCortex ~ Genotype + Age + Pathology, data=dat)
summary(lmodel)$coef

res[1,3] <- summary(lmodel)$coef["(Intercept)",1]
res[2,3] <- summary(lmodel)$coef["GenotypeTG",1]
res[3,3] <- summary(lmodel)$coef["GenotypeTG",2]
res[4,3] <- summary(lmodel)$coef["GenotypeTG",4]

res[5,3] <- summary(lmodel)$coef["Age",1]
res[6,3] <- summary(lmodel)$coef["Age",2]
res[7,3] <- summary(lmodel)$coef["Age",4]

res[8,3] <- summary(lmodel)$coef["Pathology",1]
res[9,3] <- summary(lmodel)$coef["Pathology",2]
res[10,3] <- summary(lmodel)$coef["Pathology",4]


#J20 HIP
Tissue <- "Hippocampus" 
Model  <- "J20"
dat <- samplesheet[which(samplesheet$Tissue == Tissue & samplesheet$AD_model == Model),]

lmodel <- lm(DNAmAgeClockCortex ~ Genotype + Age + Pathology, data=dat)
summary(lmodel)$coef


res[1,4] <- summary(lmodel)$coef["(Intercept)",1]
res[2,4] <- summary(lmodel)$coef["GenotypeTG",1]
res[3,4] <- summary(lmodel)$coef["GenotypeTG",2]
res[4,4] <- summary(lmodel)$coef["GenotypeTG",4]

res[5,4] <- summary(lmodel)$coef["Age",1]
res[6,4] <- summary(lmodel)$coef["Age",2]
res[7,4] <- summary(lmodel)$coef["Age",4]

res[8,4] <- summary(lmodel)$coef["Pathology",1]
res[9,4] <- summary(lmodel)$coef["Pathology",2]
res[10,4] <- summary(lmodel)$coef["Pathology",4]


write.csv(res, "DNAmAgeClockCortex_GenotypeAgePathology_results.csv")

