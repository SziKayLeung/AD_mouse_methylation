library(tidyr)
library(dplyr)
library(ggplot2)

## read data in
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame.rdat") 
#load this in to have passed QC samples only as samplesheet has failed samples
pheno <- QCmetrics[which(QCmetrics$Tissue %in% c("CortexEntorhinalis", "Hippocampus")),] #take tissue samples further


samplesheet <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Horvath/MammalianArrayNormalizationTools/Clocks/OutputMouseClocks.csv", header = T,stringsAsFactors = F)
samplesheet <- samplesheet[,c(1:50)] #filter unnecessary values
samplesheet <- samplesheet[which(samplesheet$Basename %in% pheno$Basename),] #filter to just tissue samples to join with pathology data
samplesheet$SampleID <- gsub(".*_", "",samplesheet$ExternalSampleID)
## Add pathology data in
J20_pathdata <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_coldata_VertebrateArray.csv", header = T, stringsAsFactors = F)
rTg4510_pathdata <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Tg4510_coldata_VertebrateArray.csv", header = T, stringsAsFactors = F)

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

res <- matrix(NA, nrow = 8, ncol= 4)
colnames(res) <- c("rTg4510_ECX", "rTg4510_HIP", "J20_ECX", "J20_HIP")
rownames(res) <- c("Genotype.Intercept", "Genotype.Beta", "Genotype.SE", "Genotype.P",
                   "Pathology.Intercept","Pathology.Beta", "Pathology.SE", "Pathology.P")

samplesheet$Genotype <- factor(samplesheet$Genotype, levels = c("WT","TG"))

rTg4510ECX <- samplesheet[which(samplesheet$Tissue == "CortexEntorhinalis" & samplesheet$AD_model == "rTg4510"),]
rTg4510HIP <- samplesheet[which(samplesheet$Tissue == "Hippocampus" & samplesheet$AD_model == "rTg4510"),]

rTg4510Clocks <- list(rTg4510ECX, rTg4510HIP)
names(rTg4510Clocks) <- c("ECX","HIP")

zenDir <- "/lustre/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/0_ZenOutput"
save(rTg4510Clocks, file = paste0(zenDir, "/5_epigeneticClock/rTg4510Clock.RData"))
