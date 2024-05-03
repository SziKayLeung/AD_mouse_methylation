## Comparison of raw methylation of cortex and hippocampus


setwd("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/")
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat") # unique probes
QCmetrics$SampleID <- gsub(".*_", "", QCmetrics$ExternalSampleID)



#################################### rTg4510 #################### 
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
### Hip
pheno_rTg4510_HIP <- QCmetrics[which(QCmetrics$Tissue == "Hippocampus" & QCmetrics$AD_model == "rTg4510"),]
betas_rTg4510_HIP <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno_rTg4510_HIP$Basename ]
identical(pheno_rTg4510_HIP$Basename, colnames(betas_rTg4510_HIP))
colnames(betas_rTg4510_HIP) <- pheno_rTg4510_HIP$SampleID
rownames(pheno_rTg4510_HIP) <- pheno_rTg4510_HIP$SampleID 
pheno_rTg4510_HIP <- pheno_rTg4510_HIP[rTg4510_path$Sample_ID_HIP, ] #order based on pathology samples id (NA rows will be generated)

# Match the same mouse with 2 brain regions
samples <- as.data.frame(cbind(pheno_rTg4510_ECX$SampleID, pheno_rTg4510_HIP$SampleID))
colnames(samples) <- c("ECX", "HIP")
samples <- samples[!is.na(samples$ECX),]
samples <- samples[!is.na(samples$HIP),]
samples$ECX <- as.character(samples$ECX)
samples$HIP <- as.character(samples$HIP)
## 36 matched samples!!

# Now filter these matching samples on both tissues
### ecx
pheno_rTg4510_ECX <- pheno_rTg4510_ECX[samples$ECX,]
identical(pheno_rTg4510_ECX$SampleID, samples$ECX)
betas_rTg4510_ECX <- betas_rTg4510_ECX[,pheno_rTg4510_ECX$SampleID]
identical(colnames(betas_rTg4510_ECX) , pheno_rTg4510_ECX$SampleID)
pheno_rTg4510_ECX$col <- ifelse(pheno_rTg4510_ECX$Genotype == "WT", "black", "#00AEC9")

### hip
pheno_rTg4510_HIP <- pheno_rTg4510_HIP[samples$HIP,]
identical(pheno_rTg4510_HIP$SampleID, samples$HIP)
betas_rTg4510_HIP <- betas_rTg4510_HIP[,pheno_rTg4510_HIP$SampleID]
identical(colnames(betas_rTg4510_HIP) , pheno_rTg4510_HIP$SampleID)
pheno_rTg4510_HIP$col <- ifelse(pheno_rTg4510_HIP$Genotype == "WT", "black", "#00AEC9")


############### Plots
pdf("rTg4510_ECXvsHIP.pdf")
for(i in 1:ncol(betas_rTg4510_ECX)){
plot(betas_rTg4510_ECX[,i], betas_rTg4510_HIP[,i],
     xlab = "Entorhinalis Cortex", ylab = "Hippocampus",
     main = paste(pheno_rTg4510_ECX[i,"SampleID"], pheno_rTg4510_HIP[i,"SampleID"], sep = " & "),
     cex = 0.75)
corr <- cor(betas_rTg4510_ECX[,i], betas_rTg4510_HIP[,i])
mtext(paste("Cor = ", signif(corr,3)), side = 3, adj = 1)
mtext(paste(pheno_rTg4510_ECX[i,"Genotype"]), side = 3, adj = 0)
}
dev.off()


#################################### J20 ####################

rm(list=setdiff(ls(), c("Normalised_Sesame_Betas", "QCmetrics")))

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
### Hip
pheno_J20_HIP <- QCmetrics[which(QCmetrics$Tissue == "Hippocampus" & QCmetrics$AD_model == "J20"),]
betas_J20_HIP <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno_J20_HIP$Basename ]
identical(pheno_J20_HIP$Basename, colnames(betas_J20_HIP))
colnames(betas_J20_HIP) <- pheno_J20_HIP$SampleID
rownames(pheno_J20_HIP) <- pheno_J20_HIP$SampleID 
pheno_J20_HIP <- pheno_J20_HIP[J20_path$Sample_ID_HIP, ] #order based on pathology samples id (NA rows will be generated)

# Match the same mouse with 2 brain regions
samples <- as.data.frame(cbind(pheno_J20_ECX$SampleID, pheno_J20_HIP$SampleID))
colnames(samples) <- c("ECX", "HIP")
samples <- samples[!is.na(samples$ECX),]
samples <- samples[!is.na(samples$HIP),]
samples$ECX <- as.character(samples$ECX)
samples$HIP <- as.character(samples$HIP)
## 39 matched samples!!

# Now filter these matching samples on both tissues
### ecx
pheno_J20_ECX <- pheno_J20_ECX[samples$ECX,]
identical(pheno_J20_ECX$SampleID, samples$ECX)
betas_J20_ECX <- betas_J20_ECX[,pheno_J20_ECX$SampleID]
identical(colnames(betas_J20_ECX) , pheno_J20_ECX$SampleID)
pheno_J20_ECX$col <- ifelse(pheno_J20_ECX$Genotype == "WT", "black", "#FF5A62")
### hip
pheno_J20_HIP <- pheno_J20_HIP[samples$HIP,]
identical(pheno_J20_HIP$SampleID, samples$HIP)
betas_J20_HIP <- betas_J20_HIP[,pheno_J20_HIP$SampleID]
identical(colnames(betas_J20_HIP) , pheno_J20_HIP$SampleID)


############### Plots
pdf("J20_ECXvsHIP.pdf")
for(i in 1:ncol(betas_J20_ECX)){
plot(betas_J20_ECX[,i], betas_J20_HIP[,i],
       xlab = "Entorhinalis Cortex", ylab = "Hippocampus",
       main = paste(pheno_J20_ECX[i,"SampleID"], pheno_J20_HIP[i,"SampleID"], sep = " & "),
       cex = 0.75)
corr <- cor(betas_J20_ECX[,i], betas_J20_HIP[,i])
mtext(paste("Cor = ", signif(corr,3)), side = 3, adj = 1)
mtext(paste(pheno_J20_ECX[i,"Genotype"]), side = 3, adj = 0)
}
dev.off()

