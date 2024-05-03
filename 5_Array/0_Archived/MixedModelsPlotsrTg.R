## 04/08/2020
## Run plots for mixed models effects
library(tidyverse)
library(qqman)
library(pheatmap)
source("/gpfs/ts0/projects/Research_Project-191406/Aisha/script/rrbs-ad-mice/PlotFunctions.R")
setwd("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/")


#get raw data (betas) & pheno
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat")
QCmetrics$SampleID <- gsub(".*_", "", QCmetrics$ExternalSampleID)
###################################### rTg4510  #######################################
rTg4510_path <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Tg4510_coldata_VertebrateArray.csv",
                         header = T, stringsAsFactors = F)
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
pheno <- pheno[,c("Tissue","Genotype", "Age_months", "Chip_ID", "PlateNumber", "MouseID", "Basename", "Pathology_ECX", "Pathology_HIP")]
betas <- as.matrix(betas)
rownames(pheno) <- pheno$Basename
colnames(betas) <- rownames(pheno)

##############################################################################################


## Get the results data in now
rTg4510 <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/MixedModelResults_rTg4510.csv",
                    header=T, stringsAsFactors = F)

mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)

  

colnames(rTg4510)[colnames(rTg4510) == 'X'] <- 'cpg'



# Add annotation files - chr and bp
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_mappability_name = "MusMusculus"
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #23633
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- gsub(".*:","",species_probes$MusMusculus)

rTg4510$Chr <- species_probes$Chr[match(species_probes$probeID, rTg4510[,"cpg"])]
rTg4510$Bp <- species_probes$Bp[match(species_probes$probeID, rTg4510[,"cpg"])]


# Convert x and y to numbers and remove non autosomal
rTg4510<- rTg4510[-which(rTg4510$Chr %in% c("CHR_MG51_PATCH", "CHR_MG4200_PATCH","CHR_MG3699_PATCH")),]
rTg4510$Chr <- gsub("X", 20, rTg4510$Chr)
rTg4510$Chr <- gsub("Y", 21, rTg4510$Chr)
rTg4510$Chr <- as.numeric(rTg4510$Chr)
rTg4510$Bp <- as.numeric(rTg4510$Bp)
rTg4510$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(rTg4510$cpg, mm10_Manifest$cpg)]
rownames(rTg4510) <- rTg4510$cpg
rTg4510$SNP <- rTg4510$Gene_Symbol



#manhattan plots & qqplots
threshold <- 0.05
color_Tg4510_TG <- "#00AEC9"
color_J20_TG <- "#FF5A62"

#rTg4510
pdf("MME_rTg4510_qqplots.pdf")
par(mfrow = c(2,2))
lamda <- qchisq(1-median(rTg4510$PrZ.Tissue),1)/qchisq(0.5,1)
qq(rTg4510$PrZ.Tissue, main = "rTg4510 Tissue")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(rTg4510$PrZ.GenotypeTG),1)/qchisq(0.5,1)
qq(rTg4510$PrZ.GenotypeTG, main = "rTg4510 Genotype")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(rTg4510$PrZ.TissueHippocampus.GenotypeTG),1)/qchisq(0.5,1)
qq(rTg4510$PrZ.TissueHippocampus.GenotypeTG, main = "rTg4510 Tissue*Genotype")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(rTg4510$PrZ.GenotypeTG.Age_months),1)/qchisq(0.5,1)
qq(rTg4510$PrZ.GenotypeTG.Age_months, main = "rTg4510 Genotype*Age")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)
dev.off()

rTg4510 <- rTg4510[order(rTg4510$PrZ.GenotypeTG),]
rTg4510_2 <- rTg4510[-1,]
manhattan(rTg4510_2, p = "PrZ.GenotypeTG", bp = "Bp", chr = "Chr", 
          genomewide = -log10(threshold), suggestiveline = FALSE, 
          logp=T, col = c("black", color_Tg4510_TG), main = "Genotype",
          chrlabs = c(1:19, "X", "Y"),las = 2, annotatePval = threshold,
          cex.axis = 1.5, cex.lab = 1.2)

SigrTg4510_2 <- rTg4510_2[which(rTg4510_2$FDR_adj_tissue < 0.05),] #17329

## pull out numbers etc
rTg4510$FDR_adj_GenotypeTG <- p.adjust(rTg4510[,"PrZ.GenotypeTG"], method = "fdr")
rTg4510$FDR_adj_GenotypeAge <- p.adjust(rTg4510[,"PrZ.GenotypeTG.Age_months"], method = "fdr")
rTg4510$FDR_adj_TissueGenotype <- p.adjust(rTg4510[,"PrZ.TissueHippocampus.GenotypeTG"], method = "fdr")


nrow(rTg4510[which(rTg4510$FDR_adj_GenotypeTG < 0.05),])     #75
nrow(rTg4510[which(rTg4510$FDR_adj_GenotypeAge < 0.05),])    #1971
nrow(rTg4510[which(rTg4510$FDR_adj_TissueGenotype < 0.05),]) #5229


######################## Are the effects direction and magnitude the same as within seperate tissues?
#rTg4510
#cortex 
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_rTg4510.RData")
betaResults <- betaResultsChip
# Add chr and bp before deleting rows
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]
betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
betaResults <- betaResults[,c(17,1:16)]

# Add FDR pvalues
FDR_adj_genotype <- p.adjust(betaResults[,"p.val.Genotype"], method = "fdr")
FDR_adj_age <- p.adjust(betaResults[,"p.val.Age"], method = "fdr")
FDR_adj_Interaction <- p.adjust(betaResults[,"p.val.Interaction"], method = "fdr")
betaResults <- cbind(betaResults, FDR_adj_genotype, FDR_adj_age, FDR_adj_Interaction)
betaResults <- betaResults[,c(1:10,18,11:13,19,14:17,20)]

# Convert to data frame w/out factors
betaResults2 <- as.data.frame(betaResults)
betaResults2[] <- lapply(betaResults2, as.character)
for(i in 4:20){
  betaResults2[,i] <- as.numeric(as.character(betaResults2[,i]))
}
betaResults2[,"Bp"] <- as.numeric(as.character(betaResults2[,"Bp"]))
betaResults <- betaResults2
# Convert x and y to numbers and remove non autosomal
betaResults<- betaResults[-which(betaResults$Chr %in% c("CHR_MG51_PATCH",
                                                        "CHR_MG4200_PATCH","CHR_MG3699_PATCH")),]
betaResults$Chr <- gsub("X", 20, betaResults$Chr)
betaResults$Chr <- gsub("Y", 21, betaResults$Chr)
betaResults$Chr <- as.numeric(betaResults$Chr)
betaResults$Bp <- as.numeric(betaResults$Bp)

### Add gene names in
betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]
rownames(betaResults) <- betaResults$cpg
rTg4510_cortex <- betaResults


#rTg4510
#hippocampus
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_rTg4510_HIP.RData")
betaResults <- betaResultsChip
# Add chr and bp before deleting rows
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]

betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
betaResults <- betaResults[,c(17,1:16)]
# Add FDR pvalues
FDR_adj_genotype <- p.adjust(betaResults[,"p.val.Genotype"], method = "fdr")
FDR_adj_age <- p.adjust(betaResults[,"p.val.Age"], method = "fdr")
FDR_adj_Interaction <- p.adjust(betaResults[,"p.val.Interaction"], method = "fdr")
betaResults <- cbind(betaResults, FDR_adj_genotype, FDR_adj_age, FDR_adj_Interaction)
betaResults <- betaResults[,c(1:10,18,11:13,19,14:17,20)]
# Convert to data frame w/out factors
betaResults2 <- as.data.frame(betaResults)
betaResults2[] <- lapply(betaResults2, as.character)
for(i in 4:20){
  betaResults2[,i] <- as.numeric(as.character(betaResults2[,i]))
}
betaResults2[,"Bp"] <- as.numeric(as.character(betaResults2[,"Bp"]))
betaResults <- betaResults2
# Convert x and y to numbers and remove non autosomal
betaResults<- betaResults[-which(betaResults$Chr %in% c("CHR_MG51_PATCH",
                                                        "CHR_MG4200_PATCH","CHR_MG3699_PATCH")),]
betaResults$Chr <- gsub("X", 20, betaResults$Chr)
betaResults$Chr <- gsub("Y", 21, betaResults$Chr)
betaResults$Chr <- as.numeric(betaResults$Chr)
betaResults$Bp <- as.numeric(betaResults$Bp)

### Add gene names in
betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]
rownames(betaResults) <- betaResults$cpg
rTg4510_hipp <- betaResults


############################### Comparing effect sizes - no filtering for sig sites
# Check that they are all in the same order of rows so they are comparing the same things
rTg4510 <- rTg4510[order(rTg4510$cpg),]
rTg4510_cortex <- rTg4510_cortex[order(rTg4510_cortex$cpg),]
rTg4510_hipp <- rTg4510_hipp[order(rTg4510_hipp$cpg),]

identical(rTg4510$cpg, rTg4510_cortex$cpg)
identical(rTg4510$cpg, rTg4510_hipp$cpg)

pdf("rTg_AllEffectSizes_braincortexhippocampus.pdf")
par(mfrow = c(2,2))
# cortex vs brain
cor.test(rTg4510$Betas.GenotypeTG, rTg4510_cortex$estimate.Genotype, method = "pearson")
corr <- cor(rTg4510$Betas.GenotypeTG, rTg4510_cortex$estimate.Genotype, method = "pearson")
plot(rTg4510$Betas.GenotypeTG, rTg4510_cortex$estimate.Genotype, 
     xlab = "Cortex", ylab = "Brain",
     main = "Effect sizes: rTg4510", cex = 0.7)
mtext(paste("cor = ", signif(corr,3)), side = 3, adj = 1)
abline(0,1, col="blue")
abline(h = 0)
abline(v = 0)

#hipp vs brain
cor.test(rTg4510$Betas.GenotypeTG, rTg4510_hipp$estimate.Genotype, method = "pearson")
corr <- cor(rTg4510$Betas.GenotypeTG, rTg4510_hipp$estimate.Genotype, method = "pearson")
plot(rTg4510$Betas.GenotypeTG, rTg4510_hipp$estimate.Genotype, 
     xlab = "Hippocampus", ylab = "Brain",
     main = "Effect sizes: rTg4510", cex = 0.7)
mtext(paste("cor = ", signif(corr,3)), side = 3, adj = 1)
abline(0,1, col="blue")
abline(h = 0)
abline(v = 0)

# cortex vs hipp
cor.test(rTg4510_cortex$estimate.Genotype, rTg4510_hipp$estimate.Genotype, method = "pearson")
corr <- cor(rTg4510_cortex$estimate.Genotype, rTg4510_hipp$estimate.Genotype, method = "pearson")
plot(rTg4510_cortex$estimate.Genotype, rTg4510_hipp$estimate.Genotype, 
     xlab = "Cortex", ylab = "Hippocampus",
     main = "Effect sizes: rTg4510", cex = 0.7)
mtext(paste("cor = ", signif(corr,3)), side = 3, adj = 1)
abline(0,1, col="blue")
abline(h = 0)
abline(v = 0)
dev.off()

############################### Comparing effect sizes -  filtering for sig sites
# Add FDR P value to brain
rTg4510$FDR_adj_GenotypeTG <- p.adjust(rTg4510[,"PrZ.GenotypeTG"], method = "fdr")
rTg4510_sig <- rTg4510[which(rTg4510$FDR_adj_GenotypeTG <= 0.05), "cpg"] #75 sites for brain
rTg4510_cortex_sig <- rTg4510_cortex[which(rTg4510_cortex$FDR_adj_genotype <= 0.05), "cpg"] #47 sites for cortex
rTg4510_hipp_sig <- rTg4510_hipp[which(rTg4510_hipp$FDR_adj_genotype <= 0.05), "cpg"] #373 sites forr hip

length(intersect(rTg4510_sig, rTg4510_cortex_sig)) #25 for brain & cortex
length(intersect(rTg4510_sig, rTg4510_hipp_sig)) #11 for brain & hippocampus
length(intersect(rTg4510_cortex_sig, rTg4510_hipp_sig)) #6 for cortex & hippocampus

cortexbrain <- intersect(rTg4510_sig, rTg4510_cortex_sig)
hippbrain <- intersect(rTg4510_sig, rTg4510_hipp_sig)
cortexhipp <- intersect(rTg4510_cortex_sig, rTg4510_hipp_sig)

pdf("rTg_SignificantEffectSizes_braincortexhippocampus.pdf")
par(mfrow = c(2,2))
# cortex vs brain
sig_rTg4510_cortex <- rTg4510_cortex[which(rTg4510_cortex$cpg %in% cortexbrain),]
sig_rTg4510 <- rTg4510[which(rTg4510$cpg %in% cortexbrain),]
corr <- cor(sig_rTg4510$Betas.GenotypeTG, sig_rTg4510_cortex$estimate.Genotype, method = "pearson")
plot(sig_rTg4510$Betas.GenotypeTG, sig_rTg4510_cortex$estimate.Genotype, 
     xlab = "Cortex", ylab = "Brain",
     main = "Effect sizes: rTg4510", cex = 0.7,
     xlim = c(-1.1,0.6))
mtext(paste("cor = ", signif(corr,3)), side = 3, adj = 1)
abline(0,1, col="blue")
abline(h = 0)
abline(v = 0)
text(sig_rTg4510$Betas.GenotypeTG ~ sig_rTg4510_cortex$estimate.Genotype, labels = ifelse(sig_rTg4510_cortex$FDR_adj_genotype < 7.599799e-03 ,sig_rTg4510_cortex$Gene_Symbol, NA)
     ,cex=0.8, font=2)

#hipp vs brain
sig_rTg4510_hipp <- rTg4510_hipp[which(rTg4510_hipp$cpg %in% hippbrain),]
sig_rTg4510 <- rTg4510[which(rTg4510$cpg %in% hippbrain),]
corr <- cor(sig_rTg4510$Betas.GenotypeTG, sig_rTg4510_hipp$estimate.Genotype, method = "pearson")
plot(sig_rTg4510$Betas.GenotypeTG, sig_rTg4510_hipp$estimate.Genotype, 
     xlab = "Cortex", ylab = "Brain",
     main = "Effect sizes: rTg4510", cex = 0.7)
mtext(paste("cor = ", signif(corr,3)), side = 3, adj = 1)
abline(0,1, col="blue")
abline(h = 0)
abline(v = 0)
text(sig_rTg4510$Betas.GenotypeTG ~ sig_rTg4510_hipp$estimate.Genotype,
     labels = ifelse(sig_rTg4510_hipp$FDR_adj_genotype < 8.849399e-03 ,sig_rTg4510_hipp$Gene_Symbol, NA)
     ,cex=0.8, font=2)

# cortex vs hipp
sig_rTg4510_cortex <- rTg4510_cortex[which(rTg4510_cortex$cpg %in% cortexhipp),]
sig_rTg4510_hipp <- rTg4510_hipp[which(rTg4510_hipp$cpg %in% cortexhipp),]
cor.test(sig_rTg4510_cortex$estimate.Genotype, sig_rTg4510_hipp$estimate.Genotype, method = "pearson")
corr <- cor(sig_rTg4510_cortex$estimate.Genotype, sig_rTg4510_hipp$estimate.Genotype, method = "pearson")
plot(sig_rTg4510_cortex$estimate.Genotype, sig_rTg4510_hipp$estimate.Genotype, 
      xlab = "Cortex", ylab = "Hippocampus",
      main = "Effect sizes: rTg4510", cex = 0.7)
mtext(paste("cor = ", signif(corr,3)), side = 3, adj = 1)
abline(0,1, col="blue")
abline(h = 0)
abline(v = 0)
text(sig_rTg4510_cortex$estimate.Genotype, sig_rTg4510_hipp$estimate.Genotype,
     labels = sig_rTg4510_hipp$Gene_Symbol,cex=0.8, font=2, pos = 4)

dev.off()


rTg4510<-rTg4510[order(rTg4510$FDR_adj_GenotypeTG),]
rTg4510_sig <- rTg4510[1:length(rTg4510[which(rTg4510$FDR_adj_GenotypeTG <= 0.05), "cpg"]),]
betas_sig <- betas[rownames(rTg4510_sig),]
pheno_plot <- pheno[,c("Tissue","Genotype", "Age_months")]

pdf("rTg_Heatmap_sig_genotype.pdf")
pheatmap(betas_sig,  
         annotation_col = pheno_plot,
         cluster_rows = T,
         show_rownames = F,
         show_colnames = F)
dev.off()

####################### Interactions
# how many significant interactions
rTg4510$FDR_adj_TissueGenoInt <- p.adjust(rTg4510[,"PrZ.TissueHippocampus.GenotypeTG"], method = "fdr")
sig_tissueint <- rTg4510[which(rTg4510$FDR_adj_TissueGenoInt < 0.05),] 
nrow(sig_tissueint)

# How many are significant for interaction and genotype??
length(intersect(sig_tissueint$cpg, rTg4510_cortex_sig)) #22 sites for cortex
length(intersect(sig_tissueint$cpg, rTg4510_hipp_sig)) #211 sites for hip


intcortex <- intersect(sig_tissueint$cpg, rTg4510_cortex_sig) 
inthipp <- intersect(sig_tissueint$cpg, rTg4510_hipp_sig)



sig_tissueint <- sig_tissueint[order(sig_tissueint$PrZ.TissueHippocampus.GenotypeTG),]

betaResults <- sig_tissueint
cpg = rownames(betaResults[1,])
data <- betas
coldata <- pheno
index <- cpg
ages <- pheno$Age_months
row <- betaResults[index,]
site <- row$cpg
gene <- row$Gene_Symbol
shapes = c(0,1) 
shapes <- shapes[as.numeric(coldata$Tissue)]

data_plot <- cbind(coldata[,c("Age_months","Genotype", "Tissue")], data[site,])
colnames(data_plot) <- c("Age_months","Genotype", "Tissue", site)
data_plot <- data_plot[complete.cases(data_plot), ]
ggplot(data_plot, aes(x=Genotype, y=data_plot[,4], color=Genotype)) + 
  geom_point() +
  facet_wrap(~Tissue) +
  labs(y = site, title = betaResults$Gene_Symbol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Heatmaps for genotype & tissue*genotype

#This chunk of code generates the meta data so pheatmap knows what each probe is
rTg4510<-rTg4510[order(rTg4510$PrZ.TissueHippocampus.GenotypeTG),]
rTg4510_sig <- rTg4510[1:nrow(sig_tissueint),]
betas_sig <- betas[rownames(rTg4510_sig),]
pheno_plot <- pheno[,c("Tissue","Genotype", "Age_months")]

pdf("rTg_Heatmap_sig_tissuegenotype.pdf")
pheatmap(betas_sig,  
         annotation_col = pheno_plot,
         cluster_rows = T,
         show_rownames = F,
         show_colnames = F)
dev.off()


#### Genotype*Age
rTg4510$FDR_adj_GenotypeAgeInt <- p.adjust(rTg4510[,"PrZ.GenotypeTG.Age_months"], method = "fdr")
sig_genoageint <- rTg4510[which(rTg4510$FDR_adj_GenotypeAgeInt < 0.05),] 
nrow(sig_genoageint) #1971



############################################# Pathology data ##################################################################

rTg4510 <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/MixedModelResultswPathology_rTg4510.csv",
                    header=T, stringsAsFactors = F)

mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)


colnames(rTg4510)[colnames(rTg4510) == 'X'] <- 'cpg'



# Add annotation files - chr and bp
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_mappability_name = "MusMusculus"
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #23633
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- gsub(".*:","",species_probes$MusMusculus)

rTg4510$Chr <- species_probes$Chr[match(species_probes$probeID, rTg4510[,"cpg"])]
rTg4510$Bp <- species_probes$Bp[match(species_probes$probeID, rTg4510[,"cpg"])]


# Convert x and y to numbers and remove non autosomal
rTg4510<- rTg4510[-which(rTg4510$Chr %in% c("CHR_MG51_PATCH", "CHR_MG4200_PATCH","CHR_MG3699_PATCH")),]
rTg4510$Chr <- gsub("X", 20, rTg4510$Chr)
rTg4510$Chr <- gsub("Y", 21, rTg4510$Chr)
rTg4510$Chr <- as.numeric(rTg4510$Chr)
rTg4510$Bp <- as.numeric(rTg4510$Bp)
rTg4510$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(rTg4510$cpg, mm10_Manifest$cpg)]
rownames(rTg4510) <- rTg4510$cpg
rTg4510$SNP <- rTg4510$Gene_Symbol


par(mfrow = c(2,2))
lamda <- qchisq(1-median(rTg4510$PrZ.GenotypeTG),1)/qchisq(0.5,1)
qq(rTg4510$PrZ.GenotypeTG, main = "rTg4510 Genotype")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(rTg4510$PrZ.TissueHippocampus.GenotypeTG),1)/qchisq(0.5,1)
qq(rTg4510$PrZ.TissueHippocampus.GenotypeTG, main = "rTg4510 Tissue*Genotype")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(rTg4510$PrZ.Pathology_ECX),1)/qchisq(0.5,1)
qq(rTg4510$PrZ.Pathology_ECX, main = "rTg4510 Pathology_ECX")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(rTg4510$PrZ.Pathology_HIP),1)/qchisq(0.5,1)
qq(rTg4510$PrZ.Pathology_HIP, main = "rTg4510 Pathology_HIP")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)


## pull out numbers etc
rTg4510$FDR_adj_GenotypeTG <- p.adjust(rTg4510[,"PrZ.GenotypeTG"], method = "fdr")
rTg4510$FDR_adj_TissueGenotype <- p.adjust(rTg4510[,"PrZ.TissueHippocampus.GenotypeTG"], method = "fdr")
rTg4510$FDR_adj_Pathology_ECX <- p.adjust(rTg4510[,"PrZ.Pathology_ECX"], method = "fdr")
rTg4510$FDR_adj_Pathology_HIP <- p.adjust(rTg4510[,"PrZ.Pathology_HIP"], method = "fdr")


nrow(rTg4510[which(rTg4510$FDR_adj_GenotypeTG < 0.05),])     #116
nrow(rTg4510[which(rTg4510$FDR_adj_TissueGenotype < 0.05),]) #5193
nrow(rTg4510[which(rTg4510$FDR_adj_Pathology_ECX < 0.05),])  #1
nrow(rTg4510[which(rTg4510$FDR_adj_Pathology_HIP < 0.05),])  #855



##### significant genotype
rTg4510_genosig <- rTg4510[which(rTg4510$FDR_adj_GenotypeTG < 0.05),]
rTg4510_genosig <- rTg4510_genosig[order(rTg4510_genosig$Betas.GenotypeTG),]
rTg4510_genosig_table <- rTg4510_genosig[,c("cpg", "Gene_Symbol", "Betas.GenotypeTG","SE.GenotypeTG","PrZ.GenotypeTG","FDR_adj_GenotypeTG")]
colnames(rTg4510_genosig_table) <- c("cpg", "Gene_Symbol", "Beta","SE", "Raw.Pvalues", "FDR.Adj.Pvalues")
write.csv(rTg4510_genosig_table, "rTg4510_MME_Genotype.csv")

pdf("rTg4510_Genotype_MME.pdf", onefile = T)
par(mfrow = c(2,2))
for( i in 1:6){
  par(mfrow = c(2,2))
  betaResults <- rTg4510_genosig
  cpg = rownames(betaResults[i,])
  data <- betas
  coldata <- pheno
  index <- cpg
  ages <- pheno$Age_months
  row <- betaResults[index,]
  site <- row$cpg
  gene <- row$Gene_Symbol
  shapes = c(0,1) 
  shapes <- shapes[as.numeric(coldata$Tissue)]
  data_plot <- cbind(coldata[,c("Age_months","Genotype", "Tissue")], data[site,])
  colnames(data_plot) <- c("Age_months","Genotype", "Tissue", site)
  data_plot <- data_plot[complete.cases(data_plot), ]
  
  p <- ggplot(data_plot, aes(x=Genotype, y=data_plot[,4], color=Genotype)) + 
    geom_point() +
    facet_wrap(~Tissue) +
    scale_color_manual(values=c("black", "#00AEC9")) +
    labs(y = "Methylation ", title = paste(betaResults$Gene_Symbol[i], betaResults$cpg[i], sep = " - ")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(p)
}
dev.off()


##### significant genotype*tissue
rTg4510_genotisssig <- rTg4510[which(rTg4510$FDR_adj_TissueGenotype < 0.05),]
rTg4510_genotisssig <- rTg4510_genotisssig[order(rTg4510_genotisssig$Betas.TissueHippocampus.GenotypeTG),]
rTg4510_genotisssig_table <- rTg4510_genotisssig[,c( "cpg", "Gene_Symbol", "Betas.TissueHippocampus.GenotypeTG",
                                                     "SE.TissueHippocampus.GenotypeTG","PrZ.TissueHippocampus.GenotypeTG",
                                                     "FDR_adj_TissueGenotype")]
colnames(rTg4510_genotisssig_table) <- c("cpg", "Gene_Symbol", "Beta","SE", "Raw.Pvalues", "FDR.Adj.Pvalues")
write.csv(rTg4510_genotisssig_table, "rTg4510_MME_TissueGenotypeInteraction.csv")

pdf("rTg4510_GenotypeTissue_MME.pdf", onefile = T)
par(mfrow = c(2,2))
for( i in 1:6){
  par(mfrow = c(2,2))
  betaResults <- rTg4510_genotisssig
  cpg = rownames(betaResults[i,])
  data <- betas
  coldata <- pheno
  index <- cpg
  ages <- pheno$Age_months
  row <- betaResults[index,]
  site <- row$cpg
  gene <- row$Gene_Symbol
  shapes = c(0,1) 
  shapes <- shapes[as.numeric(coldata$Tissue)]
  data_plot <- cbind(coldata[,c("Age_months","Genotype", "Tissue")], data[site,])
  colnames(data_plot) <- c("Age_months","Genotype", "Tissue", site)
  data_plot <- data_plot[complete.cases(data_plot), ]
  
  p <- ggplot(data_plot, aes(x=Genotype, y=data_plot[,4], color=Genotype)) + 
    geom_point() +
    facet_wrap(~Tissue) +
    scale_color_manual(values=c("black", "#00AEC9")) +
    labs(y = "Methylation ", title = paste(betaResults$Gene_Symbol[i], betaResults$cpg[i], sep = " - ")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(p)
}
dev.off()



# how many significant pathology - ECX

rTg4510_pathecxsig <- rTg4510[which(rTg4510$FDR_adj_Pathology_ECX < 0.05),]
rTg4510_pathecxsig <- rTg4510_pathecxsig[order(rTg4510_pathecxsig$Betas.Pathology_ECX),]
rTg4510_pathecxsig_table <- rTg4510_pathecxsig[,c( "cpg", "Gene_Symbol", "Betas.Pathology_ECX","SE.Pathology_ECX","PrZ.Pathology_ECX", "FDR_adj_Pathology_ECX")]
colnames(rTg4510_pathecxsig_table) <- c("cpg", "Gene_Symbol", "Beta","SE", "Raw.Pvalues", "FDR.Adj.Pvalues")
write.csv(rTg4510_pathecxsig_table, "rTg4510_MME_Pathology_ECX.csv")

pdf("rTg4510_PathECX_MME.pdf", onefile = T)
par(mfrow = c(2,2))
for( i in 1:6){
  par(mfrow = c(2,2))
  betaResults <- rTg4510_pathecxsig
  cpg = rownames(betaResults[i,])
  data <- betas
  coldata <- pheno
  index <- cpg
  ages <- pheno$Age_months
  row <- betaResults[index,]
  site <- row$cpg
  gene <- row$Gene_Symbol
  shapes = c(0,1) 
  shapes <- shapes[as.numeric(coldata$Tissue)]
  data_plot <- cbind(coldata[,c("Age_months","Genotype", "Tissue","Pathology_ECX", "Pathology_HIP")], data[site,])
  colnames(data_plot) <- c("Age_months","Genotype", "Tissue","Pathology_ECX", "Pathology_HIP", site)
  data_plot <- data_plot[complete.cases(data_plot), ]
  
  p <- ggplot(data_plot, aes(x=Pathology_ECX, y=data_plot[,6], color=Genotype)) + 
    geom_point() +
    facet_wrap(~Tissue) +
    scale_color_manual(values=c("black", "#00AEC9")) +
    labs(y = "Methylation ", title = paste(betaResults$Gene_Symbol[i], betaResults$cpg[i], sep = " - ")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(p)
}
dev.off()


# how many significant pathology - HIP
rTg4510_pathhipsig <- rTg4510[which(rTg4510$FDR_adj_Pathology_HIP < 0.05),]
rTg4510_pathhipsig <- rTg4510_pathhipsig[order(rTg4510_pathhipsig$Betas.Pathology_HIP),]
rTg4510_pathhipsig_table <- rTg4510_pathhipsig[,c( "cpg", "Gene_Symbol", "Betas.Pathology_HIP","SE.Pathology_HIP","PrZ.Pathology_HIP", "FDR_adj_Pathology_HIP")]
colnames(rTg4510_pathhipsig_table) <- c("cpg", "Gene_Symbol", "Beta","SE", "Raw.Pvalues", "FDR.Adj.Pvalues")
write.csv(rTg4510_pathhipsig_table, "rTg4510_MME_Pathology_HIP.csv")

pdf("rTg4510_PathHIP_MME.pdf", onefile = T)
par(mfrow = c(2,2))
for( i in 1:6){
  par(mfrow = c(2,2))
  betaResults <- rTg4510_pathhipsig
  cpg = rownames(betaResults[i,])
  data <- betas
  coldata <- pheno
  index <- cpg
  ages <- pheno$Age_months
  row <- betaResults[index,]
  site <- row$cpg
  gene <- row$Gene_Symbol
  shapes = c(0,1) 
  shapes <- shapes[as.numeric(coldata$Tissue)]
  data_plot <- cbind(coldata[,c("Age_months","Genotype", "Tissue","Pathology_ECX", "Pathology_HIP")], data[site,])
  colnames(data_plot) <- c("Age_months","Genotype", "Tissue","Pathology_ECX", "Pathology_HIP", site)
  data_plot <- data_plot[complete.cases(data_plot), ]
  
  p <- ggplot(data_plot, aes(x=Pathology_HIP, y=data_plot[,6], color=Genotype)) + 
    geom_point() +
    facet_wrap(~Tissue) +
    scale_color_manual(values=c("black", "#00AEC9")) +
    labs(y = "Methylation ", title = paste(betaResults$Gene_Symbol[i], betaResults$cpg[i], sep = " - ")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(p)
}
dev.off()




