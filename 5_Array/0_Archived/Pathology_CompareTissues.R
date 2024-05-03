## Aisha Dahir A.N.Dahir@exeter.ac.uk
## 14.07.20
## This script compares the results of the two tissues Cortex and Hippocampus for Pathology

library(tidyverse)
################### Pathology ################




print("#################################### rTg4510 ################################")
#################################### rTg4510 ################################




###  CORTEX 
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_rTg4510.RData")
betaResults <- betaResultsPathology
rm(list=setdiff(ls(), c("betaResults"))) 

#Load coordinates matrix
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #24048 probes 
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- species_probes$MusMusculus
species_probes$Bp <- gsub("*.:","", species_probes$Bp)
species_probes <- species_probes[,-2]
species_probes$Position <- paste("chr", paste(species_probes$Chr, species_probes$Bp, sep = ":"), sep ="")
Coordinates <- species_probes

# Add chr and bp before deleting rows
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]
betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
betaResults <- betaResults[,c(10,1:9)]
# Add FDR pvalues
FDR_adj_Pathology <- p.adjust(betaResults[,"p.val.Pathology"], method = "fdr")
betaResults <- cbind(betaResults, FDR_adj_Pathology)
betaResults <- betaResults[,c(1:4,11,5:10)]
# Convert to data frame w/out factors
betaResults2 <- as.data.frame(betaResults)
betaResults2[] <- lapply(betaResults2, as.character)
for(i in 4:11){
  betaResults2[,i] <- as.numeric(as.character(betaResults2[,i]))
}
betaResults2[,"Bp"] <- as.numeric(as.character(betaResults2[,"Bp"]))
betaResults2[,"cpg"] <- as.character(betaResults2[,"cpg"]) 
betaResults <- betaResults2
# Convert x and y to numbers and remove non autosomal
betaResults<- betaResults[-which(betaResults$Chr %in% c("CHR_MG51_PATCH",
                                                      "CHR_MG4200_PATCH","CHR_MG3699_PATCH")),]
betaResults$Chr <- gsub("X", 20, betaResults$Chr)
betaResults$Chr <- gsub("Y", 21, betaResults$Chr)
betaResults$Chr <- as.numeric(betaResults$Chr)
betaResults$Bp <- as.numeric(betaResults$Bp)
### Add gene names in
mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)
betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]
betaResultsCortex <- betaResults

### Hippocampus
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_rTg4510_HIP.RData")
betaResults <- betaResultsPathology

# Add chr and bp before deleting rows
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]
betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
betaResults <- betaResults[,c(10,1:9)]
# Add FDR pvalues
FDR_adj_Pathology <- p.adjust(betaResults[,"p.val.Pathology"], method = "fdr")
betaResults <- cbind(betaResults, FDR_adj_Pathology)
betaResults <- betaResults[,c(1:4,11,5:10)]
# Convert to data frame w/out factors
betaResults2 <- as.data.frame(betaResults)
betaResults2[] <- lapply(betaResults2, as.character)
for(i in 4:11){
  betaResults2[,i] <- as.numeric(as.character(betaResults2[,i]))
}
betaResults2[,"Bp"] <- as.numeric(as.character(betaResults2[,"Bp"]))
betaResults2[,"cpg"] <- as.character(betaResults2[,"cpg"]) 
betaResults <- betaResults2
# Convert x and y to numbers and remove non autosomal
betaResults<- betaResults[-which(betaResults$Chr %in% c("CHR_MG51_PATCH",
                                                      "CHR_MG4200_PATCH","CHR_MG3699_PATCH")),]
betaResults$Chr <- gsub("X", 20, betaResults$Chr)
betaResults$Chr <- gsub("Y", 21, betaResults$Chr)
betaResults$Chr <- as.numeric(betaResults$Chr)
betaResults$Bp <- as.numeric(betaResults$Bp)
### Add gene names in
mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)
betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]

betaResultsHip <- betaResults
rm(list=setdiff(ls(), c("betaResultsCortex", "betaResultsHip"))) 


# Table
#Make a table that lists the number os significant sites and how many are shared in both tissues
betaResultsCortexSig <- betaResultsCortex[which(betaResultsCortex$FDR_adj_Pathology <= 0.05),]
betaResultsHipSig <- betaResultsHip[which(betaResultsHip$FDR_adj_Pathology <= 0.05),]

PathologyTbl <- matrix(NA, nrow = 2, ncol = 3)
PathologyTbl[1,1] <- nrow(betaResultsCortex)
PathologyTbl[1,2] <- nrow(betaResultsCortexSig)
PathologyTbl[1,3] <- length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg))
PathologyTbl[2,1] <- nrow(betaResultsHip)
PathologyTbl[2,2] <- nrow(betaResultsHipSig)
PathologyTbl[2,3] <- length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg))

PathologyTbl <- as.data.frame(PathologyTbl)
rownames(PathologyTbl) <- c("Entorhinal Cortex", "Hippocampus")
colnames(PathologyTbl) <- c("sites", "significant sites", "shared significant sites")

library(gridExtra)
png("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_PathologyCortexHip.png")
p<-tableGrob(PathologyTbl)
grid.arrange(p)
dev.off()

## Pull out the shared significant sites....
sharesigsites <- intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)
CorSharedSigSites  <- betaResultsCortex[which(betaResultsCortex$cpg %in% sharesigsites),]
HipSharedSigSites  <- betaResultsHip[which(betaResultsHip$cpg %in% sharesigsites),]

SharedSigSites <- cbind(CorSharedSigSites, HipSharedSigSites)
SharedSigSites <- SharedSigSites[,c(1,12, 2,3,10,5,22,17)]
colnames(SharedSigSites) <- c("cpg","Gene_Symbol", "Chr","Bp","Cortex.estimate","Cortex.FDR_adj_Pathology",
                              "Hip.estimate","Hip.FDR_adj_Pathology")
SharedSigSites <- SharedSigSites[order(SharedSigSites$Cortex.estimate, decreasing = T),]

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_PathologySharedSigSites.pdf", w = 15, h = 50)
p<-tableGrob(SharedSigSites[1:150,])
grid.arrange(p)
dev.off()

write.csv(SharedSigSites, "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_PathologySharedSigSites.csv")


##Plots

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_PathologyEffectSizesCortexHip.pdf")
plot(betaResultsCortex$estimate.Pathology, betaResultsHip$estimate.Pathology, cex = 0.7,
	xlab="Entorhinal Cortex", ylab="Hippocampus", main="All Sites Effect Sizes")
abline(v=0, h=0)
dev.off()

#filter for hippocampus significant sites then plot the effect sizes again ....
betaResultsCortexSig2 <- betaResultsCortex[which(betaResultsCortex$cpg %in% betaResultsHipSig$cpg),]
identical(betaResultsCortexSig2$cpg, betaResultsHipSig$cpg)
#plot effect sizes but with hippocampus number of sig as there will be more points. Color red where sites are sig in cortex
if(nrow(betaResultsHipSig) > 0){
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_PathologyEffectSizesCortexHip_HipSignificant.pdf")
plot(betaResultsCortexSig2$estimate.Pathology, betaResultsHipSig$estimate.Pathology, cex = 0.7,
	col=ifelse(betaResultsCortexSig2$estimate.Pathology %in% betaResultsCortexSig$estimate.Pathology, "red", "black"),
	 main="Hippocampus Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
dev.off()
} else {
	print("Hippocampus has no sig hits so will not create hippocampus sig sites effect size plots")
}

# Cortex has no hits so cannot run this section
#filter for cortex significant sites then plot the effect sizes again ....
betaResultsHipSig2 <- betaResultsHip[which(betaResultsHip$cpg %in% betaResultsCortexSig$cpg),]
identical(betaResultsHipSig2$cpg, betaResultsCortexSig$cpg)
#plot effect sizes but with cortex number of sig. Color red where sig sites are shared
if(nrow(betaResultsCortexSig) > 0){
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_PathologyEffectSizesCortexHip_CortexSignificant.pdf")
plot(betaResultsCortexSig$estimate.Pathology, betaResultsHipSig2$estimate.Pathology, cex = 0.7,
	col=ifelse(betaResultsHipSig2$estimate.Pathology %in% betaResultsHipSig$estimate.Pathology, "red", "black"),
	 main="Cortex Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
dev.off()
} else {
	print("Cortex has no sig hits so will not create cortex sig sites effect size plots")
}


sharedsigsites <- intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)
sharedsigsites <- betaResultsCortex[which(betaResultsCortex$cpg %in% sharedsigsites), c("cpg", "Gene_Symbol", "estimate.Pathology", "FDR_adj_Pathology")]
sharedsigsites <- sharedsigsites[order(sharedsigsites$FDR_adj_Pathology),]
if(nrow(sharedsigsites) > 0){
png("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_PathologySharedSigSites.png")
p<-tableGrob(sharedsigsites)
grid.arrange(p)
dev.off()
} else{
	print("No rows for shared sites information")
}

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_PathologyPvals.pdf", onefile = T)
plot(betaResultsCortex$FDR_adj_Pathology, betaResultsHip$FDR_adj_Pathology)
plot(-log10(betaResultsCortex$FDR_adj_Pathology),-log10(betaResultsHip$FDR_adj_Pathology))
dev.off()


# Make one figure with all effect sizes in one pfd; 4 panels a) All sites b) Cortex sig sites c) Hip sig sites d) shares sig sites
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_PathologyEffectSizes.pdf", onefile = T)
par(mfrow = c(2,2))
# a) All sites
plot(betaResultsCortex$estimate.Pathology, betaResultsHip$estimate.Pathology, cex = 0.7,
	xlab="Entorhinal Cortex", ylab="Hippocampus", main="All Sites Effect Sizes")
abline(v=0, h=0)

# b) Cortex Sig sites
if(nrow(betaResultsCortexSig) > 0){
plot(betaResultsCortexSig$estimate.Pathology, betaResultsHipSig2$estimate.Pathology, cex = 0.7,
	col=ifelse(betaResultsHipSig2$estimate.Pathology %in% betaResultsHipSig$estimate.Pathology, "red", "black"),
	 main="Cortex Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
} else {
	print("Cortex has no sig hits so will not create cortex sig sites effect size plots")
}

# c) Hip Sig sites
if(nrow(betaResultsHipSig) > 0){
plot(betaResultsCortexSig2$estimate.Pathology, betaResultsHipSig$estimate.Pathology, cex = 0.7,
	col=ifelse(betaResultsCortexSig2$estimate.Pathology %in% betaResultsCortexSig$estimate.Pathology, "red", "black"),
	 main="Hippocampus Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
} else {
	print("Hippocampus has no sig hits so will not create hippocampus sig sites effect size plots")
}
# d) Sig only sites
if(length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)) > 0){
betaResultsCortexSig3 <- betaResultsCortex[which(betaResultsCortex$cpg %in% intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)),]
betaResultsHipSig3 <- betaResultsHip[which(betaResultsHip$cpg %in% intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)),]
plot(betaResultsCortexSig3$estimate.Pathology, betaResultsHipSig3$estimate.Pathology, cex = 0.7,
	 main="Shared Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
} else {
	print("4th plot of significant sites not plotted")
}
dev.off()








print("#################################### J20 ################################")
#################################### J20 ################################















###  CORTEX 
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_J20.RData")
betaResults <- betaResultsPathology
rm(list=setdiff(ls(), c("betaResults"))) 

#Load coordinates matrix
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #24048 probes 
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- species_probes$MusMusculus
species_probes$Bp <- gsub("*.:","", species_probes$Bp)
species_probes <- species_probes[,-2]
species_probes$Position <- paste("chr", paste(species_probes$Chr, species_probes$Bp, sep = ":"), sep ="")
Coordinates <- species_probes

# Add chr and bp before deleting rows
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]
betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
betaResults <- betaResults[,c(10,1:9)]
# Add FDR pvalues
FDR_adj_Pathology <- p.adjust(betaResults[,"p.val.Pathology"], method = "fdr")
betaResults <- cbind(betaResults, FDR_adj_Pathology)
betaResults <- betaResults[,c(1:4,11,5:10)]
# Convert to data frame w/out factors
betaResults2 <- as.data.frame(betaResults)
betaResults2[] <- lapply(betaResults2, as.character)
for(i in 4:11){
  betaResults2[,i] <- as.numeric(as.character(betaResults2[,i]))
}
betaResults2[,"Bp"] <- as.numeric(as.character(betaResults2[,"Bp"]))
betaResults2[,"cpg"] <- as.character(betaResults2[,"cpg"]) 
betaResults <- betaResults2
# Convert x and y to numbers and remove non autosomal
betaResults<- betaResults[-which(betaResults$Chr %in% c("CHR_MG51_PATCH",
                                                      "CHR_MG4200_PATCH","CHR_MG3699_PATCH")),]
betaResults$Chr <- gsub("X", 20, betaResults$Chr)
betaResults$Chr <- gsub("Y", 21, betaResults$Chr)
betaResults$Chr <- as.numeric(betaResults$Chr)
betaResults$Bp <- as.numeric(betaResults$Bp)
### Add gene names in
mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)
betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]
betaResultsCortex <- betaResults

### Hippocampus
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_J20_HIP.RData")
betaResults <- betaResultsPathology
# Add chr and bp before deleting rows
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]
betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
betaResults <- betaResults[,c(10,1:9)]
# Add FDR pvalues
FDR_adj_Pathology <- p.adjust(betaResults[,"p.val.Pathology"], method = "fdr")
betaResults <- cbind(betaResults, FDR_adj_Pathology)
betaResults <- betaResults[,c(1:4,11,5:10)]
# Convert to data frame w/out factors
betaResults2 <- as.data.frame(betaResults)
betaResults2[] <- lapply(betaResults2, as.character)
for(i in 4:11){
  betaResults2[,i] <- as.numeric(as.character(betaResults2[,i]))
}
betaResults2[,"Bp"] <- as.numeric(as.character(betaResults2[,"Bp"]))
betaResults2[,"cpg"] <- as.character(betaResults2[,"cpg"]) 
betaResults <- betaResults2
# Convert x and y to numbers and remove non autosomal
betaResults<- betaResults[-which(betaResults$Chr %in% c("CHR_MG51_PATCH",
                                                      "CHR_MG4200_PATCH","CHR_MG3699_PATCH")),]
betaResults$Chr <- gsub("X", 20, betaResults$Chr)
betaResults$Chr <- gsub("Y", 21, betaResults$Chr)
betaResults$Chr <- as.numeric(betaResults$Chr)
betaResults$Bp <- as.numeric(betaResults$Bp)
### Add gene names in
mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)
betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]

betaResultsHip <- betaResults

rm(list=setdiff(ls(), c("betaResultsCortex", "betaResultsHip"))) 


# Table
#Make a table that lists the number os significant sites and how many are shared in both tissues
betaResultsCortexSig <- betaResultsCortex[which(betaResultsCortex$FDR_adj_Pathology <= 0.05),]
betaResultsHipSig <- betaResultsHip[which(betaResultsHip$FDR_adj_Pathology <= 0.05),]

PathologyTbl <- matrix(NA, nrow = 2, ncol = 3)
PathologyTbl[1,1] <- nrow(betaResultsCortex)
PathologyTbl[1,2] <- nrow(betaResultsCortexSig)
PathologyTbl[1,3] <- length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg))
PathologyTbl[2,1] <- nrow(betaResultsHip)
PathologyTbl[2,2] <- nrow(betaResultsHipSig)
PathologyTbl[2,3] <- length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg))

PathologyTbl <- as.data.frame(PathologyTbl)
rownames(PathologyTbl) <- c("Entorhinal Cortex", "Hippocampus")
colnames(PathologyTbl) <- c("sites", "significant sites", "shared significant sites")

library(gridExtra)
png("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_PathologyCortexHip.png")
p<-tableGrob(PathologyTbl)
grid.arrange(p)
dev.off()

## Pull out the shared significant sites....
sharesigsites <- intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)
CorSharedSigSites  <- betaResultsCortex[which(betaResultsCortex$cpg %in% sharesigsites),]
HipSharedSigSites  <- betaResultsHip[which(betaResultsHip$cpg %in% sharesigsites),]

SharedSigSites <- cbind(CorSharedSigSites, HipSharedSigSites)
SharedSigSites <- SharedSigSites[,c(1,12, 2,3,10,5,22,17)]
colnames(SharedSigSites) <- c("cpg","Gene_Symbol", "Chr","Bp","Cortex.estimate","Cortex.FDR_adj_Pathology",
                              "Hip.estimate","Hip.FDR_adj_Pathology")

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_PathologySharedSigSites.pdf", w = 15)
p<-tableGrob(SharedSigSites)
grid.arrange(p)
dev.off()

write.csv(SharedSigSites, "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_PathologySharedSigSites.csv")


##Plots

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_PathologyEffectSizesCortexHip.pdf")
plot(betaResultsCortex$estimate.Pathology, betaResultsHip$estimate.Pathology, cex = 0.7,
	xlab="Entorhinal Cortex", ylab="Hippocampus", main="All Sites Effect Sizes")
abline(v=0, h=0)
dev.off()

#filter for hippocampus significant sites then plot the effect sizes again ....
betaResultsCortexSig2 <- betaResultsCortex[which(betaResultsCortex$cpg %in% betaResultsHipSig$cpg),]
identical(betaResultsCortexSig2$cpg, betaResultsHipSig$cpg)
betaResultsHipSig2 <- betaResultsHip[which(betaResultsHip$cpg %in% betaResultsCortexSig$cpg),]
identical(betaResultsHipSig2$cpg, betaResultsCortexSig$cpg)
#plot effect sizes but with hippocampus number of sig as there will be more points. Color red where sites are sig in cortex
if(nrow(betaResultsCortexSig) > 0){
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_PathologyEffectSizesCortexHip_CortexSignificant.pdf")
plot(betaResultsCortexSig$estimate.Pathology, betaResultsHipSig2$estimate.Pathology, cex = 0.7,
	col=ifelse(betaResultsHipSig2$estimate.Pathology %in% betaResultsHipSig$estimate.Pathology, "red", "black"),
	 main="Cortex Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
dev.off()
} else {
	print("Cortex has no sig hits so will not create cortex sig sites effect size plots")
}

#filter for cortex significant sites then plot the effect sizes again ....
#plot effect sizes but with cortex number of sig. Color red where sig sites are shared
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_PathologyEffectSizesCortexHip_CortexSignificant.pdf")
plot(betaResultsCortexSig$estimate.Pathology, betaResultsHipSig2$estimate.Pathology, cex = 0.7,
	col=ifelse(betaResultsHipSig2$estimate.Pathology %in% betaResultsHipSig$estimate.Pathology, "red", "black"),
	 main="Cortex Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
dev.off()

sharedsigsites <- intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)
sharedsigsites <- betaResultsCortex[which(betaResultsCortex$cpg %in% sharedsigsites), c("cpg", "Gene_Symbol", "estimate.Pathology", "FDR_adj_Pathology")]
sharedsigsites <- sharedsigsites[order(sharedsigsites$FDR_adj_Pathology),]
if(nrow(sharedsigsites) > 0){
png("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_PathologySharedSigSites.png")
p<-tableGrob(sharedsigsites)
grid.arrange(p)
dev.off()
} else{
	print("No rows for shared sites information")
}

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_PathologyPvals.pdf", onefile = T)
plot(betaResultsCortex$FDR_adj_Pathology, betaResultsHip$FDR_adj_Pathology)
plot(-log10(betaResultsCortex$FDR_adj_Pathology),-log10(betaResultsHip$FDR_adj_Pathology))
dev.off()


# Make one figure with all effect sizes in one pfd; 4 panels a) All sites b) Cortex sig sites c) Hip sig sites d) shares sig sites
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_PathologyEffectSizes.pdf", onefile = T)
par(mfrow = c(2,2))
# a) All sites
plot(betaResultsCortex$estimate.Pathology, betaResultsHip$estimate.Pathology, cex = 0.7,
	xlab="Entorhinal Cortex", ylab="Hippocampus", main="All Sites Effect Sizes")
abline(v=0, h=0)

# b) Cortex Sig sites
if(nrow(betaResultsCortexSig) > 0){
plot(betaResultsCortexSig$estimate.Pathology, betaResultsHipSig2$estimate.Pathology, cex = 0.7,
	col=ifelse(betaResultsHipSig2$estimate.Pathology %in% betaResultsHipSig$estimate.Pathology, "red", "black"),
	 main="Cortex Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
} else {
	print("Cortex has no sig hits so will not create cortex sig sites effect size plots")
}

# c) Hip Sig sites
if(nrow(betaResultsHipSig) > 0){
plot(betaResultsCortexSig2$estimate.Pathology, betaResultsHipSig$estimate.Pathology, cex = 0.7,
	col=ifelse(betaResultsCortexSig2$estimate.Pathology %in% betaResultsCortexSig$estimate.Pathology, "red", "black"),
	 main="Hippocampus Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
} else {
	print("Hippocampus has no sig hits so will not create hippocampus sig sites effect size plots")
}

# d) Sig only sites
if(length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)) > 0){
betaResultsCortexSig3 <- betaResultsCortex[which(betaResultsCortex$cpg %in% intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)),]
betaResultsHipSig3 <- betaResultsHip[which(betaResultsHip$cpg %in% intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)),]
plot(betaResultsCortexSig3$estimate.Pathology, betaResultsHipSig3$estimate.Pathology, cex = 0.7,
	 main="Shared Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
} else {
	print("4th plot of significant sites not plotted")
}
dev.off()