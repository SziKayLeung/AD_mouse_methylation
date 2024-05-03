## Aisha Dahir A.N.Dahir@exeter.ac.uk
## 14.07.20
## This script compares the results of the two tissues Cortex and Hippocampus for Interaction

library(tidyverse)
################### Interaction ################




print("#################################### rTg4510 ################################")
#################################### rTg4510 ################################




###  CORTEX 
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_rTg4510.RData")
betaResults <- betaResultsChip
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

#Convert the betaresults to the right format
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]
betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
# Remove NA rows
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
mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)
betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]

betaResultsCortex <- betaResults

### Hippocampus
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_rTg4510_HIP.RData")
betaResults <- betaResultsChip

#Convert the betaresults to the right format
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]
betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
# Remove NA rows
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
mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)
betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]

betaResultsHip <- betaResults

rm(list=setdiff(ls(), c("betaResultsCortex", "betaResultsHip"))) 


# Table
#Make a table that lists the number os significant sites and how many are shared in both tissues
betaResultsCortexSig <- betaResultsCortex[which(betaResultsCortex$FDR_adj_Interaction <= 0.05),]
betaResultsHipSig <- betaResultsHip[which(betaResultsHip$FDR_adj_Interaction <= 0.05),]

InteractionTbl <- matrix(NA, nrow = 2, ncol = 3)
InteractionTbl[1,1] <- nrow(betaResultsCortex)
InteractionTbl[1,2] <- nrow(betaResultsCortexSig)
InteractionTbl[1,3] <- length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg))
InteractionTbl[2,1] <- nrow(betaResultsHip)
InteractionTbl[2,2] <- nrow(betaResultsHipSig)
InteractionTbl[2,3] <- length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg))

InteractionTbl <- as.data.frame(InteractionTbl)
rownames(InteractionTbl) <- c("Entorhinal Cortex", "Hippocampus")
colnames(InteractionTbl) <- c("sites", "significant sites", "shared significant sites")

library(gridExtra)
png("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_InteractionCortexHip.png")
p<-tableGrob(InteractionTbl)
grid.arrange(p)
dev.off()

## Pull out the shared significant sites....
sharesigsites <- intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)
CorSharedSigSites  <- betaResultsCortex[which(betaResultsCortex$cpg %in% sharesigsites),]
HipSharedSigSites  <- betaResultsHip[which(betaResultsHip$cpg %in% sharesigsites),]

SharedSigSites <- cbind(CorSharedSigSites, HipSharedSigSites)
SharedSigSites <- SharedSigSites[,c(1,42, 2,3,17,20,38,41)]
colnames(SharedSigSites) <- c("cpg","Gene_Symbol", "Chr","Bp","Cortex.estimate","Cortex.FDR_adj_Interaction",
                              "Hip.estimate","Hip.FDR_adj_Interaction")
SharedSigSites <- SharedSigSites[order(SharedSigSites$Cortex.estimate, decreasing = T),]

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_InteractionSharedSigSites.pdf", w = 15, h = 75)
p<-tableGrob(SharedSigSites)
grid.arrange(p)
dev.off()

write.csv(SharedSigSites, "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_InteractionSharedSigSites.csv")


##Plots

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_InteractionEffectSizesCortexHip.pdf")
plot(betaResultsCortex$estimate.Interaction, betaResultsHip$estimate.Interaction, cex = 0.7,
	xlab="Entorhinal Cortex", ylab="Hippocampus", main="All Sites Effect Sizes")
abline(v=0, h=0)
dev.off()

#filter for hippocampus significant sites then plot the effect sizes again ....
betaResultsCortexSig2 <- betaResultsCortex[which(betaResultsCortex$cpg %in% betaResultsHipSig$cpg),]
identical(betaResultsCortexSig2$cpg, betaResultsHipSig$cpg)
#plot effect sizes but with hippocampus number of sig as there will be more points. Color red where sites are sig in cortex
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_InteractionEffectSizesCortexHip_HipSignificant.pdf")
plot(betaResultsCortexSig2$estimate.Interaction, betaResultsHipSig$estimate.Interaction, cex = 0.7,
	col=ifelse(betaResultsCortexSig2$estimate.Interaction %in% betaResultsCortexSig$estimate.Interaction, "red", "black"),
	 main="Hippocampus Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
dev.off()

#filter for cortex significant sites then plot the effect sizes again ....
betaResultsHipSig2 <- betaResultsHip[which(betaResultsHip$cpg %in% betaResultsCortexSig$cpg),]
identical(betaResultsHipSig2$cpg, betaResultsCortexSig$cpg)
#plot effect sizes but with cortex number of sig. Color red where sig sites are shared
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_InteractionEffectSizesCortexHip_CortexSignificant.pdf")
plot(betaResultsCortexSig$estimate.Interaction, betaResultsHipSig2$estimate.Interaction, cex = 0.7,
	col=ifelse(betaResultsHipSig2$estimate.Interaction %in% betaResultsHipSig$estimate.Interaction, "red", "black"),
	 main="Cortex Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
dev.off()


sharedsigsites <- intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)
sharedsigsites <- betaResultsCortex[which(betaResultsCortex$cpg %in% sharedsigsites), c("cpg", "Gene_Symbol", "estimate.Interaction", "FDR_adj_Interaction")]
sharedsigsites <- sharedsigsites[order(sharedsigsites$FDR_adj_Interaction),]
if(nrow(sharedsigsites) > 0){
png("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_InteractionSharedSigSites.png")
p<-tableGrob(sharedsigsites)
grid.arrange(p)
dev.off()
} else{
	print("No rows for shared sites information")
}

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_InteractionPvals.pdf", onefile = T)
plot(betaResultsCortex$FDR_adj_Interaction, betaResultsHip$FDR_adj_Interaction)
plot(-log10(betaResultsCortex$FDR_adj_Interaction),-log10(betaResultsHip$FDR_adj_Interaction))
dev.off()


# Make one figure with all effect sizes in one pfd; 4 panels a) All sites b) Cortex sig sites c) Hip sig sites d) shares sig sites
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/rTg_InteractionEffectSizes.pdf", onefile = T)
par(mfrow = c(2,2))
# a) All sites
plot(betaResultsCortex$estimate.Interaction, betaResultsHip$estimate.Interaction, cex = 0.7,
	xlab="Entorhinal Cortex", ylab="Hippocampus", main="All Sites Effect Sizes")
abline(v=0, h=0)

# b) Cortex Sig sites
plot(betaResultsCortexSig$estimate.Interaction, betaResultsHipSig2$estimate.Interaction, cex = 0.7,
	col=ifelse(betaResultsHipSig2$estimate.Interaction %in% betaResultsHipSig$estimate.Interaction, "red", "black"),
	 main="Cortex Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)

# c) Hip Sig sites
plot(betaResultsCortexSig2$estimate.Interaction, betaResultsHipSig$estimate.Interaction, cex = 0.7,
	col=ifelse(betaResultsCortexSig2$estimate.Interaction %in% betaResultsCortexSig$estimate.Interaction, "red", "black"),
	 main="Hippocampus Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)

# d) Sig only sites
if(length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)) > 0){
betaResultsCortexSig3 <- betaResultsCortex[which(betaResultsCortex$cpg %in% intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)),]
betaResultsHipSig3 <- betaResultsHip[which(betaResultsHip$cpg %in% intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)),]
plot(betaResultsCortexSig3$estimate.Interaction, betaResultsHipSig3$estimate.Interaction, cex = 0.7,
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
betaResults <- betaResultsChip
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

#Convert the betaresults to the right format
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]
betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
# Remove NA rows
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
mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)
betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]

betaResultsCortex <- betaResults

### Hippocampus
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_J20_HIP.RData")
betaResults <- betaResultsChip

#Convert the betaresults to the right format
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]
betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
# Remove NA rows
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
mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)
betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]

betaResultsHip <- betaResults

rm(list=setdiff(ls(), c("betaResultsCortex", "betaResultsHip"))) 


# Table
#Make a table that lists the number os significant sites and how many are shared in both tissues
betaResultsCortexSig <- betaResultsCortex[which(betaResultsCortex$FDR_adj_Interaction <= 0.05),]
betaResultsHipSig <- betaResultsHip[which(betaResultsHip$FDR_adj_Interaction <= 0.05),]

InteractionTbl <- matrix(NA, nrow = 2, ncol = 3)
InteractionTbl[1,1] <- nrow(betaResultsCortex)
InteractionTbl[1,2] <- nrow(betaResultsCortexSig)
InteractionTbl[1,3] <- length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg))
InteractionTbl[2,1] <- nrow(betaResultsHip)
InteractionTbl[2,2] <- nrow(betaResultsHipSig)
InteractionTbl[2,3] <- length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg))

InteractionTbl <- as.data.frame(InteractionTbl)
rownames(InteractionTbl) <- c("Entorhinal Cortex", "Hippocampus")
colnames(InteractionTbl) <- c("sites", "significant sites", "shared significant sites")

library(gridExtra)
png("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_InteractionCortexHip.png")
p<-tableGrob(InteractionTbl)
grid.arrange(p)
dev.off()


##Plots

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_InteractionEffectSizesCortexHip.pdf")
plot(betaResultsCortex$estimate.Interaction, betaResultsHip$estimate.Interaction, cex = 0.7,
	xlab="Entorhinal Cortex", ylab="Hippocampus", main="All Sites Effect Sizes")
abline(v=0, h=0)
dev.off()

#filter for hippocampus significant sites then plot the effect sizes again ....
betaResultsCortexSig2 <- betaResultsCortex[which(betaResultsCortex$cpg %in% betaResultsHipSig$cpg),]
identical(betaResultsCortexSig2$cpg, betaResultsHipSig$cpg)
#plot effect sizes but with hippocampus number of sig as there will be more points. Color red where sites are sig in cortex
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_InteractionEffectSizesCortexHip_HipSignificant.pdf")
plot(betaResultsCortexSig2$estimate.Interaction, betaResultsHipSig$estimate.Interaction, cex = 0.7,
	col=ifelse(betaResultsCortexSig2$estimate.Interaction %in% betaResultsCortexSig$estimate.Interaction, "red", "black"),
	 main="Hippocampus Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
dev.off()

#filter for cortex significant sites then plot the effect sizes again ....
betaResultsHipSig2 <- betaResultsHip[which(betaResultsHip$cpg %in% betaResultsCortexSig$cpg),]
identical(betaResultsHipSig2$cpg, betaResultsCortexSig$cpg)
#plot effect sizes but with cortex number of sig. Color red where sig sites are shared
if(nrow(betaResultsCortexSig) > 0){
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_InteractionEffectSizesCortexHip_CortexSignificant.pdf")
plot(betaResultsCortexSig$estimate.Interaction, betaResultsHipSig2$estimate.Interaction, cex = 0.7,
	col=ifelse(betaResultsHipSig2$estimate.Interaction %in% betaResultsHipSig$estimate.Interaction, "red", "black"),
	 main="Cortex Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
dev.off()
} else {
	print("No cortex significant sites so no plot for cortex")
}

sharedsigsites <- intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)
sharedsigsites <- betaResultsCortex[which(betaResultsCortex$cpg %in% sharedsigsites), c("cpg", "Gene_Symbol", "estimate.Interaction", "FDR_adj_Interaction")]
sharedsigsites <- sharedsigsites[order(sharedsigsites$FDR_adj_Interaction),]
if(nrow(sharedsigsites) > 0){
png("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_InteractionSharedSigSites.png")
p<-tableGrob(sharedsigsites)
grid.arrange(p)
dev.off()
} else{
	print("No rows for shared sites information")
}

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_InteractionPvals.pdf", onefile = T)
plot(betaResultsCortex$FDR_adj_Interaction, betaResultsHip$FDR_adj_Interaction)
plot(-log10(betaResultsCortex$FDR_adj_Interaction),-log10(betaResultsHip$FDR_adj_Interaction))
dev.off()



pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_InteractionEffectSizes.pdf", onefile = T)
par(mfrow = c(2,2))
# a) All sites
plot(betaResultsCortex$estimate.Interaction, betaResultsHip$estimate.Interaction, cex = 0.7,
	xlab="Entorhinal Cortex", ylab="Hippocampus", main="All Sites Effect Sizes")
abline(v=0, h=0)

# b) Cortex Sig sites
if(nrow(betaResultsCortexSig) > 0){
plot(betaResultsCortexSig$estimate.Interaction, betaResultsHipSig2$estimate.Interaction, cex = 0.7,
	col=ifelse(betaResultsHipSig2$estimate.Interaction %in% betaResultsHipSig$estimate.Interaction, "red", "black"),
	 main="Cortex Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
} else {
	print("No cortex significant sites so no plot for cortex")
}

# c) Hip Sig sites
plot(betaResultsCortexSig2$estimate.Interaction, betaResultsHipSig$estimate.Interaction, cex = 0.7,
	col=ifelse(betaResultsCortexSig2$estimate.Interaction %in% betaResultsCortexSig$estimate.Interaction, "red", "black"),
	 main="Hippocampus Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)

# d) Sig only sites
if(length(intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)) > 0){
betaResultsCortexSig3 <- betaResultsCortex[which(betaResultsCortex$cpg %in% intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)),]
betaResultsHipSig3 <- betaResultsHip[which(betaResultsHip$cpg %in% intersect(betaResultsCortexSig$cpg, betaResultsHipSig$cpg)),]
plot(betaResultsCortexSig3$estimate.Interaction, betaResultsHipSig3$estimate.Interaction, cex = 0.7,
	 main="Shared Significant Sites Effect Sizes",
	 xlab="Entorhinal Cortex", ylab="Hippocampus")
abline(v=0, h=0)
} else {
	print("4th plot of significant sites not plotted")
}

dev.off()