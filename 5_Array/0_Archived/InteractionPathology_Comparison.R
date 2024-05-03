## Aisha Dahir A.N.Dahir@exeter.ac.uk
## 15.07.20
## This script is to compare the interaction and pathology results - effect sizes?


library(tidyverse)
######################## Load data ########################

print("#################################### Interaction ################################")
print("#################################### rTg4510 ################################")


###  CORTEX 
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_rTg4510.RData")
betaResults <- betaResultsChip

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

betaResultsrTgInteractionCortex <- betaResultsCortex
betaResultsrTgInteractionHip <- betaResultsHip

rm(list=setdiff(ls(), c("betaResultsrTgInteractionCortex", "betaResultsrTgInteractionHip"))) 


print("#################################### Interaction ################################")
print("#################################### J20 ################################")



###  CORTEX 
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_J20.RData")
betaResults <- betaResultsChip

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

betaResultsJ20InteractionCortex <- betaResultsCortex
betaResultsJ20InteractionHip <- betaResultsHip

rm(list=setdiff(ls(), c("betaResultsrTgInteractionCortex", "betaResultsrTgInteractionHip",
						"betaResultsJ20InteractionCortex", "betaResultsJ20InteractionHip"))) 


print("#################################### Pathology ################################")
print("#################################### rTg4510 ################################")

###  CORTEX 
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_rTg4510.RData")
betaResults <- betaResultsPathology

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

betaResultsrTgPathologyCortex <- betaResultsCortex
betaResultsrTgPathologyHip <- betaResultsHip

rm(list=setdiff(ls(), c("betaResultsrTgInteractionCortex", "betaResultsrTgInteractionHip",
						"betaResultsJ20InteractionCortex", "betaResultsJ20InteractionHip", 
						"betaResultsrTgPathologyCortex", "betaResultsrTgPathologyHip"))) 
 

print("#################################### Pathology ################################")
print("#################################### J20 ################################")

###  CORTEX 
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_J20.RData")
betaResults <- betaResultsPathology

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

betaResultsJ20PathologyCortex <- betaResultsCortex
betaResultsJ20PathologyHip <- betaResultsHip

rm(list=setdiff(ls(), c("betaResultsrTgInteractionCortex", "betaResultsrTgInteractionHip",
						"betaResultsJ20InteractionCortex", "betaResultsJ20InteractionHip", 
						"betaResultsrTgPathologyCortex", "betaResultsrTgPathologyHip",
						"betaResultsJ20PathologyCortex", "betaResultsJ20PathologyHip")))


########################   Plots   ########################



pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Pvals_InteractionPathology.pdf")
par(mfrow = c(2,2))

# a) Cortex rTg
plot(-log10(betaResultsrTgInteractionCortex$p.val.Interaction),
     -log10(betaResultsrTgPathologyCortex$p.val.Pathology),
     xlab = "Interaction", ylab = "Pathology",
     main = "rTg Cortex", cex = 0.7)

# b) Cortex J20
plot(-log10(betaResultsJ20InteractionCortex$p.val.Interaction), 
     -log10(betaResultsJ20PathologyCortex$p.val.Pathology),
     xlab = "Interaction", ylab = "Pathology",
     main = "J20 Cortex", cex = 0.7)

# c) Hip rTg
plot(-log10(betaResultsrTgInteractionHip$p.val.Interaction),
     -log10(betaResultsrTgPathologyHip$p.val.Pathology),
     xlab = "Interaction", ylab = "Pathology",
     main = "rTg Hippocampus", cex = 0.7)

# d) Hip J20
plot(-log10(betaResultsJ20InteractionHip$p.val.Interaction),
     -log10(betaResultsJ20PathologyHip$p.val.Pathology),
     xlab = "Interaction", ylab = "Pathology",
     main = "J20 Hippocampus", cex = 0.7)

dev.off()



##################### Interaction effect size = Interaction + Age
pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/EffectSizes_InteractionPathology.pdf")
par(mfrow = c(2,2))
# a) Cortex rTg
betaResultsrTgInteractionCortexSig <- betaResultsrTgInteractionCortex[which(betaResultsrTgInteractionCortex$FDR_adj_Interaction < 0.05),]
nrow(betaResultsrTgInteractionCortexSig)
betaResultsrTgPathologyCortexSig <- betaResultsrTgPathologyCortex[which(betaResultsrTgPathologyCortex$FDR_adj_Pathology < 0.05),]
nrow(betaResultsrTgPathologyCortexSig)

#pull out common sites that exist in pathology and interaction
betaResultsrTgPathologyCortexSig3 <- betaResultsrTgPathologyCortex[which(betaResultsrTgPathologyCortex$cpg %in% 
                                                                           intersect(betaResultsrTgPathologyCortexSig$cpg, betaResultsrTgInteractionCortexSig$cpg)),]
betaResultsrTgInteractionCortexSig3 <- betaResultsrTgInteractionCortex[which(betaResultsrTgInteractionCortex$cpg %in% 
                                                                           intersect(betaResultsrTgPathologyCortexSig$cpg, betaResultsrTgInteractionCortexSig$cpg)),]

#TGAge <- betaResultsrTgInteractionCortexSig3$estimate.Genotype + betaResultsrTgInteractionCortexSig3$estimate.Age
IntAge <- betaResultsrTgInteractionCortexSig3$estimate.Age + betaResultsrTgInteractionCortexSig3$estimate.Interaction
plot(IntAge,
     betaResultsrTgPathologyCortexSig3$estimate.Pathology, 
     xlab = "Interaction + Age", ylab = "Pathology",
     main = "rTg Cortex", cex = 0.7)
abline(h = 0)
abline(v = 0)


# b) Cortex J20
betaResultsJ20InteractionCortexSig <- betaResultsJ20InteractionCortex[which(betaResultsJ20InteractionCortex$FDR_adj_Interaction < 0.05),]
nrow(betaResultsJ20InteractionCortexSig)
betaResultsJ20PathologyCortexSig <- betaResultsJ20PathologyCortex[which(betaResultsJ20PathologyCortex$FDR_adj_Pathology < 0.05),]
nrow(betaResultsJ20PathologyCortexSig)

betaResultsJ20PathologyCortexSig3 <- betaResultsJ20PathologyCortex[which(betaResultsJ20PathologyCortex$cpg %in% 
                                                                           intersect(betaResultsJ20PathologyCortexSig$cpg, betaResultsJ20InteractionCortexSig$cpg)),]
betaResultsJ20InteractionCortexSig3 <- betaResultsJ20InteractionCortex[which(betaResultsJ20InteractionCortex$cpg %in% 
                                                                               intersect(betaResultsJ20PathologyCortexSig$cpg, betaResultsJ20InteractionCortexSig$cpg)),]

#TGAge <- betaResultsJ20InteractionCortexSig3$estimate.Genotype + betaResultsJ20InteractionCortexSig3$estimate.Age
IntAge <- betaResultsJ20InteractionCortexSig3$estimate.Age + betaResultsJ20InteractionCortexSig3$estimate.Interaction
plot(IntAge,
     betaResultsJ20PathologyCortexSig3$estimate.Pathology, 
     xlab = "Interaction + Age", ylab = "Pathology",
     main = "J20 Cortex", cex = 0.7)
abline(h = 0)
abline(v = 0)
# c) Hip rTg
betaResultsrTgInteractionHipSig <- betaResultsrTgInteractionHip[which(betaResultsrTgInteractionHip$FDR_adj_Interaction < 0.05),]
nrow(betaResultsrTgInteractionHipSig)
betaResultsrTgPathologyHipSig <- betaResultsrTgPathologyHip[which(betaResultsrTgPathologyHip$FDR_adj_Pathology < 0.05),]
nrow(betaResultsrTgPathologyHipSig)

#pull out common sites that exist in pathology and interaction
betaResultsrTgPathologyHipSig3 <- betaResultsrTgPathologyHip[which(betaResultsrTgPathologyHip$cpg %in% 
                                                                           intersect(betaResultsrTgPathologyHipSig$cpg, betaResultsrTgInteractionHipSig$cpg)),]
betaResultsrTgInteractionHipSig3 <- betaResultsrTgInteractionHip[which(betaResultsrTgInteractionHip$cpg %in% 
                                                                               intersect(betaResultsrTgPathologyHipSig$cpg, betaResultsrTgInteractionHipSig$cpg)),]

#TGAge <- betaResultsrTgInteractionHipSig3$estimate.Genotype + betaResultsrTgInteractionHipSig3$estimate.Age
IntAge <- betaResultsrTgInteractionHipSig3$estimate.Age + betaResultsrTgInteractionHipSig3$estimate.Interaction
plot(IntAge,
     betaResultsrTgPathologyHipSig3$estimate.Pathology, 
     xlab = "Interaction + Age", ylab = "Pathology",
     main = "rTg Hip", cex = 0.7)
abline(h = 0)
abline(v = 0)

# d) Hip J20

betaResultsJ20InteractionHipSig <- betaResultsJ20InteractionHip[which(betaResultsJ20InteractionHip$FDR_adj_Interaction < 0.05),]
nrow(betaResultsJ20InteractionHipSig)
betaResultsJ20PathologyHipSig <- betaResultsJ20PathologyHip[which(betaResultsJ20PathologyHip$FDR_adj_Pathology < 0.05),]
nrow(betaResultsJ20PathologyHipSig)

#pull out common sites that exist in pathology and interaction
betaResultsJ20PathologyHipSig3 <- betaResultsJ20PathologyHip[which(betaResultsJ20PathologyHip$cpg %in% 
                                                                     intersect(betaResultsJ20PathologyHipSig$cpg, betaResultsJ20InteractionHipSig$cpg)),]
betaResultsJ20InteractionHipSig3 <- betaResultsJ20InteractionHip[which(betaResultsJ20InteractionHip$cpg %in% 
                                                                         intersect(betaResultsJ20PathologyHipSig$cpg, betaResultsJ20InteractionHipSig$cpg)),]

#TGAge <- betaResultsJ20InteractionHipSig3$estimate.Genotype + betaResultsJ20InteractionHipSig3$estimate.Age
IntAge <- betaResultsJ20InteractionHipSig3$estimate.Age + betaResultsJ20InteractionHipSig3$estimate.Interaction
plot(IntAge,
     betaResultsJ20PathologyHipSig3$estimate.Pathology, 
     xlab = "Interaction + Age", ylab = "Pathology",
     main = "J20 Hip", cex = 0.7)
abline(h = 0)
abline(v = 0)

dev.off()


