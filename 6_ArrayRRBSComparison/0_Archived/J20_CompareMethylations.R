## Aisha Dahir A.N.Dahir@exeter.ac.uk

##  J20 samples only!!

## Three parts to the script :
## 1) Correlation between smoothed RRBS and array
## 2) Correlation between non-smoothed RRBS and array
## 3) Comparison of correlations

setwd("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/CompareTech/J20/")
#######################################################################################################
####################################       1)Smoothed RRBS                #############################
#######################################################################################################
# Load the  library & data in
library(ggplot2)
library(tidyverse)
library(reshape2)
library(BiSeq)
library(ggpubr)
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat") # unique probes
load("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/biseq_predictedMeth_J20.RData")


############# Data Cleaning ###############

# Pull out betas matrices in both datasets
RRBS_betas <- as.data.frame(methLevel(predictedMeth))

if(!identical(colnames(Normalised_Sesame_Betas), QCmetrics$Basename)) {
  print("Colnames of beta matrix and pheno basenames do not match")
}
sampleid <- gsub(".*_", "", QCmetrics$ExternalSampleID) # Pull out the sample id so that the array sample id match the rrbs data 
QCmetrics$sampleid <- sampleid
colnames(Normalised_Sesame_Betas) <- sampleid
Array_betas <- Normalised_Sesame_Betas[,QCmetrics[which(QCmetrics$Tissue == "CortexEntorhinalis" & QCmetrics$AD_model == "J20"), "sampleid"]] #pull out samples that exist in the rrbs


# Pull out probes that exist in both datasets
dim(Array_betas) #23633 probes
dim(RRBS_betas) #1,809,190 probes

# Pull out the genetic location in both datasets
# Array
# No genetic location provided in the array rdat so will have to get from the hovath csv file
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #24048 probes 
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- species_probes$MusMusculus
species_probes$Bp <- gsub("*.:","", species_probes$Bp)
species_probes <- species_probes[,-2]
species_probes$Position <- paste("chr", paste(species_probes$Chr, species_probes$Bp, sep = ":"), sep ="")
Array_coordinates <- species_probes

#RRBS
RRBS_coordinates <-  as.data.frame(predictedMeth@rowRanges)
RRBS_coordinates$Position <- paste(RRBS_coordinates$seqnames, RRBS_coordinates$start, sep = ":")


# Filter on coordinats that match then collapse the betas to match the beta matrix
sharedsites <- intersect(RRBS_coordinates$Position, Array_coordinates$Position) 
sharedsitesArray <- sharedsites
length(sharedsites) #1233
SharedArray_coordinates <- Array_coordinates[which(Array_coordinates$Position %in% sharedsites),]
SharedRRBS_coordinates <- RRBS_coordinates[which(RRBS_coordinates$Position %in% sharedsites),]
if(anyDuplicated(SharedRRBS_coordinates$Position) != 0){
  print("Duplicate sites in RRBS")
}
if(anyDuplicated(SharedArray_coordinates$Position) != 0){
  print("Duplicate sites in Array")
}

# NOTE - THERE ARE SOME PROBES ON THE ARRAY THAT ARE DUPLICATES - Remove one duplicate
SharedArray_coordinates$Duplicates <- duplicated(SharedArray_coordinates$Position)
SharedArray_coordinates <- SharedArray_coordinates[which(SharedArray_coordinates$Duplicates == FALSE),] # Remove one duplicate

# Filter these sites on the beta matrices and make sure they are in the same order
SharedArray_betas <- Array_betas[rownames(Array_betas) %in% SharedArray_coordinates$probeID,]
SharedRRBS_betas <- RRBS_betas[rownames(RRBS_betas) %in% rownames(SharedRRBS_coordinates),]
rownames(SharedArray_betas) <- SharedArray_coordinates$Position
rownames(SharedRRBS_betas) <- SharedRRBS_coordinates$Position
SharedArray_betas <- SharedArray_betas[order(rownames(SharedArray_betas)),]
SharedRRBS_betas <- SharedRRBS_betas[order(rownames(SharedRRBS_betas)),]
identical(rownames(SharedRRBS_betas), rownames(SharedArray_betas))

# Now sites match... need to select the right samples to compare
# Filter for same samples (40 samples)
sharedsamples <- intersect(colnames(SharedRRBS_betas),colnames(SharedArray_betas))
SharedArray_betas <- SharedArray_betas[,sharedsamples]
SharedRRBS_betas <- SharedRRBS_betas[,sharedsamples]
identical(colnames(SharedRRBS_betas), colnames(SharedArray_betas))


############## Plotting #############

pdf("UniqueProbes_ComparingSmoothedMethylationToArray.pdf",w = 12)

# Make empty vector so that in the loop the correlation (of sites) for each samples is stored in - should have 40 after loop
cor_matrix <- c()
for(i in 1:ncol(SharedArray_betas)){
  cor_matrix <- cbind(cor_matrix, cor(SharedArray_betas[,i], SharedRRBS_betas[,i], use = "complete.obs"))
}

# Make vector into data frame to plot using ggplot
cor_matrix <- as.data.frame(cor_matrix)
colnames(cor_matrix) <- colnames(SharedArray_betas)
cor_matrix <- t(cor_matrix)
cor_matrix <- as.data.frame(cor_matrix)
colnames(cor_matrix) <- "Correlation"
cor_matrix$Samples <- rownames(cor_matrix)

# Plot 1 - correlation dot plot for each sample
ggplot(cor_matrix, aes( x = cor_matrix$Samples,y= cor_matrix$Correlation, color = )) +
  geom_point() +
  theme_classic() +
  labs(x = "rTg Samples", y = "Correlation", 
       title = paste("Correlation of Smoothed Methylation and Array",  paste(length(sharedsites), " sites", sep = ""), sep = " - ")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(c(0.5,1))


# Plot 2 - Heatmap of sample correlations
corr <- cor(SharedArray_betas, SharedRRBS_betas, use = "complete.obs")
# change legenf
melted_corr <- melt(corr, na.rm = TRUE)
ggplot(data = melted_corr, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.75, limit = c(0.5,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  labs(x = "", y = "") 

# Plot 3 - For loop to plot scatter plot for each sample with pearson correlation

for(i in 1:ncol(SharedArray_betas)){
  
  comb_df <- cbind(SharedArray_betas[,i], SharedRRBS_betas[,i])
  comb_df <- as.data.frame(comb_df)
  colnames(comb_df) <- c("Array", "RRBS")
  
  p <- ggplot(comb_df , aes(x = Array, y = RRBS))+
    geom_point() +
    geom_smooth(method=lm, se=FALSE) +
    theme_classic() +
    stat_cor(method="pearson") +
    labs(title = colnames(SharedArray_betas)[i], x = "Array Beta Methylation", y = "Smoothed Methylation") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  print(p)
}

dev.off()


write.csv(cor_matrix, file = "smoothed_correlation.csv")
#######################################################################################################
#######################################        2)Non-smoothed       ###################################
#######################################################################################################
# Load the  library & data in
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
load("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/biseq_J20.RData")
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat")

############# Data Cleaning ################

# Beta matrices
# Pull out betas in both datasets - USE RRBS.CLUST.LIM AS IT LIMITS COVERAGE WITH 90% QUANTILE - THIS IS USED PRE SMOOTHING
rrbs.clust.unlim.rel <- rawToRel(rrbs.clust.unlim)
RRBS_betas <- as.data.frame(methLevel(rrbs.clust.unlim.rel))
identical(colnames(Normalised_Sesame_Betas), QCmetrics$Basename)
sampleid <- gsub(".*_", "", QCmetrics$ExternalSampleID)
QCmetrics$sampleid <- sampleid
colnames(Normalised_Sesame_Betas) <- sampleid
Array_betas <- Normalised_Sesame_Betas[,QCmetrics[which(QCmetrics$Tissue == "CortexEntorhinalis" & QCmetrics$AD_model == "J20"), "sampleid"]]


# Pull out probes that exist in both datasets
dim(Array_betas) #24048 probes
dim(RRBS_betas) #1,761,997 probes

# Array coordinates
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #24048 probes 
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- species_probes$MusMusculus
species_probes$Bp <- gsub("*.:","", species_probes$Bp)
species_probes <- species_probes[,-2]
species_probes$Position <- paste("chr", paste(species_probes$Chr, species_probes$Bp, sep = ":"), sep ="")
Array_coordinates <- species_probes

#RRBS
RRBS_coordinates <-  as.data.frame(rrbs.clust.unlim.rel@rowRanges)
RRBS_coordinates$Position <- paste(RRBS_coordinates$seqnames, RRBS_coordinates$start, sep = ":")


# Filter on coordinats that match then collapse the betas to match the beta matrix
sharedsites <- sharedsitesArray 
length(sharedsites) #1201
SharedArray_coordinates <- Array_coordinates[which(Array_coordinates$Position %in% sharedsites),]
SharedRRBS_coordinates <- RRBS_coordinates[which(RRBS_coordinates$Position %in% sharedsites),]
if(anyDuplicated(SharedRRBS_coordinates$Position) != 0){
  print("Duplicate sites in RRBS")
}
if(anyDuplicated(SharedArray_coordinates$Position) != 0){
  print("Duplicate sites in Array")
}


# Filter these sites on the beta matrices and make sure they are in the same order
SharedArray_betas <- Array_betas[rownames(Array_betas) %in% SharedArray_coordinates$probeID,]
SharedRRBS_betas <- RRBS_betas[rownames(RRBS_betas) %in% rownames(SharedRRBS_coordinates),]
rownames(SharedArray_betas) <- SharedArray_coordinates$Position
rownames(SharedRRBS_betas) <- SharedRRBS_coordinates$Position
SharedArray_betas <- SharedArray_betas[order(rownames(SharedArray_betas)),]
SharedRRBS_betas <- SharedRRBS_betas[order(rownames(SharedRRBS_betas)),]
identical(rownames(SharedRRBS_betas), rownames(SharedArray_betas))

# Now sites match... need to select the right samples to compare
# Filter for same samples (40 samples)
sharedsamples <- intersect(colnames(SharedRRBS_betas),colnames(SharedArray_betas))
SharedArray_betas <- SharedArray_betas[,sharedsamples]
SharedRRBS_betas <- SharedRRBS_betas[,sharedsamples]
identical(colnames(SharedRRBS_betas), colnames(SharedArray_betas))


############## Plotting #############

pdf("UniqueProbes_ComparingNonSmoothedMethylationToArray.pdf",w = 12)

# Make empty vector so that in the loop the correlation (of sites) for each samples is stored in - should have 40 after loop
cor_matrix <- c()
for(i in 1:ncol(SharedArray_betas)){
  cor_matrix <- cbind(cor_matrix, cor(SharedArray_betas[,i], SharedRRBS_betas[,i], use = "complete.obs"))
}

# Make vector into data frame to plot using ggplot
cor_matrix <- as.data.frame(cor_matrix)
colnames(cor_matrix) <- colnames(SharedArray_betas)
cor_matrix <- t(cor_matrix)
cor_matrix <- as.data.frame(cor_matrix)
colnames(cor_matrix) <- "Correlation"
cor_matrix$Samples <- rownames(cor_matrix)

# Plot 1 - correlation dot plot for each sample
ggplot(cor_matrix, aes( x = cor_matrix$Samples,y= cor_matrix$Correlation, color = )) +
  geom_point() +
  theme_classic() +
  labs(x = "rTg Samples", y = "Correlation", 
       title = paste("Correlation of Non Smoothed Methylation and Array",  paste(length(sharedsites), " sites", sep = ""), sep = " - ")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(c(0,1))


# Plot 2 - Heatmap of sample correlations
corr <- cor(SharedArray_betas, SharedRRBS_betas, use = "complete.obs")
# change legenf
melted_corr <- melt(corr, na.rm = TRUE)
ggplot(data = melted_corr, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  labs(x = "", y = "") 

# Plot 3 - For loop to plot scatter plot for each sample with pearson correlation

for(i in 1:ncol(SharedArray_betas)){
  
  comb_df <- cbind(SharedArray_betas[,i], SharedRRBS_betas[,i])
  comb_df <- as.data.frame(comb_df)
  colnames(comb_df) <- c("Array", "RRBS")
  
  p <- ggplot(comb_df , aes(x = Array, y = RRBS))+
    geom_point() +
    geom_smooth(method=lm, se=FALSE) +
    theme_classic() +
    stat_cor(method="pearson") +
    labs(title = colnames(SharedArray_betas)[i], x = "Array Beta Methylation", y = "Smoothed Methylation") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  print(p)
}

dev.off()

write.csv(cor_matrix, file = "nonsmoothed_correlation.csv")


#######################################################################################################
#########################           3) Compare the correlations           #############################
#######################################################################################################
nonsmoothedcorr <- read.csv("nonsmoothed_correlation.csv", header= T, row.names = 1)
smoothedcorr <- read.csv("smoothed_correlation.csv", header= T, row.names = 1)

CorrelationMatrices <- cbind(smoothedcorr, nonsmoothedcorr[,1])
colnames(CorrelationMatrices) <- c("Smoothed", "SampleID", "NonSmoothed")
CorrelationMatrices <- CorrelationMatrices[,c(2,1,3)]

ggplot(CorrelationMatrices, aes( x = CorrelationMatrices$SampleID, y= cor_matrix$Correlation, color = )) +
  geom_point() +
  theme_classic() +
  labs(x = "J20 Samples", y = "Correlation", 
       title = paste("Correlation of Non Smoothed Methylation and Array",  paste(length(sharedsites), " sites", sep = ""), sep = " - ")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(c(0,1))


pdf("Comparison_Smoothed_NonSmoothed.pdf", w = 12)
ggplot(CorrelationMatrices, aes(x = SampleID)) +
  geom_point(aes(y= Smoothed, color= "#00BFC4")) +
  geom_point(aes(y= NonSmoothed, color= "#F8766D")) +
  scale_color_manual(name="RRBS",
                     labels=c("Smoothed","Non Smoothed"),
                     values=c("#00BFC4","#F8766D")) +
  theme_classic() +
  labs(x = "J20 Samples", y = "Correlation", 
       title = "1233 sites", " sites", sep = "")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(c(0.5,1))
  
ggplot(CorrelationMatrices , aes(x = Smoothed, y = NonSmoothed))+
  geom_point() +
  geom_smooth(method=lm, se=FALSE) +
  theme_classic() +
  stat_cor(method="pearson") +
  labs(title = "1233 sites", x = "Smoothed", y = "Non-Smoothed") +
  theme(plot.title = element_text(hjust = 0.5))  
  
dev.off()
  

  