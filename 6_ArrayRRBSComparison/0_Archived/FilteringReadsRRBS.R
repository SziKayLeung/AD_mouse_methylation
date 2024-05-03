### Compare the 1201 sites in rTg4510 and plot the coverage with it
library(tidyverse)
library(BiSeq)
library(ggplot2)
library(ggpubr)
library(reshape2)

### Load data
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat") # unique probes
load(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/biseq_rTg4510.RData")


############# Data Cleaning ###############

# Pull out betas matrices in both datasets
RRBS_betas <- as.data.frame(methLevel(predictedMeth))
rrbs_betas_totalreads <- as.data.frame(totalReads(rrbs.clust.lim))
# rrbs_betas_methreads <- as.data.frame(methReads(rrbs.clust.lim))
# Object of class SimpleList of two matrices, named totalReads and methReads. 
# The matrix totalReads contains the number of reads spanning a CpG-site. The rows represent the CpG sites in rowRanges and the columns represent the samples in colData.
# The matrix methReads contains the number of methylated reads spanning a CpG-site.
## We are interested in how many reads there are so we will use total reads for the meantime

if(!identical(colnames(Normalised_Sesame_Betas), QCmetrics$Basename)) {
  print("Colnames of beta matrix and pheno basenames do not match")
}
sampleid <- gsub(".*_", "", QCmetrics$ExternalSampleID) # Pull out the sample id so that the array sample id match the rrbs data 
QCmetrics$sampleid <- sampleid
colnames(Normalised_Sesame_Betas) <- sampleid
Array_betas <- Normalised_Sesame_Betas[,QCmetrics[which(QCmetrics$Tissue == "CortexEntorhinalis" & QCmetrics$AD_model == "rTg4510"), "sampleid"]] #pull out samples that exist in the rrbs


# Pull out probes that exist in both datasets
dim(Array_betas) #24048 probes
dim(RRBS_betas) #1,761,997 probes

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
rrbs_betas_totalreads_cor <- as.data.frame(rrbs.clust.lim@rowRanges)
identical(rownames(rrbs_betas_totalreads), rownames(rrbs_betas_totalreads_cor))
rrbs_betas_totalreads_cor$Position <- paste(rrbs_betas_totalreads_cor$seqnames, rrbs_betas_totalreads_cor$start, sep = ":")

# Filter on coordinats that match then collapse the betas to match the beta matrix
sharedsites <- intersect(RRBS_coordinates$Position, Array_coordinates$Position) 
sharedsitesArray <- sharedsites
length(sharedsites) #1201
SharedArray_coordinates <- Array_coordinates[which(Array_coordinates$Position %in% sharedsites),]
SharedRRBS_coordinates <- RRBS_coordinates[which(RRBS_coordinates$Position %in% sharedsites),]
Sharedrrbs_coordiantes <- rrbs_betas_totalreads_cor[which(rrbs_betas_totalreads_cor$Position %in% sharedsites),]

# Filter these sites on the beta matrices and make sure they are in the same order
SharedArray_betas <- Array_betas[rownames(Array_betas) %in% SharedArray_coordinates$probeID,]
SharedRRBS_betas <- RRBS_betas[rownames(RRBS_betas) %in% rownames(SharedRRBS_coordinates),]
Sharedrrbs_totalreads <- rrbs_betas_totalreads[rownames(rrbs_betas_totalreads) %in% rownames(Sharedrrbs_coordiantes),]
rownames(SharedArray_betas) <- SharedArray_coordinates$Position
rownames(SharedRRBS_betas) <- SharedRRBS_coordinates$Position
rownames(Sharedrrbs_totalreads) <- Sharedrrbs_coordiantes$Position
SharedArray_betas <- SharedArray_betas[order(rownames(SharedArray_betas)),]
SharedRRBS_betas <- SharedRRBS_betas[order(rownames(SharedRRBS_betas)),]
Sharedrrbs_totalreads <- Sharedrrbs_totalreads[order(rownames(Sharedrrbs_totalreads)),]
identical(rownames(SharedRRBS_betas), rownames(SharedArray_betas))
identical(rownames(SharedRRBS_betas), rownames(Sharedrrbs_totalreads))
# Now sites match... need to select the right samples to compare
# Filter for same samples (40 samples)
sharedsamples <- intersect(colnames(SharedRRBS_betas),colnames(SharedArray_betas))
SharedArray_betas <- SharedArray_betas[,sharedsamples]
SharedRRBS_betas <- SharedRRBS_betas[,sharedsamples]
Sharedrrbs_totalreads <- Sharedrrbs_totalreads[,sharedsamples]
identical(colnames(SharedRRBS_betas), colnames(SharedArray_betas))
identical(colnames(SharedRRBS_betas), colnames(Sharedrrbs_totalreads))

save(Sharedrrbs_totalreads, file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/CompareTech/rTg4510/Reads.RData")


############## Plotting #############
# For loop to plot scatter plot for each sample with pearson correlation

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/CompareTech/rTg4510/RRBSvsArray_wTotalReads.pdf")
for(i in 1:ncol(SharedArray_betas)){
  
  comb_df <- cbind(SharedArray_betas[,i], SharedRRBS_betas[,i], Sharedrrbs_totalreads[,i])
  comb_df <- as.data.frame(comb_df)
  colnames(comb_df) <- c("Array", "RRBS", "RRBS_TotalReads")
  
  p <- ggplot(comb_df , aes(x = Array, y = RRBS, colour = RRBS_TotalReads))+
    geom_point() +
    geom_smooth(method=lm, se=FALSE) +
    theme_classic() +
    stat_cor(method="pearson") +
    labs(title = colnames(SharedArray_betas)[i], x = "Array Beta Methylation", y = "Smoothed Methylation") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_gradient(high = "blue",low = "red")
  
  q <- ggplot(comb_df , aes(RRBS_TotalReads))+
    geom_density(aes(fill = "red", alpha  =0.5)) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(ggarrange(q,p, ncol = 1, nrow = 2,
                  heights = c(0.4,0.8)))
}
dev.off()



######### Meth reads

# Pull out betas matrices in both datasets
RRBS_betas <- as.data.frame(methLevel(predictedMeth))
rrbs_betas_methreads <- as.data.frame(methReads(rrbs.clust.lim))
# rrbs_betas_methreads <- as.data.frame(methReads(rrbs.clust.lim))
# Object of class SimpleList of two matrices, named methReads and methReads. 
# The matrix methReads contains the number of reads spanning a CpG-site. The rows represent the CpG sites in rowRanges and the columns represent the samples in colData.
# The matrix methReads contains the number of methylated reads spanning a CpG-site.
## We are interested in how many reads there are so we will use meth reads for the meantime

if(!identical(colnames(Normalised_Sesame_Betas), QCmetrics$Basename)) {
  print("Colnames of beta matrix and pheno basenames do not match")
}
sampleid <- gsub(".*_", "", QCmetrics$ExternalSampleID) # Pull out the sample id so that the array sample id match the rrbs data 
QCmetrics$sampleid <- sampleid
colnames(Normalised_Sesame_Betas) <- sampleid
Array_betas <- Normalised_Sesame_Betas[,QCmetrics[which(QCmetrics$Tissue == "CortexEntorhinalis" & QCmetrics$AD_model == "rTg4510"), "sampleid"]] #pull out samples that exist in the rrbs


# Pull out probes that exist in both datasets
dim(Array_betas) #24048 probes
dim(RRBS_betas) #1,761,997 probes

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
rrbs_betas_methreads_cor <- as.data.frame(rrbs.clust.lim@rowRanges)
identical(rownames(rrbs_betas_methreads), rownames(rrbs_betas_methreads_cor))
rrbs_betas_methreads_cor$Position <- paste(rrbs_betas_methreads_cor$seqnames, rrbs_betas_methreads_cor$start, sep = ":")

# Filter on coordinats that match then collapse the betas to match the beta matrix
sharedsites <- intersect(RRBS_coordinates$Position, Array_coordinates$Position) 
sharedsitesArray <- sharedsites
length(sharedsites) #1201
SharedArray_coordinates <- Array_coordinates[which(Array_coordinates$Position %in% sharedsites),]
SharedRRBS_coordinates <- RRBS_coordinates[which(RRBS_coordinates$Position %in% sharedsites),]
Sharedrrbs_coordiantes <- rrbs_betas_methreads_cor[which(rrbs_betas_methreads_cor$Position %in% sharedsites),]

# Filter these sites on the beta matrices and make sure they are in the same order
SharedArray_betas <- Array_betas[rownames(Array_betas) %in% SharedArray_coordinates$probeID,]
SharedRRBS_betas <- RRBS_betas[rownames(RRBS_betas) %in% rownames(SharedRRBS_coordinates),]
Sharedrrbs_methreads <- rrbs_betas_methreads[rownames(rrbs_betas_methreads) %in% rownames(Sharedrrbs_coordiantes),]
rownames(SharedArray_betas) <- SharedArray_coordinates$Position
rownames(SharedRRBS_betas) <- SharedRRBS_coordinates$Position
rownames(Sharedrrbs_methreads) <- Sharedrrbs_coordiantes$Position
SharedArray_betas <- SharedArray_betas[order(rownames(SharedArray_betas)),]
SharedRRBS_betas <- SharedRRBS_betas[order(rownames(SharedRRBS_betas)),]
Sharedrrbs_methreads <- Sharedrrbs_methreads[order(rownames(Sharedrrbs_methreads)),]
identical(rownames(SharedRRBS_betas), rownames(SharedArray_betas))
identical(rownames(SharedRRBS_betas), rownames(Sharedrrbs_methreads))
# Now sites match... need to select the right samples to compare
# Filter for same samples (40 samples)
sharedsamples <- intersect(colnames(SharedRRBS_betas),colnames(SharedArray_betas))
SharedArray_betas <- SharedArray_betas[,sharedsamples]
SharedRRBS_betas <- SharedRRBS_betas[,sharedsamples]
Sharedrrbs_methreads <- Sharedrrbs_methreads[,sharedsamples]
identical(colnames(SharedRRBS_betas), colnames(SharedArray_betas))
identical(colnames(SharedRRBS_betas), colnames(Sharedrrbs_methreads))

############## Plotting #############
# For loop to plot scatter plot for each sample with pearson correlation

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/CompareTech/rTg4510/RRBSvsArray_wMethReads.pdf")
for(i in 1:ncol(SharedArray_betas)){
  
  comb_df <- cbind(SharedArray_betas[,i], SharedRRBS_betas[,i], Sharedrrbs_methreads[,i])
  comb_df <- as.data.frame(comb_df)
  colnames(comb_df) <- c("Array", "RRBS", "RRBS_MethReads")
  
  p <- ggplot(comb_df , aes(x = Array, y = RRBS, colour = RRBS_MethReads))+
    geom_point() +
    geom_smooth(method=lm, se=FALSE) +
    theme_classic() +
    stat_cor(method="pearson") +
    labs(title = colnames(SharedArray_betas)[i], x = "Array Beta Methylation", y = "Smoothed Methylation") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_gradient(high = "blue",low = "red")
  
  q <- ggplot(comb_df , aes(RRBS_MethReads))+
    geom_density(aes(fill = "red", alpha  =0.5)) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(ggarrange(q,p, ncol = 1, nrow = 2,
                  heights = c(0.4,0.8)))
}
dev.off()

resave(SharedArray_betas, SharedRRBS_betas, Sharedrrbs_methreads, file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/CompareTech/rTg4510/Reads.RData" )


########################### Load at coverage of these 1201 sites ##################

load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/CompareTech/rTg4510/Reads.RData")

## plot x = meth reads, y = methylations
ggplot(comb_df, aes(x = RRBS_MethReads, y = RRBS)) +
  geom_point()+
  theme_classic()


identical(rownames(Sharedrrbs_methreads), rownames(Sharedrrbs_totalreads))
identical(colnames(Sharedrrbs_methreads), colnames(Sharedrrbs_totalreads))

i = 1
comb_df <- cbind(SharedArray_betas[,i], SharedRRBS_betas[,i], Sharedrrbs_methreads[,i], Sharedrrbs_totalreads[,i])
comb_df <- as.data.frame(comb_df)
colnames(comb_df) <- c("Array", "RRBS", "RRBS_MethReads", "RRBS_TotalReads")


plot(comb_df$RRBS_TotalReads, comb_df$RRBS_MethReads)
ggplot(comb_df, aes(x = RRBS_TotalReads,
                    y = RRBS_MethReads, col = Array)) +
  geom_point()+
  theme_classic()
plot(comb_df$RRBS_MethReads/comb_df$RRBS_TotalReads, comb_df$RRBS,
     xlab = "MethReads/TotalReads", ylab = "RRBS Smoothed Value",
     main = paste(colnames(SharedArray_betas)[i], nrow(comb_df), sep = "  Sites: "))
corr <- cor(comb_df$RRBS_MethReads/comb_df$RRBS_TotalReads, comb_df$RRBS,use = "complete.obs")
mtext(paste("Lamda = ", signif(corr,3)), side = 3, adj = 1)

ggplot(comb_df, aes(x = RRBS_MethReads/RRBS_TotalReads,
                    y = RRBS, col = RRBS_TotalReads)) +
  geom_point()+
  theme_classic() +
  labs(title = paste(colnames(SharedArray_betas)[i], nrow(comb_df), sep = "  Sites:"), 
       x = "MethReads/TotalReads", y = "RRBS Smoothed Value", color = "Total Reads") +
  geom_smooth(method=lm, se=FALSE, col = "black") +
  stat_cor(method="pearson")


comb_df2 <- matrix(NA, nrow = nrow(SharedRRBS_betas), ncol = ncol(SharedRRBS_betas))
comb_df2 <- as.data.frame(comb_df2)
for(i in 1:ncol(SharedRRBS_betas)){
  comb_df <- cbind(SharedArray_betas[,i], SharedRRBS_betas[,i], Sharedrrbs_methreads[,i], Sharedrrbs_totalreads[,i])
  comb_df <- as.data.frame(comb_df)
  colnames(comb_df) <- c("Array", "RRBS", "RRBS_MethReads", "RRBS_TotalReads")
  comb_df <- comb_df[order(comb_df$RRBS_TotalReads, decreasing = T),]
  comb_df2[,i] <- comb_df$RRBS_TotalReads
}
colnames(comb_df2) <- colnames(SharedRRBS_betas)
rownames(comb_df2) <- rownames(SharedRRBS_betas)

plot(c(1:nrow(comb_df2)), comb_df2[,1], type="l",main = "Total reads per site",
     xlab = "Sites", ylab = "Total Reads")
for ( i in 2:ncol(comb_df2)){
  lines(c(1:nrow(comb_df2)), comb_df2[,i])
}

## filter for sites with over 10 reads.
df <- subset(comb_df2, comb_df2 > 10)
df <- na.omit(df)
range(df)


plot(c(1:nrow(df)), df[,1], type="l",main = "Total reads per site",
     xlab = "Sites", ylab = "Total Reads")
for ( i in 2:ncol(df)){
  lines(c(1:nrow(df)), df[,i])
}

## Plot correlation Array to RRBS - does the correaltion perform better?

# SharedArray_betas <- SharedArray_betas[!rownames(SharedArray_betas) %in% rownames(df),]
# SharedRRBS_betas <- SharedRRBS_betas[!rownames(SharedRRBS_betas) %in% rownames(df),]
# Sharedrrbs_totalreads <- Sharedrrbs_totalreads[!rownames(Sharedrrbs_totalreads) %in% rownames(df),]

thresholds <- c(5,10,15,20,25)

for(t in 1:length(thresholds)){
  filenames <- paste("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/CompareTech/rTg4510/FilteredRRBSvsArray_wTotalReads", thresholds[t], sep = "")
  pdf(paste(filenames, ".pdf", sep = ""))
  for(i in 1:ncol(SharedArray_betas)){
    
    comb_df <- cbind(SharedArray_betas[,i], SharedRRBS_betas[,i], Sharedrrbs_methreads[,i], Sharedrrbs_totalreads[,i])
    comb_df <- as.data.frame(comb_df)
    colnames(comb_df) <- c("Array", "RRBS", "RRBS_MethReads", "RRBS_TotalReads")
    
    #filter for rows in total reads > 10
    
    comb_df <- comb_df[which(comb_df$RRBS_TotalReads > thresholds[t]),]
  
    
    p <- ggplot(comb_df , aes(x = Array, y = RRBS, colour = RRBS_TotalReads))+
      geom_point() +
      geom_smooth(method=lm, se=FALSE) +
      theme_classic() +
      stat_cor(method="pearson") +
      labs(title = paste(colnames(SharedArray_betas)[i], nrow(comb_df), sep = "  Sites:"), x = "Array Beta Methylation", y = "Smoothed Methylation") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_colour_gradient(high = "blue",low = "red")
    
    q <- ggplot(comb_df , aes(RRBS_TotalReads))+
      geom_density(aes(fill = "red", alpha  =0.5)) +
      theme_classic() +
      theme(legend.position = "none") 
    
    print(ggarrange(q,p, ncol = 1, nrow = 2,
                    heights = c(0.4,0.8)))
  }
  dev.off()
}

i = 1
comb_df <- cbind(SharedArray_betas[,i], SharedRRBS_betas[,i], Sharedrrbs_methreads[,i], Sharedrrbs_totalreads[,i])
comb_df <- as.data.frame(comb_df)
colnames(comb_df) <- c("Array", "RRBS", "RRBS_MethReads", "RRBS_TotalReads")

#filter for rows in total reads > 10
comb_df <- comb_df[which(comb_df$RRBS_TotalReads > 10),]

ggplot(comb_df, aes(x = RRBS_MethReads/RRBS_TotalReads,
                    y = RRBS, col = RRBS_TotalReads)) +
  geom_point()+
  theme_classic() +
  labs(title = paste(colnames(SharedArray_betas)[i], nrow(comb_df), sep = "  Sites:"), 
       x = "MethReads/TotalReads", y = "RRBS Smoothed Value", color = "Total Reads") +
  geom_smooth(method=lm, se=FALSE, col = "black") +
  stat_cor(method="pearson")


################################### Pre-smoothing & post-smoothing RRBS ###############################

## pre-smoothed data

SharedRRBS_betas
Sharedrrbs_methreads