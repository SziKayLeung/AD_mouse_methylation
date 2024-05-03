## This script is to compare the results from the  pathology model

library(tidyverse)

################ Files needed

# Add annotation files - chr and bp
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_mappability_name = "MusMusculus"
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #23633
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- gsub(".*:","",species_probes$MusMusculus)

#Gene names
mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)



### Cleaning up betaResults into a function
cleantable <- function(betaResults, species_probes,mm10_Manifest){
  # Add chr and bp before deleting rows
  Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
  Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]
  betaResults <- cbind(Bp, betaResults)
  betaResults <- cbind(Chr, betaResults)
  betaResults <- betaResults[,c(10,1:9)]
  # Add FDR pvalues
  FDR_adj_pathology <- p.adjust(betaResults[,"p.val.Pathology"], method = "fdr")
  betaResults <- cbind(betaResults, FDR_adj_pathology)
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
  betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]
  
  return(betaResults)
}
pathplots <- function(betaResults, Pathology_col, col){
  betaResults <- betaResults[order(betaResults$FDR_adj_pathology),]
  rownames(betaResults) <- betaResults$cpg
  pheno$col <-  ifelse(pheno$Genotype == "WT","#000000", col )
  
  if(!identical(colnames(betas), pheno$Basename)){
    print("betas and pheno not in the same order!")
  }
  
  par(mfrow = c(2,3)) ## will plot 6 figures per page in 2 rows by 3 columns
  for(i in 1:6){
    plot(Pathology_col, as.numeric(betas[rownames(betaResults)[i],]*100), 
         xlab = "Pathology", ylab = "DNA methylation (%)", pch = 16, 
         main = paste(rownames(betas)[i], betaResults[i,"Gene_Symbol"], sep = " "),
         col = pheno$col) ## we multiple by 100 to plot of a % scale rather than proportion
    ## to add line of best fit
    mtext(paste("P = ", signif(betaResults[i, "FDR_adj_pathology"],3)), side = 3, adj = 0, cex = .7) ## add P value to plot rounded to 3 sf
    par(xpd  = TRUE)
    legend(x ="topright",legend = c("WT","TG"), 
           col = c("#000000",col), pch = 19, cex = .7,
           inset=c(-0.02,-0.25))
  }
  
}

pdf("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyModels.pdf", onefile = T, w = 10)

#############################   rTg4510   ECX   #########################
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_ECX_rTg4510.RData")

#load data
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat")
# Filter for rTg samples
pheno <- QCmetrics[which(QCmetrics$Tissue == "CortexEntorhinalis" & QCmetrics$AD_model == "rTg4510"),]
betas <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno$Basename ]
pheno$AgeMonths = ifelse(pheno$AgeDays < 100, 2, ifelse(pheno$AgeDays < 150, 4, ifelse(pheno$AgeDays < 200, 6, 8)))
pheno$Genotype <- as.factor(pheno$Genotype)
pheno$Sample_ID_ECX <- gsub(".*_", "", pheno$ExternalSampleID)
pathology <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Tg4510_coldata_VertebrateArray.csv",
                      header = T, stringsAsFactors = F)
phenowPath <- merge(pheno,pathology, by = "Sample_ID_ECX")
colnames(phenowPath)[colnames(phenowPath) == 'Genotype.x'] <- 'Genotype'
pheno <- phenowPath
betas <- betas[,order(colnames(betas))]
pheno <- pheno[order(pheno$Basename),]
identical(colnames(betas), pheno$Basename)

### model 1 : Pathology
ResultsPathologyCleaned <- cleantable(betaResults = ResultsPathology, species_probes,mm10_Manifest)

### model 2 : Pathology+Age
ResultsPathologyAgeCleaned <- cleantable(betaResults = ResultsPathologyAge, species_probes,mm10_Manifest)

### model 3 : Pathology+Age+Genotyp
ResultsPathologyAgeGenotypeCleaned <- cleantable(betaResults = ResultsPathologyAgeGenotype, species_probes,mm10_Manifest)

## Tabulate number of significant hits
PathTbl <- c(nrow(ResultsPathologyCleaned[which(ResultsPathologyCleaned$FDR_adj_pathology <= 0.05),]),
             nrow(ResultsPathologyAgeCleaned[which(ResultsPathologyAgeCleaned$FDR_adj_pathology <= 0.05),]),
             nrow(ResultsPathologyAgeGenotypeCleaned[which(ResultsPathologyAgeGenotypeCleaned$FDR_adj_pathology <= 0.05),]))
PathTbl <- as.data.frame(t(PathTbl))
colnames(PathTbl) <- c("Pathology Only", "Pathology + Age", "Pathology + Age + Genotype")
p<-tableGrob(PathTbl)
grid.arrange(top="rTg4510 ECX", p)



pathplots(betaResults = ResultsPathologyCleaned, Pathology_col = pheno$Pathology_ECX, col = "#00AEC9")
pathplots(betaResults = ResultsPathologyAgeCleaned, Pathology_col = pheno$Pathology_ECX, col = "#00AEC9")
pathplots(betaResults = ResultsPathologyAgeGenotypeCleaned, Pathology_col = pheno$Pathology_ECX, col = "#00AEC9")

## how many in path+age+geno in pathology only

ResultsPathologyAgeGenotypeCleanedSig <- ResultsPathologyAgeGenotypeCleaned[which(ResultsPathologyAgeGenotypeCleaned$FDR_adj_pathology <= 0.05),]
ResultsPathologyCleanedSig <- ResultsPathologyCleaned[which(ResultsPathologyCleaned$FDR_adj_pathology <= 0.05),]
length(intersect(ResultsPathologyAgeGenotypeCleanedSig$cpg, ResultsPathologyCleanedSig$cpg)) #68
sharedsites <- intersect(ResultsPathologyAgeGenotypeCleanedSig$cpg, ResultsPathologyCleanedSig$cpg) 

sharedsitespathfull<- ResultsPathologyAgeGenotypeCleanedSig[which(ResultsPathologyAgeGenotypeCleanedSig$cpg %in% sharedsites),] 
sharedsitespathonly<- ResultsPathologyCleanedSig[which(ResultsPathologyCleanedSig$cpg %in% sharedsites),] 
identical(sharedsitespathonly$cpg, sharedsitespathfull$cpg)

par(mfrow = c(1,2))
compPvals <- cbind(sharedsitespathonly$FDR_adj_pathology, sharedsitespathfull$FDR_adj_pathology)
compPvals <- as.data.frame(compPvals)
colnames(compPvals) <- c("Pathology only", "Pathology + Age + Genotype")
plot(x = -log(compPvals$`Pathology only`), y = -log(compPvals$`Pathology + Age + Genotype`),
     xlab= "Pathology only", ylab = "Pathology + Age + Genotype", 
     main = "rTg4510 ECX -log10 FDR Pvalues")

compEffSize <- cbind(sharedsitespathonly$estimate.Pathology, sharedsitespathfull$estimate.Pathology)
compEffSize <- as.data.frame(compEffSize)
colnames(compEffSize) <- c("Pathology only", "Pathology + Age + Genotype")
plot(x = compEffSize$`Pathology only`, y = compEffSize$`Pathology + Age + Genotype`,
     xlab= "Pathology only", ylab = "Pathology + Age + Genotype", 
     main = "rTg4510 ECX Effect sizes")
abline(h = 0)
abline(v = 0)


#############################   rTg4510   HIP  #########################
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_HIP_rTg4510.RData")

#load data
pheno <- QCmetrics[which(QCmetrics$Tissue == "Hippocampus" & QCmetrics$AD_model == "rTg4510"),]
betas <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno$Basename ]
pheno$AgeMonths = ifelse(pheno$AgeDays < 100, 2, ifelse(pheno$AgeDays < 150, 4, ifelse(pheno$AgeDays < 200, 6, 8)))
pheno$Genotype <- as.factor(pheno$Genotype)
pathology_pheno <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Tg4510_coldata_VertebrateArray.csv", header = T, stringsAsFactors = F)
pheno$Sample_ID_HIP <- gsub(".*_", "", pheno$ExternalSampleID)
phenowPath <- merge(pheno,pathology_pheno, by = "Sample_ID_HIP")
colnames(phenowPath)[colnames(phenowPath) == 'Genotype.x'] <- 'Genotype'
pheno <- phenowPath
betas <- betas[,order(colnames(betas))]
pheno <- pheno[order(pheno$Basename),]
identical(colnames(betas), pheno$Basename)


### model 1 : Pathology
ResultsPathologyCleaned <- cleantable(betaResults = ResultsPathology, species_probes,mm10_Manifest)

### model 2 : Pathology+Age
ResultsPathologyAgeCleaned <- cleantable(betaResults = ResultsPathologyAge, species_probes,mm10_Manifest)

### model 3 : Pathology+Age+Genotyp
ResultsPathologyAgeGenotypeCleaned <- cleantable(betaResults = ResultsPathologyAgeGenotype, species_probes,mm10_Manifest)

## Tabulate number of significant hits
PathTbl <- c(nrow(ResultsPathologyCleaned[which(ResultsPathologyCleaned$FDR_adj_pathology <= 0.05),]),
             nrow(ResultsPathologyAgeCleaned[which(ResultsPathologyAgeCleaned$FDR_adj_pathology <= 0.05),]),
             nrow(ResultsPathologyAgeGenotypeCleaned[which(ResultsPathologyAgeGenotypeCleaned$FDR_adj_pathology <= 0.05),]))
PathTbl <- as.data.frame(t(PathTbl))
colnames(PathTbl) <- c("Pathology Only", "Pathology + Age", "Pathology + Age + Genotype")
p<-tableGrob(PathTbl)
grid.arrange(top="rTg4510 HIP",p)

pathplots(betaResults = ResultsPathologyCleaned, Pathology_col = pheno$Pathology_HIP, col = "#00AEC9")
pathplots(betaResults = ResultsPathologyAgeCleaned, Pathology_col = pheno$Pathology_HIP, col = "#00AEC9")
pathplots(betaResults = ResultsPathologyAgeGenotypeCleaned, Pathology_col = pheno$Pathology_HIP, col = "#00AEC9")

## how many in path+age+geno in pathology only
ResultsPathologyAgeGenotypeCleanedSig <- ResultsPathologyAgeGenotypeCleaned[which(ResultsPathologyAgeGenotypeCleaned$FDR_adj_pathology <= 0.05),]
ResultsPathologyCleanedSig <- ResultsPathologyCleaned[which(ResultsPathologyCleaned$FDR_adj_pathology <= 0.05),]
length(intersect(ResultsPathologyAgeGenotypeCleanedSig$cpg, ResultsPathologyCleanedSig$cpg)) #3476
sharedsites <- intersect(ResultsPathologyAgeGenotypeCleanedSig$cpg, ResultsPathologyCleanedSig$cpg) 

sharedsitespathfull<- ResultsPathologyAgeGenotypeCleanedSig[which(ResultsPathologyAgeGenotypeCleanedSig$cpg %in% sharedsites),] 
sharedsitespathonly<- ResultsPathologyCleanedSig[which(ResultsPathologyCleanedSig$cpg %in% sharedsites),] 
identical(sharedsitespathonly$cpg, sharedsitespathfull$cpg)

par(mfrow = c(1,2))
compPvals <- cbind(sharedsitespathonly$FDR_adj_pathology, sharedsitespathfull$FDR_adj_pathology)
compPvals <- as.data.frame(compPvals)
colnames(compPvals) <- c("Pathology only", "Pathology + Age + Genotype")
plot(x = -log(compPvals$`Pathology only`), y = -log(compPvals$`Pathology + Age + Genotype`),
     xlab= "Pathology only", ylab = "Pathology + Age + Genotype", 
     main = "rTg4510 Hip -log10 FDR Pvalues")

compEffSize <- cbind(sharedsitespathonly$estimate.Pathology, sharedsitespathfull$estimate.Pathology)
compEffSize <- as.data.frame(compEffSize)
colnames(compEffSize) <- c("Pathology only", "Pathology + Age + Genotype")
plot(x = compEffSize$`Pathology only`, y = compEffSize$`Pathology + Age + Genotype`,
     xlab= "Pathology only", ylab = "Pathology + Age + Genotype", 
     main = "rTg4510 Hip Effect sizes")
abline(v = 0)
abline(h = 0)


#############################   J20   ECX   #########################
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_ECX_J20.RData")

##load data
# Filter for J20 samples
pheno <- QCmetrics[which(QCmetrics$Tissue == "CortexEntorhinalis" & QCmetrics$AD_model == "J20"),]
betas <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno$Basename ]
pheno$Genotype <- as.factor(pheno$Genotype)
# Load pathology data to set ages
pathology_pheno <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_ECX_Pathology.csv", header = T, stringsAsFactors = F)
colnames(pathology_pheno)[colnames(pathology_pheno) == 'X...Sample_ID_ECX'] <- 'Sample_ID_ECX'
pheno$Sample_ID_ECX <- gsub(".*_", "", pheno$ExternalSampleID)
phenowPath <- merge(pheno,pathology_pheno, by = "Sample_ID_ECX")
colnames(phenowPath)[colnames(phenowPath) == 'Genotype.x'] <- 'Genotype'
pheno <- phenowPath
betas <- betas[,order(colnames(betas))]
pheno <- pheno[order(pheno$Basename),]
identical(colnames(betas), pheno$Basename)


### model 1 : Pathology
ResultsPathologyCleaned <- cleantable(betaResults = ResultsPathology, species_probes,mm10_Manifest)

### model 2 : Pathology+Age
ResultsPathologyAgeCleaned <- cleantable(betaResults = ResultsPathologyAge, species_probes,mm10_Manifest)

### model 3 : Pathology+Age+Genotyp
ResultsPathologyAgeGenotypeCleaned <- cleantable(betaResults = ResultsPathologyAgeGenotype, species_probes,mm10_Manifest)

## Tabulate number of significant hits
PathTbl <- c(nrow(ResultsPathologyCleaned[which(ResultsPathologyCleaned$FDR_adj_pathology <= 0.05),]),
             nrow(ResultsPathologyAgeCleaned[which(ResultsPathologyAgeCleaned$FDR_adj_pathology <= 0.05),]),
             nrow(ResultsPathologyAgeGenotypeCleaned[which(ResultsPathologyAgeGenotypeCleaned$FDR_adj_pathology <= 0.05),]))
PathTbl <- as.data.frame(t(PathTbl))
colnames(PathTbl) <- c("Pathology Only", "Pathology + Age", "Pathology + Age + Genotype")
p<-tableGrob(PathTbl)
grid.arrange(top="J20 ECX",p)


pathplots(betaResults = ResultsPathologyCleaned, Pathology_col = pheno$Pathology_ECX, col = "#FF5A62")
pathplots(betaResults = ResultsPathologyAgeCleaned, Pathology_col = pheno$Pathology_ECX, col = "#FF5A62")
pathplots(betaResults = ResultsPathologyAgeGenotypeCleaned, Pathology_col = pheno$Pathology_ECX, col = "#FF5A62")

## how many in path+age+geno in pathology only
ResultsPathologyAgeGenotypeCleanedSig <- ResultsPathologyAgeGenotypeCleaned[which(ResultsPathologyAgeGenotypeCleaned$FDR_adj_pathology <= 0.05),]
ResultsPathologyCleanedSig <- ResultsPathologyCleaned[which(ResultsPathologyCleaned$FDR_adj_pathology <= 0.05),]
length(intersect(ResultsPathologyAgeGenotypeCleanedSig$cpg, ResultsPathologyCleanedSig$cpg)) #3476
sharedsites <- intersect(ResultsPathologyAgeGenotypeCleanedSig$cpg, ResultsPathologyCleanedSig$cpg) 

sharedsitespathfull<- ResultsPathologyAgeGenotypeCleanedSig[which(ResultsPathologyAgeGenotypeCleanedSig$cpg %in% sharedsites),] 
sharedsitespathonly<- ResultsPathologyCleanedSig[which(ResultsPathologyCleanedSig$cpg %in% sharedsites),] 
identical(sharedsitespathonly$cpg, sharedsitespathfull$cpg)

par(mfrow = c(1,2))
compPvals <- cbind(sharedsitespathonly$FDR_adj_pathology, sharedsitespathfull$FDR_adj_pathology)
compPvals <- as.data.frame(compPvals)
colnames(compPvals) <- c("Pathology only", "Pathology + Age + Genotype")
plot(x = -log(compPvals$`Pathology only`), y = -log(compPvals$`Pathology + Age + Genotype`),
     xlab= "Pathology only", ylab = "Pathology + Age + Genotype", 
     main = "J20 ECX -log10 FDR Pvalues")

compEffSize <- cbind(sharedsitespathonly$estimate.Pathology, sharedsitespathfull$estimate.Pathology)
compEffSize <- as.data.frame(compEffSize)
colnames(compEffSize) <- c("Pathology only", "Pathology + Age + Genotype")
plot(x = compEffSize$`Pathology only`, y = compEffSize$`Pathology + Age + Genotype`,
     xlab= "Pathology only", ylab = "Pathology + Age + Genotype", 
     main = "J20 ECX Effect sizes")
abline(v = 0)
abline(h = 0)

#############################   J20   HIP  #########################
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/PathologyResults_HIP_J20.RData")

#load data
pheno <- QCmetrics[which(QCmetrics$Tissue == "Hippocampus" & QCmetrics$AD_model == "J20"),]
betas <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno$Basename ]
pheno$Genotype <- as.factor(pheno$Genotype)
# Load pathology data to set ages
pathology_pheno <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/J20_coldata_VertebrateArray.csv", header = T, stringsAsFactors = F)
pheno$Sample_ID_HIP <- gsub(".*_", "", pheno$ExternalSampleID) 
phenowPath <- merge(pheno,pathology_pheno, by = "Sample_ID_HIP")
colnames(phenowPath)[colnames(phenowPath) == 'Genotype.x'] <- 'Genotype'
pheno <- phenowPath
betas <- betas[,order(colnames(betas))]
pheno <- pheno[order(pheno$Basename),]
identical(colnames(betas), pheno$Basename)



### model 1 : Pathology
ResultsPathologyCleaned <- cleantable(betaResults = ResultsPathology, species_probes,mm10_Manifest)

### model 2 : Pathology+Age
ResultsPathologyAgeCleaned <- cleantable(betaResults = ResultsPathologyAge, species_probes,mm10_Manifest)

### model 3 : Pathology+Age+Genotyp
ResultsPathologyAgeGenotypeCleaned <- cleantable(betaResults = ResultsPathologyAgeGenotype, species_probes,mm10_Manifest)

## Tabulate number of significant hits
PathTbl <- c(nrow(ResultsPathologyCleaned[which(ResultsPathologyCleaned$FDR_adj_pathology <= 0.05),]),
             nrow(ResultsPathologyAgeCleaned[which(ResultsPathologyAgeCleaned$FDR_adj_pathology <= 0.05),]),
             nrow(ResultsPathologyAgeGenotypeCleaned[which(ResultsPathologyAgeGenotypeCleaned$FDR_adj_pathology <= 0.05),]))
PathTbl <- as.data.frame(t(PathTbl))
colnames(PathTbl) <- c("Pathology Only", "Pathology + Age", "Pathology + Age + Genotype")
p<-tableGrob(PathTbl)
grid.arrange(top="J20 HIP",p)


pathplots(betaResults = ResultsPathologyCleaned, Pathology_col = pheno$Pathology_HIP, col = "#FF5A62")
pathplots(betaResults = ResultsPathologyAgeCleaned, Pathology_col = pheno$Pathology_HIP, col = "#FF5A62")
pathplots(betaResults = ResultsPathologyAgeGenotypeCleaned, Pathology_col = pheno$Pathology_HIP, col = "#FF5A62")

## how many in path+age+geno in pathology only
ResultsPathologyAgeGenotypeCleanedSig <- ResultsPathologyAgeGenotypeCleaned[which(ResultsPathologyAgeGenotypeCleaned$FDR_adj_pathology <= 0.05),]
ResultsPathologyCleanedSig <- ResultsPathologyCleaned[which(ResultsPathologyCleaned$FDR_adj_pathology <= 0.05),]
length(intersect(ResultsPathologyAgeGenotypeCleanedSig$cpg, ResultsPathologyCleanedSig$cpg)) #10
sharedsites <- intersect(ResultsPathologyAgeGenotypeCleanedSig$cpg, ResultsPathologyCleanedSig$cpg) 

sharedsitespathfull<- ResultsPathologyAgeGenotypeCleanedSig[which(ResultsPathologyAgeGenotypeCleanedSig$cpg %in% sharedsites),] 
sharedsitespathonly<- ResultsPathologyCleanedSig[which(ResultsPathologyCleanedSig$cpg %in% sharedsites),] 
identical(sharedsitespathonly$cpg, sharedsitespathfull$cpg)

par(mfrow = c(1,2))
compPvals <- cbind(sharedsitespathonly$FDR_adj_pathology, sharedsitespathfull$FDR_adj_pathology)
compPvals <- as.data.frame(compPvals)
colnames(compPvals) <- c("Pathology only", "Pathology + Age + Genotype")
plot(x = -log(compPvals$`Pathology only`), y = -log(compPvals$`Pathology + Age + Genotype`),
     xlab= "Pathology only", ylab = "Pathology + Age + Genotype", 
     main = "J20 HIP -log10 FDR Pvalues")

compEffSize <- cbind(sharedsitespathonly$estimate.Pathology, sharedsitespathfull$estimate.Pathology)
compEffSize <- as.data.frame(compEffSize)
colnames(compEffSize) <- c("Pathology only", "Pathology + Age + Genotype")
plot(x = compEffSize$`Pathology only`, y = compEffSize$`Pathology + Age + Genotype`,
     xlab= "Pathology only", ylab = "Pathology + Age + Genotype", 
     main = "J20 HIP Effect sizes")
abline(v = 0)
abline(h = 0)

dev.off()


