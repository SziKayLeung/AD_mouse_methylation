## 02.07.2020
## There are 1201 sites shared between RRBS data and Array data. Can we compare the stats are similiar at these sites
setwd("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/CompareTech/J20/")
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(cgwtools)

############################################  rTg #####################################

#################### Genotype ################

print(" -------------  rTg Genotype & Interaction ----------------")
rm(list=ls())

## load array
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_J20.RData")
rm(betaResultsPathology,betaResults, betaResults2,RRBSBetaResults_Shared, ArrayBetaResults_Shared)
ls()
# keep : betaResultsChip

## load rrbs
load("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/biseqResultsGenotype_J20.RData")
# keep : "DMRsGenotype"        "betaResultsGenotype"

ls()
#"DMRsGenotype"        "betaResultsChip"     "betaResultsGenotype"

# Clean betaResults chip files into clean results tables
#Load coordinates table
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #24048 probes 
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- species_probes$MusMusculus
species_probes$Bp <- gsub("*.:","", species_probes$Bp)
species_probes <- species_probes[,-2]
species_probes$Position <- paste("chr", paste(species_probes$Chr, species_probes$Bp, sep = ":"), sep ="")

#Clean Array table to genotype only
betaResults <- betaResultsChip
# Add chr and bp 
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
mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)
betaResults$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(betaResults$cpg, mm10_Manifest$cpg)]
colnames(betaResults)

#  [1] "cpg"                   "Chr"                   "Bp"
#  [4] "p.val.Genotype"        "meth.group1.WT"        "meth.group2.TG"
#  [7] "meth.diff.Genotype"    "pseudo.R.sqrt"         "estimate.Genotype"
# [10] "std.error.Genotype"    "FDR_adj_genotype"      "p.val.Age"
# [13] "estimate.Age"          "std.error.Age"         "FDR_adj_age"
# [16] "p.val.Interaction"     "estimate.Interaction"  "std.error.Interaction"
# [19] "p.val.modelLRT"        "FDR_adj_Interaction"   "Gene_Symbol"

ArrayGenotypeResults <- betaResults[,c(1:11,21)]
RRBSGenotypeResults <- betaResultsGenotype


# Load Interaction data
load("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/biseqResultsInteraction_J20.RData")
ArrayInteractionResults <- betaResults[,c(1:3,16:21)]
RRBSInteractionResults <- betaResultsInteraction
rm(list=setdiff(ls(), c("ArrayGenotypeResults","RRBSGenotypeResults",
						"ArrayInteractionResults","RRBSInteractionResults")))


# keep matching chr and pos rows  in both dfs
#clean rrbs chr col
RRBSGenotypeResults$chr <- as.character(RRBSGenotypeResults$chr)
RRBSGenotypeResults$chr<- gsub(".*chr", "", RRBSGenotypeResults$chr)
RRBSGenotypeResults$chr <- gsub("X", 20, RRBSGenotypeResults$chr)
RRBSGenotypeResults$chr <- gsub("Y", 21, RRBSGenotypeResults$chr)
RRBSGenotypeResults <- RRBSGenotypeResults[which(RRBSGenotypeResults$chr %in% c(as.character(1:21))),]
RRBSGenotypeResults$chr <- as.numeric(as.character(RRBSGenotypeResults$chr))

RRBSInteractionResults$chr <- as.character(RRBSInteractionResults$chr)
RRBSInteractionResults$chr<- gsub(".*chr", "", RRBSInteractionResults$chr)
RRBSInteractionResults$chr <- gsub("X", 20, RRBSInteractionResults$chr)
RRBSInteractionResults$chr <- gsub("Y", 21, RRBSInteractionResults$chr)
RRBSInteractionResults <- RRBSInteractionResults[which(RRBSInteractionResults$chr %in% c(as.character(1:21))),]
RRBSInteractionResults$chr <- as.numeric(as.character(RRBSInteractionResults$chr))


#change col names with A for Array and R for RRBS
colnames(ArrayGenotypeResults) <- paste("A", colnames(ArrayGenotypeResults), sep = ".")
colnames(RRBSGenotypeResults) <- paste("R", colnames(RRBSGenotypeResults), sep = ".")

colnames(ArrayInteractionResults) <- paste("A", colnames(ArrayInteractionResults), sep = ".")
colnames(RRBSInteractionResults) <- paste("R", colnames(RRBSInteractionResults), sep = ".")
# paste chr and bp to 'position' col

RRBSGenotypeResults$position <- paste(RRBSGenotypeResults$R.chr, RRBSGenotypeResults$R.pos, sep = ":")
ArrayGenotypeResults$position <- paste(ArrayGenotypeResults$A.Chr, ArrayGenotypeResults$A.Bp, sep = ":")
str(RRBSGenotypeResults)
str(ArrayGenotypeResults)

RRBSInteractionResults$position <- paste(RRBSInteractionResults$R.chr, RRBSInteractionResults$R.pos, sep = ":")
ArrayInteractionResults$position <- paste(ArrayInteractionResults$A.Chr, ArrayInteractionResults$A.Bp, sep = ":")
str(RRBSInteractionResults)
str(ArrayInteractionResults)


GenotypeResults <- merge(ArrayGenotypeResults, RRBSGenotypeResults, by = "position")
dim(GenotypeResults)
# 1201   23


InteractionResults <- merge(ArrayInteractionResults, RRBSInteractionResults, by = "position")
dim(InteractionResults)
# 1201   20


# Need to select effect sizes which are significant for genotype
# Add RRBS FDR GEnotype
R.FDR_adj_genotype <- p.adjust(GenotypeResults[,"R.p.val"], method = "fdr")
GenotypeResults$R.FDR_adj_genotype <- R.FDR_adj_genotype
rm(R.FDR_adj_genotype)

R.FDR_adj_Interaction <- p.adjust(InteractionResults[,"R.p.val"], method = "fdr")
InteractionResults$R.FDR_adj_Interaction <- R.FDR_adj_Interaction
rm(R.FDR_adj_Interaction)

# Select the FDR values
A.SignGenotypeResults <- GenotypeResults[which(GenotypeResults$A.FDR_adj_genotype <= 0.05),]
R.SignGenotypeResults <- GenotypeResults[which(GenotypeResults$R.FDR_adj_genotype <= 0.05),]
dim(A.SignGenotypeResults) #5 sites
dim(R.SignGenotypeResults) #0 sites

A.SignInteractionResults <- InteractionResults[which(InteractionResults$A.FDR_adj_Interaction <= 0.05),]
R.SignInteractionResults <- InteractionResults[which(InteractionResults$R.FDR_adj_Interaction <= 0.05),]
dim(A.SignInteractionResults) #26 sites
dim(R.SignInteractionResults) #0 sites
## Array tends to be more significant for both genotype and interaction.... what about pathology?


save(GenotypeResults, InteractionResults ,file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/CompareTech/J20/TechComparison_J20.RData")


##########################################   Pathology     ################


print(" -------------  rTg Pathology ----------------")
rm(list=ls())
## load array
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_J20.RData")
rm(betaResultsChip,betaResults, betaResults2,RRBSBetaResults_Shared, ArrayBetaResults_Shared)
ls()
betaResults <- betaResultsPathology
# keep : betaResultsPathology
## load rrbs
load("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/biseqResultsPathology_J20.RData")
ls()
rm(DMRsPathology)
# betaresults = array
# betaresultspathology = rrbs


# Clean betaResults chip files into clean results tables
#Load coordinates table
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #24048 probes 
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- species_probes$MusMusculus
species_probes$Bp <- gsub("*.:","", species_probes$Bp)
species_probes <- species_probes[,-2]
species_probes$Position <- paste("chr", paste(species_probes$Chr, species_probes$Bp, sep = ":"), sep ="")

#Clean Array table 
# Add chr and bp 
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]

betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
betaResults <- betaResults[,c(10,1:9)]

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
colnames(betaResults)

#  [1] "cpg"                   "Chr"                   "Bp"
#  [4] "p.val.Genotype"        "meth.group1.WT"        "meth.group2.TG"
#  [7] "meth.diff.Genotype"    "pseudo.R.sqrt"         "estimate.Genotype"
# [10] "std.error.Genotype"    "FDR_adj_genotype"      "p.val.Age"
# [13] "estimate.Age"          "std.error.Age"         "FDR_adj_age"
# [16] "p.val.Interaction"     "estimate.Interaction"  "std.error.Interaction"
# [19] "p.val.modelLRT"        "FDR_adj_Interaction"   "Gene_Symbol"

ArrayPathologyResults <- betaResults
RRBSPathologyResults <- betaResultsPathology

rm(list=setdiff(ls(), c("ArrayPathologyResults","RRBSPathologyResults")))


# keep matching chr and pos rows  in both dfs
#clean rrbs chr col
RRBSPathologyResults$chr <- as.character(RRBSPathologyResults$chr)
RRBSPathologyResults$chr<- gsub(".*chr", "", RRBSPathologyResults$chr)
RRBSPathologyResults$chr <- gsub("X", 20, RRBSPathologyResults$chr)
RRBSPathologyResults$chr <- gsub("Y", 21, RRBSPathologyResults$chr)
RRBSPathologyResults <- RRBSPathologyResults[which(RRBSPathologyResults$chr %in% c(as.character(1:21))),]
RRBSPathologyResults$chr <- as.numeric(as.character(RRBSPathologyResults$chr))

#change col names with A for Array and R for RRBS
colnames(ArrayPathologyResults) <- paste("A", colnames(ArrayPathologyResults), sep = ".")
colnames(RRBSPathologyResults) <- paste("R", colnames(RRBSPathologyResults), sep = ".")

# paste chr and bp to 'position' col
RRBSPathologyResults$position <- paste(RRBSPathologyResults$R.chr, RRBSPathologyResults$R.pos, sep = ":")
ArrayPathologyResults$position <- paste(ArrayPathologyResults$A.Chr, ArrayPathologyResults$A.Bp, sep = ":")
str(RRBSPathologyResults)
str(ArrayPathologyResults)

PathologyResults <- merge(ArrayPathologyResults, RRBSPathologyResults, by = "position")
dim(PathologyResults)
# 1201   23

R.FDR_adj_Pathology <- p.adjust(PathologyResults[,"R.p.val"], method = "fdr")
PathologyResults$R.FDR_adj_Pathology <- R.FDR_adj_Pathology
rm(R.FDR_adj_Pathology)

resave(PathologyResults ,file = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/CompareTech/J20/TechComparison_J20.RData")

