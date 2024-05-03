## compare FDR and bonferroni to see if there is the same number of significant hits
## Normally array uses BonF but can we use FDR is the sample of signigicants are similiar ??
setwd("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/")
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_rTg4510.RData")


betaResults <- betaResultsChip
# Add chr and bp before deleting rows
# Add annotation files - chr and bp
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_mappability_name = "MusMusculus"
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #23633
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- gsub(".*:","",species_probes$MusMusculus)
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


########### Number of FDR #####

fdr_sig_genotype <- betaResults[which(betaResults[,"FDR_adj_genotype"] < 0.05),]
fdr_sig_age <- betaResults[which(betaResults[,"FDR_adj_age"] < 0.05),]
fdr_sig_interaction <- betaResults[which(betaResults[,"FDR_adj_Interaction"] < 0.05),]

########## Number of Bonf #####
bonfP <- 0.05/nrow(betaResultsChip) #use species probes as that includes non autosomal
bonf_sig_genotype <- betaResults[which(betaResults[,"p.val.Genotype"] < bonfP),]
bonf_sig_age <- betaResults[which(betaResults[,"p.val.Age"] < bonfP),]
bonf_sig_interaction <- betaResults[which(betaResults[,"p.val.Interaction"] < bonfP),]


res <- matrix(NA, nrow = 3, ncol = 2)
colnames(res) <- c("FDR", "Bonferroni")
rownames(res) <- c("Genotype", "Age", "Interaction")

res[1,1] = nrow(fdr_sig_genotype)
res[2,1] = nrow(fdr_sig_age)
res[3,1] = nrow(fdr_sig_interaction)
res[1,2] = nrow(bonf_sig_genotype)
res[2,2] = nrow(bonf_sig_age)
res[3,2] = nrow(bonf_sig_interaction)


#### Do the same for pathology

betaResults <- betaResultsPathology
# Add chr and bp before deleting rows
Chr <- species_probes$Chr[match(species_probes$probeID, betaResults[,"cpg"])]
Bp <- species_probes$Bp[match(species_probes$probeID, betaResults[,"cpg"])]
betaResults <- cbind(Bp, betaResults)
betaResults <- cbind(Chr, betaResults)
# Remove NA rows
# na_df <- as.data.frame(betaResults[rowSums(is.na(betaResults)) > 0,]) #54 probes has no data in columns
# betaResults <- betaResults[-which(betaResults[,"cpg"] %in% na_df$cpg),] #remove these NAs rows
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

########## Number of FDR & Bonf  #####
fdr_sig_pathology <- betaResults[which(betaResults[,"FDR_adj_pathology"] < 0.05),]
bonf_sig_pathology <- betaResults[which(betaResults[,"p.val.Pathology"] < bonfP),]

res_path <- c(nrow(fdr_sig_pathology), nrow(bonf_sig_pathology))
res <- rbind(res, res_path)
rownames(res) <-c("Genotype", "Age", "Interaction", "Pathology")
print(res)

library(gridExtra)
png("FDRvsBonf.png")
p<-tableGrob(res)
grid.arrange(p)
dev.off()


