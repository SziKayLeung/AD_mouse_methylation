## MME Results 

#raw data
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat")
QCmetrics$SampleID <- gsub(".*_", "", QCmetrics$ExternalSampleID)
###################################### J20  #######################################

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
pheno_J20_ECX <- cbind(pheno_J20_ECX, J20_path)
### Hip
pheno_J20_HIP <- QCmetrics[which(QCmetrics$Tissue == "Hippocampus" & QCmetrics$AD_model == "J20"),]
betas_J20_HIP <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno_J20_HIP$Basename ]
identical(pheno_J20_HIP$Basename, colnames(betas_J20_HIP))
colnames(betas_J20_HIP) <- pheno_J20_HIP$SampleID
rownames(pheno_J20_HIP) <- pheno_J20_HIP$SampleID 
pheno_J20_HIP <- pheno_J20_HIP[J20_path$Sample_ID_HIP, ] #order based on pathology samples id (NA rows will be generated)
pheno_J20_HIP <- cbind(pheno_J20_HIP, J20_path)

# Match the same mouse with 2 brain regions
samples <- as.data.frame(cbind(pheno_J20_ECX$SampleID, pheno_J20_HIP$SampleID))
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
#pheno_J20_ECX <- pheno_J20_ECX[samples$ECX,]
identical(pheno_J20_ECX$SampleID, samples$ECX)
pheno_J20_ECX$MouseID <- samples$MouseID[match(samples$ECX, pheno_J20_ECX$SampleID)] #Now the orders are the same - add in the mouseIDs
#betas_J20_ECX <- betas_J20_ECX[,pheno_J20_ECX$SampleID]
identical(colnames(betas_J20_ECX) , pheno_J20_ECX$SampleID)
pheno_J20_ECX <- pheno_J20_ECX[colnames(betas_J20_ECX), ]
identical(colnames(betas_J20_ECX) , pheno_J20_ECX$SampleID)


### hip
#pheno_J20_HIP <- pheno_J20_HIP[samples$HIP,]
identical(pheno_J20_HIP$SampleID, samples$HIP)
pheno_J20_HIP$MouseID <- samples$MouseID[match(samples$HIP, pheno_J20_HIP$SampleID)] #Now the orders are the same - add in the mouseIDs
#betas_J20_HIP <- betas_J20_HIP[,pheno_J20_HIP$SampleID]
identical(colnames(betas_J20_HIP) , pheno_J20_HIP$SampleID)
pheno_J20_HIP <- pheno_J20_HIP[colnames(betas_J20_HIP), ]
identical(colnames(betas_J20_HIP) , pheno_J20_HIP$SampleID)


pheno <- rbind(pheno_J20_ECX, pheno_J20_HIP)
betas <- cbind(betas_J20_ECX, betas_J20_HIP)
identical(pheno$SampleID, colnames(betas))


pheno$Genotype <- factor(pheno$Genotype, levels = c("WT","TG"))
pheno$Tissue <- factor(pheno$Tissue, levels = c("CortexEntorhinalis","Hippocampus"))
pheno$Chip_ID <- as.factor(pheno$Chip_ID)
pheno$MouseID <- as.factor(pheno$MouseID)


betas <- as.matrix(betas)

pheno_J20 <- pheno
betas_J20 <- betas

pheno <- pheno[,c("Tissue","Genotype", "Age_months", "Chip_ID", "PlateNumber", "MouseID", "Basename", "Pathology_ECX", "Pathology_HIP")]
betas <- as.matrix(betas)
rownames(pheno) <- pheno$Basename
colnames(betas) <- rownames(pheno)

##########



mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)

J20 <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/MixedModelResults_J20.csv",
                header=T, stringsAsFactors = F)
colnames(J20)[colnames(J20) == 'X'] <- 'cpg'

# Add annotation files - chr and bp
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_mappability_name = "MusMusculus"
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #23633
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- gsub(".*:","",species_probes$MusMusculus)

J20$Chr <- species_probes$Chr[match(species_probes$probeID, J20[,"cpg"])]
J20$Bp <- species_probes$Bp[match(species_probes$probeID, J20[,"cpg"])]

J20<- J20[-which(J20$Chr %in% c("CHR_MG51_PATCH",
                                "CHR_MG4200_PATCH","CHR_MG3699_PATCH")),]
J20$Chr <- gsub("X", 20, J20$Chr)
J20$Chr <- gsub("Y", 21, J20$Chr)
J20$Chr <- as.numeric(J20$Chr)
J20$Bp <- as.numeric(J20$Bp)
J20$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(J20$cpg, mm10_Manifest$cpg)]
rownames(J20) <- J20$cpg
J20$SNP <- J20$Gene_Symbol


#qqplot and manhattan
pdf("MME_J20_qqplot.pdf")
par(mfrow = c(2,2))
lamda <- qchisq(1-median(J20$PrZ.Tissue),1)/qchisq(0.5,1)
qq(J20$PrZ.Tissue, main = "J20 Tissue")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(J20$PrZ.GenotypeTG),1)/qchisq(0.5,1)
qq(J20$PrZ.GenotypeTG, main = "J20 Genotype")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(J20$PrZ.TissueHippocampus.GenotypeTG),1)/qchisq(0.5,1)
qq(J20$PrZ.TissueHippocampus.GenotypeTG, main = "J20 Tissue*Genotype")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(J20$PrZ.GenotypeTG.Age_months),1)/qchisq(0.5,1)
qq(J20$PrZ.GenotypeTG.Age_months, main = "J20 Genotype*Age")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)
dev.off()

manhattan(J20, p = "FDR_adj_tissue", bp = "Bp", chr = "Chr", 
          genomewide = -log10(threshold), suggestiveline = FALSE, 
          logp=T, col = c("black", color_J20_TG), main = "Tissue",
          chrlabs = c(1:19, "X", "Y"),las = 2, annotatePval = threshold,
          cex.axis = 1.5, cex.lab = 1.2)


SigJ20 <- J20[which(J20$FDR_adj_tissue < 0.05),] #17471 significant hits

## pull out numbers etc
J20$FDR_adj_GenotypeTG <- p.adjust(J20[,"PrZ.GenotypeTG"], method = "fdr")
J20$FDR_adj_GenotypeAge <- p.adjust(J20[,"PrZ.GenotypeTG.Age_months"], method = "fdr")
J20$FDR_adj_TissueGenotype <- p.adjust(J20[,"PrZ.TissueHippocampus.GenotypeTG"], method = "fdr")


nrow(J20[which(J20$FDR_adj_GenotypeTG < 0.05),])     #75
nrow(J20[which(J20$FDR_adj_GenotypeAge < 0.05),])    #1971
nrow(J20[which(J20$FDR_adj_TissueGenotype < 0.05),]) #5229

######################## Are the effects direction and magnitude the same as within seperate tissues?
#J20
#cortex 
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_J20.RData")
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
J20_cortex <- betaResults


#J20
#hippocampus
load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/ArrayBetaRegResults_J20_HIP.RData")
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
J20_hipp <- betaResults


############################### Comparing effect sizes - no filtering for sig sites
# Check that they are all in the same order of rows so they are comparing the same things
J20 <- J20[order(J20$cpg),]
J20_cortex <- J20_cortex[order(J20_cortex$cpg),]
J20_hipp <- J20_hipp[order(J20_hipp$cpg),]

identical(J20$cpg, J20_cortex$cpg)
identical(J20$cpg, J20_hipp$cpg)

pdf("J20_AllEffectSizes_braincortexhippocampus.pdf")
par(mfrow = c(2,2))
# cortex vs brain
cor.test(J20$Betas.GenotypeTG, J20_cortex$estimate.Genotype, method = "pearson")
corr <- cor(J20$Betas.GenotypeTG, J20_cortex$estimate.Genotype, method = "pearson")
plot(J20$Betas.GenotypeTG, J20_cortex$estimate.Genotype, 
     xlab = "Cortex", ylab = "Brain",
     main = "Effect sizes: J20", cex = 0.7)
mtext(paste("cor = ", signif(corr,3)), side = 3, adj = 1)
abline(0,1, col="blue")
abline(h = 0)
abline(v = 0)

#hipp vs brain
cor.test(J20$Betas.GenotypeTG, J20_hipp$estimate.Genotype, method = "pearson")
corr <- cor(J20$Betas.GenotypeTG, J20_hipp$estimate.Genotype, method = "pearson")
plot(J20$Betas.GenotypeTG, J20_hipp$estimate.Genotype, 
     xlab = "Hippocampus", ylab = "Brain",
     main = "Effect sizes: J20", cex = 0.7)
mtext(paste("cor = ", signif(corr,3)), side = 3, adj = 1)
abline(0,1, col="blue")
abline(h = 0)
abline(v = 0)

# cortex vs hipp
cor.test(J20_cortex$estimate.Genotype, J20_hipp$estimate.Genotype, method = "pearson")
corr <- cor(J20_cortex$estimate.Genotype, J20_hipp$estimate.Genotype, method = "pearson")
plot(J20_cortex$estimate.Genotype, J20_hipp$estimate.Genotype, 
     xlab = "Cortex", ylab = "Hippocampus",
     main = "Effect sizes: J20", cex = 0.7)
mtext(paste("cor = ", signif(corr,3)), side = 3, adj = 1)
abline(0,1, col="blue")
abline(h = 0)
abline(v = 0)
dev.off()

############################### Comparing effect sizes -  filtering for sig sites
# Add FDR P value to brain
J20$FDR_adj_GenotypeTG <- p.adjust(J20[,"PrZ.GenotypeTG"], method = "fdr")
J20_sig <- J20[which(J20$FDR_adj_GenotypeTG <= 0.05), "cpg"] #2 sites for brain
J20_cortex_sig <- J20_cortex[which(J20_cortex$FDR_adj_genotype <= 0.05), "cpg"] #12 sites for cortex
J20_hipp_sig <- J20_hipp[which(J20_hipp$FDR_adj_genotype <= 0.05), "cpg"] #72 sites forr hip

length(intersect(J20_sig, J20_cortex_sig)) #1 for brain & cortex
length(intersect(J20_sig, J20_hipp_sig)) #0 for brain & hippocampus
length(intersect(J20_cortex_sig, J20_hipp_sig)) #0 for cortex & hippocampus

cortexbrain <- intersect(J20_sig, J20_cortex_sig)
hippbrain <- intersect(J20_sig, J20_hipp_sig)
cortexhipp <- intersect(J20_cortex_sig, J20_hipp_sig)

### Only one site for MME and J20 Cortex - Zbtb20


J20<-J20[order(J20$FDR_adj_GenotypeTG),]
J20_sig <- J20[1:100,]
betas_sig <- betas[rownames(J20_sig),]
pheno_plot <- pheno[,c("Tissue","Genotype", "Age_months")]

pdf("J20_Heatmap_sig_genotype.pdf")
pheatmap(betas_sig,  
         annotation_col = pheno_plot,
         cluster_rows = T,
         show_rownames = F,
         show_colnames = F)

J20_sig <- J20[1:2,]
betas_sig <- betas[rownames(J20_sig),]
pheatmap(betas_sig,  
         annotation_col = pheno_plot,
         cluster_rows = T,
         show_rownames = F,
         show_colnames = F)
dev.off()

####################### Interactions
# how many significant interactions
J20$FDR_adj_TissueGenoInt <- p.adjust(J20[,"PrZ.TissueHippocampus.GenotypeTG"], method = "fdr")
sig_tissueint <- J20[which(J20$FDR_adj_TissueGenoInt < 0.05),] 
nrow(sig_tissueint) #248

# How many are significant for interaction and genotype??
length(intersect(sig_tissueint$cpg, J20_cortex_sig)) #1 sites for cortex
length(intersect(sig_tissueint$cpg, J20_hipp_sig)) #9 sites for hip


intcortex <- intersect(sig_tissueint$cpg, J20_cortex_sig) 
inthipp <- intersect(sig_tissueint$cpg, J20_hipp_sig)



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
J20<-J20[order(J20$PrZ.TissueHippocampus.GenotypeTG),]
J20_sig <- J20[1:nrow(sig_tissueint),]
betas_sig <- betas[rownames(J20_sig),]
pheno_plot <- pheno[,c("Tissue","Genotype", "Age_months")]

pdf("rTg_Heatmap_sig_tissuegenotype.pdf")
pheatmap(betas_sig,  
         annotation_col = pheno_plot,
         cluster_rows = T,
         show_rownames = F,
         show_colnames = F)
dev.off()



#### Genotype*Age
J20$FDR_adj_GenotypeAgeInt <- p.adjust(J20[,"PrZ.GenotypeTG.Age_months"], method = "fdr")
sig_genoageint <- J20[which(J20$FDR_adj_GenotypeAgeInt < 0.05),] 
nrow(sig_genoageint) #0


############################################# Pathology data ##################################################################

J20 <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/MixedModelResultswPathology_J20.csv",
                    header=T, stringsAsFactors = F)

mm10_Manifest <- read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)

colnames(J20)[colnames(J20) == 'X'] <- 'cpg'



# Add annotation files - chr and bp
probe_mapping_file <- read_csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.0.csv", col_types = cols(.default = "c"))
species_mappability_name = "MusMusculus"
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')]
species_probes <- as.data.frame(species_probes[!is.na(species_probes$MusMusculus),]) #23633
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$Bp <- gsub(".*:","",species_probes$MusMusculus)

J20$Chr <- species_probes$Chr[match(species_probes$probeID, J20[,"cpg"])]
J20$Bp <- species_probes$Bp[match(species_probes$probeID, J20[,"cpg"])]


# Convert x and y to numbers and remove non autosomal
J20<- J20[-which(J20$Chr %in% c("CHR_MG51_PATCH", "CHR_MG4200_PATCH","CHR_MG3699_PATCH")),]
J20$Chr <- gsub("X", 20, J20$Chr)
J20$Chr <- gsub("Y", 21, J20$Chr)
J20$Chr <- as.numeric(J20$Chr)
J20$Bp <- as.numeric(J20$Bp)
J20$Gene_Symbol <- mm10_Manifest$Gene_Symbol[match(J20$cpg, mm10_Manifest$cpg)]
rownames(J20) <- J20$cpg
J20$SNP <- J20$Gene_Symbol


par(mfrow = c(2,2))
lamda <- qchisq(1-median(J20$PrZ.GenotypeTG),1)/qchisq(0.5,1)
qq(J20$PrZ.GenotypeTG, main = "J20 Genotype")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(J20$PrZ.TissueHippocampus.GenotypeTG),1)/qchisq(0.5,1)
qq(J20$PrZ.TissueHippocampus.GenotypeTG, main = "J20 Tissue*Genotype")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(J20$PrZ.Pathology_ECX),1)/qchisq(0.5,1)
qq(J20$PrZ.Pathology_ECX, main = "J20 Pathology_ECX")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)

lamda <- qchisq(1-median(J20$PrZ.Pathology_HIP),1)/qchisq(0.5,1)
qq(J20$PrZ.Pathology_HIP, main = "J20 Pathology_HIP")
mtext(paste("Lamda = ", signif(lamda,3)), side = 3, adj = 1)


## pull out numbers etc
J20$FDR_adj_GenotypeTG <- p.adjust(J20[,"PrZ.GenotypeTG"], method = "fdr")
J20$FDR_adj_TissueGenotype <- p.adjust(J20[,"PrZ.TissueHippocampus.GenotypeTG"], method = "fdr")
J20$FDR_adj_Pathology_ECX <- p.adjust(J20[,"PrZ.Pathology_ECX"], method = "fdr")
J20$FDR_adj_Pathology_HIP <- p.adjust(J20[,"PrZ.Pathology_HIP"], method = "fdr")


nrow(J20[which(J20$FDR_adj_GenotypeTG < 0.05),])     #26
nrow(J20[which(J20$FDR_adj_TissueGenotype < 0.05),]) #157
nrow(J20[which(J20$FDR_adj_Pathology_ECX < 0.05),])  #0
nrow(J20[which(J20$FDR_adj_Pathology_HIP < 0.05),])  #0

##### significant pathology sites across tissues


##### significant genotype
J20_genosig <- J20[which(J20$FDR_adj_GenotypeTG < 0.05),]
J20_genosig <- J20_genosig[order(J20_genosig$Betas.GenotypeTG),]
J20_genosig_table <- J20_genosig[,c("cpg", "Gene_Symbol", "Betas.GenotypeTG","SE.GenotypeTG","PrZ.GenotypeTG","FDR_adj_GenotypeTG")]
colnames(J20_genosig_table) <- c("cpg", "Gene_Symbol", "Beta","SE", "Raw.Pvalues", "FDR.Adj.Pvalues")
write.csv(J20_genosig_table, "J20_MME_Genotype.csv")


pdf("J20_Genotype_MME.pdf", onefile = T)
par(mfrow = c(2,2))
for( i in 1:6){
  par(mfrow = c(2,2))
  betaResults <- J20_genosig
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
    scale_color_manual(values=c("black", "#FF5A62")) +
    labs(y = "Methylation", title = paste(betaResults$Gene_Symbol[i], betaResults$cpg[i], sep = " - ")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(p)
}
dev.off()


##### significant genotype*tissue
J20_genotisssig <- J20[which(J20$FDR_adj_TissueGenotype < 0.05),]
J20_genotisssig <- J20_genotisssig[order(J20_genotisssig$Betas.TissueHippocampus.GenotypeTG),]
J20_genotisssig_table <- J20_genotisssig[,c( "cpg", "Gene_Symbol", "Betas.TissueHippocampus.GenotypeTG",
                                                     "SE.TissueHippocampus.GenotypeTG","PrZ.TissueHippocampus.GenotypeTG",
                                                     "FDR_adj_TissueGenotype")]
colnames(J20_genotisssig_table) <- c("cpg", "Gene_Symbol", "Beta","SE", "Raw.Pvalues", "FDR.Adj.Pvalues")
write.csv(J20_genotisssig_table, "J20_MME_TissueGenotypeInteraction.csv")

pdf("J20_GenotypeTissue_MME.pdf", onefile = T)
par(mfrow = c(2,2))
for( i in 1:6){
  par(mfrow = c(2,2))
  betaResults <- J20_genotisssig
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
    scale_color_manual(values=c("black", "#FF5A62")) +
    labs(y = "Methylation", title = paste(betaResults$Gene_Symbol[i], betaResults$cpg[i], sep = " - ")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(p)
}
dev.off()



# how many significant pathology - ECX

J20_pathecxsig <- J20[which(J20$FDR_adj_Pathology_ECX < 0.05),]
J20_pathecxsig <- J20_pathecxsig[order(J20_pathecxsig$FDR_adj_Pathology_ECX),]
J20_pathecxsig_table <- J20_pathecxsig[,c( "cpg", "Gene_Symbol", "Betas.Pathology_ECX","SE.Pathology_ECX","PrZ.Pathology_ECX", "FDR_adj_Pathology_ECX")]
colnames(J20_pathecxsig_table) <- c("cpg", "Gene_Symbol", "Beta","SE", "Raw.Pvalues", "FDR.Adj.Pvalues")
write.csv(J20_pathecxsig_table, "J20_MME_Pathology_ECX.csv")

pdf("J20_PathECX_MME.pdf", onefile = T)
par(mfrow = c(2,2))
for( i in 1:6){
  par(mfrow = c(2,2))
  betaResults <- J20_pathecxsig
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
    scale_color_manual(values=c("black", "#FF5A62")) +
    labs(y = "Methylation", title = paste(betaResults$Gene_Symbol[i], betaResults$cpg[i], sep = " - ")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(p)
}
dev.off()


# how many significant pathology - HIP
J20_pathhipsig <- J20[which(J20$FDR_adj_Pathology_HIP < 0.05),]
J20_pathhipsig <- J20_pathhipsig[order(J20_pathhipsig$FDR_adj_Pathology_HIP),]
J20_pathhipsig_table <- J20_pathhipsig[,c( "cpg", "Gene_Symbol", "Betas.Pathology_HIP","SE.Pathology_HIP","PrZ.Pathology_HIP", "FDR_adj_Pathology_HIP")]
colnames(J20_pathhipsig_table) <- c("cpg", "Gene_Symbol", "Beta","SE", "Raw.Pvalues", "FDR.Adj.Pvalues")
write.csv(J20_pathhipsig_table, "J20_MME_Pathology_HIP.csv")

pdf("J20_PathHIP_MME.pdf", onefile = T)
par(mfrow = c(2,2))
for( i in 1:6){
  par(mfrow = c(2,2))
  betaResults <- J20_pathhipsig
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
    scale_color_manual(values=c("black", "#FF5A62")) +
    labs(y = "Methylation", title = paste(betaResults$Gene_Symbol[i], betaResults$cpg[i], sep = " - ")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(p)
}
dev.off()


