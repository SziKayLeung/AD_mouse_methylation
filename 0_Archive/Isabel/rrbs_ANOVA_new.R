#Isabel Castanho I.Castanho@exeter.ac.uk

# RRBS ANOVA using filtered matrix from Dorothea ################## Nothing is significant

# /usr/lib64/R3.5.1/bin/R

# setwd
setwd("/mnt/data1/isabel/RRBS/")


##### Tg4510 #####

# column data / phenotypes
coldata <- read.csv("/mnt/data1/isabel/RRBS/Tg4510_phenotype_RRBS.csv", row.names=1, stringsAsFactors=FALSE)
colnames(coldata)
rownames(coldata)
#removed T24 ######### ASK THEA WHY IT'S NOT THERE ###############################################################
coldata <- coldata[1:61,]
coldata$Age_months <- as.factor(coldata$Age_months)
coldata$Genotype <- as.factor(coldata$Genotype)
coldata$Genotype <- relevel(coldata$Genotype, "WT")
levels(coldata$Genotype)

# get methylation data
load(file="/mnt/data1/Thea/IsabelSamples/data/rTg4510MatrixFiltered.Rdata")
#ls()
methylation_data <- as.data.frame(RRBS)
rownames(methylation_data) <- RRBS[,1]
methylation_data <- methylation_data[,2:62]
colnames(methylation_data) <- rownames(coldata)

### STATISTICS ###
# number of sites -> important for downstream code
n_rows <- nrow(methylation_data)

## FULL MODEL
stats_results <- matrix(data=NA, ncol=12, nrow=n_rows)

colnames(stats_results) <- c("t-value_Genotype", "Effect-size_Genotype", "P-value_Genotype", "FDR_Genotype", 
	"F-value_Age", "Effect-size_Age", "P-value_Age", "FDR_Age", 
	"F-value_Interaction", "Effect-size_Interaction", "P-value_Interaction", "FDR_Interaction")

for (i in 1:nrow(methylation_data)) {
	H1 <- lm(as.numeric(methylation_data[i,]) ~ coldata$Genotype + coldata$Age_months + coldata$Genotype *coldata$Age_months)
	H0 <- lm(as.numeric(methylation_data[i,]) ~ coldata$Genotype + coldata$Age_months)
	H0_2 <- lm(as.numeric(methylation_data[i,]) ~ coldata$Genotype)
	stats_results[i,1] <- summary(H0)$coefficients[2,3] # t-value genotype 
	stats_results[i,2] <- coef(H0)[2] # Effect size genotype 
	stats_results[i,3] <- summary(H0)$coefficients[2,4] # p-value genotype 
	# stats_results[i,4] will have FDR genotype
	stats_results[i,5] <- anova(H0,H0_2)[2,5] # F-value age
	stats_results[i,6] <- coef(H0)[5] # Effect size age (time-point 4 compared to baseline)
	stats_results[i,7] <- anova(H0,H0_2)[2,6] # p-value age
	# stats_results[i,8] will have FDR age
	stats_results[i,9] <- anova(H0,H1)[2,5] # F-value interaction 
	stats_results[i,10] <- coef(H1)[8] # Effect size interaction (time-point 4 compared to baseline)
	stats_results[i,11] <- anova(H0,H1)[2,6] # p-value interaction
	# stats_results[i,12] will have FDR interaction
}

rownames(stats_results) <- rownames(methylation_data)

FDR_adj_genotype <- p.adjust(stats_results[,"P-value_Genotype"], method = "fdr")
FDR_adj_age <- p.adjust(stats_results[,"P-value_Age"], method = "fdr")
FDR_adj_interaction <- p.adjust(stats_results[,"P-value_Interaction"], method = "fdr")

stats_table <-cbind(
	FDR_adj_genotype, as.data.frame(stats_results[,c("P-value_Genotype", "Effect-size_Genotype", "t-value_Genotype")]), 
	FDR_adj_age, as.data.frame(stats_results[,c("P-value_Age", "Effect-size_Age", "F-value_Age")]),	
	FDR_adj_interaction, as.data.frame(stats_results[,c("P-value_Interaction", "Effect-size_Interaction", "F-value_Interaction")])
	)
rownames(stats_table) <- rownames(methylation_data)
write.csv(stats_table, file = "RRBS_ANOVAstats_table.csv")

threshold <- 0.05
sig_genotype <- stats_table[which(stats_table[,"FDR_adj_genotype"]<threshold),]
#write.csv(sig_genotype, file = "Tg4510_RRBS_sig_genotype.csv")
sig_age <- stats_table[which(stats_table[,"FDR_adj_age"]<threshold),]
#write.csv(sig_age, file = "Tg4510_RRBS sig_age.csv")
sig_interaction <- stats_table[which(stats_table[,"FDR_adj_interaction"]<threshold),]
#write.csv(sig_interaction, file = "Tg4510_RRBS_sig_interaction.csv")


## GENOTYPE
stats_results <- matrix(data=NA, ncol=4, nrow=n_rows)

colnames(stats_results) <- c("t-value_Genotype", "Effect-size_Genotype", "P-value_Genotype", "FDR_Genotype")

for (i in 1:nrow(methylation_data)) {
	H1 <- lm(as.numeric(methylation_data[i,]) ~ coldata$Genotype) # 58 degrees of freedom (1 and 58 DF)
	stats_results[i,1] <- summary(H1)$coefficients[2,3] # t-value genotype 
	stats_results[i,2] <- coef(H1)[2] # Effect size genotype 
	stats_results[i,3] <- summary(H1)$coefficients[2,4] # p-value genotype 
}

rownames(stats_results) <- rownames(methylation_data)

FDR_adj_genotype <- p.adjust(stats_results[,"P-value_Genotype"], method = "fdr")

#stats_table <-cbind(FDR_adj_genotype, as.data.frame(stats_results[,c("P-value_Genotype", "Effect-size_Genotype", "t-value_Genotype")]))
stats_table <-cbind(stats_results, FDR_adj_genotype)

rownames(stats_table) <- rownames(methylation_data)
write.csv(stats_table, file = "RRBS_ANOVA_GenotypeStatsTable.csv")

threshold <- 0.05
sig_genotype <- stats_table[which(stats_table[,"FDR_adj_genotype"]<threshold),]
write.csv(sig_genotype, file = "Tg4510_RRBS_sig_genotype.csv")