#Isabel Castanho I.Castanho@exeter.ac.uk

# RRBS heatmaps

# R 3.4.3

setwd("/mnt/data1/isabel/RRBS/")

library("DESeq2")
library("pheatmap")
library("WGCNA")
#library(vioplot)

color_Tg4510_TG <- "#00AEC9"
color_J20_TG <- "#FF5A62"
age1 <- "#FFB000"
age2 <- "#FE6100"
age3 <- "#DC267F"
age4 <- "#785EF0"

##### Tg4510 #####
# column data / phenotypes
coldata <- read.csv("/mnt/data1/isabel/RRBS/Tg4510_phenotype_RRBS.csv", row.names=1, stringsAsFactors=FALSE)
colnames(coldata)
rownames(coldata)
coldata$Age_months <- as.factor(coldata$Age_months)
coldata$Genotype <- as.factor(coldata$Genotype)
coldata$Genotype <- relevel(coldata$Genotype, "WT")
levels(coldata$Genotype)

# get methylation data
load(file="/mnt/data1/Thea/IsabelSamples/data/rTg4510GenotypeMatrixFiltered.Rdata") # object: RRBS
methylation <- RRBS[,2:63]
rownames(methylation) <- RRBS[,1]
colnames(methylation) <- sub("_m.*", "", colnames(methylation))
identical(colnames(methylation), rownames(coldata))
save(methylation, file = "/mnt/data1/isabel/RRBS/rTg4510methylation.RData")

# get gene list
load("/mnt/data1/isabel/RRBS/rTg4510GenotypeStatsManhattan.RData") # object: genotype_stats_clean

threshold <- 0.05
sig_genotype <- genotype_stats_clean[which(genotype_stats_clean[,"FDR"]<threshold),]
write.csv(sig_genotype, file = "Tg4510_RRBS_sig_genotype.csv")
#genotype_stats_FDR <- read.csv("/mnt/data1/isabel/RRBS/rTg4510GenotypeStats.csv", row.names=1, stringsAsFactors=FALSE) # not filtered for transgene-associated genes!

# code colours
ann_colors = list(
	Genotype=c(WT="black", TG=color_Tg4510_TG),
	Age_months=c("2"=age1, "4"=age2, "6"=age3, "8"=age4))

## Genotype
# clustered, methylation
sites_to_plot <- rownames(sig_genotype)
plotdata <- methylation[sites_to_plot,]# 20023 sites (rows)
plotdata <- plotdata[complete.cases(plotdata), ] # 3625 sites (rows)

pdf(file = "/mnt/data1/isabel/RRBS/Tg4510_genotype_heatmapclustered_methylation.pdf", width=10, height=15)
pheatmap(
	plotdata[, order(coldata$Genotype, coldata$Age_months)],
	color = blueWhiteRed(100)[0:100],
	cluster_rows=TRUE,
	cluster_cols=TRUE,
	annotation_col=coldata[ order(coldata$Genotype, coldata$Age_months),c("Genotype", "Age_months")],
	annotation_colors=ann_colors,
	show_colnames=FALSE,
	fontsize_row = 5,
	show_rownames = FALSE,
	scale="row")
dev.off()

##### J20 #####
# column data / phenotypes
coldata <- read.csv("/mnt/data1/isabel/RRBS/J20_phenotype_RRBS.csv", row.names=1, stringsAsFactors=FALSE)
colnames(coldata)
rownames(coldata)
coldata$Age_months <- as.factor(coldata$Age_months)
coldata$Genotype <- as.factor(coldata$Genotype)
coldata$Genotype <- relevel(coldata$Genotype, "WT")
levels(coldata$Genotype)

# get methylation data
load(file="/mnt/data1/Thea/IsabelSamples/data/J20GenotypeMatrixFiltered.Rdata") # object: RRBS
methylation <- RRBS[,2:64]
rownames(methylation) <- RRBS[,1]
colnames(methylation) <- sub("_m.*", "", colnames(methylation))
identical(colnames(methylation), rownames(coldata))
save(methylation, file = "/mnt/data1/isabel/RRBS/J20methylation.RData")

# get gene list
load("/mnt/data1/isabel/RRBS/J20GenotypeStatsManhattan.RData") # object: genotype_stats

threshold <- 0.05
sig_genotype <- genotype_stats[which(genotype_stats[,"FDR"]<threshold),]
write.csv(sig_genotype, file = "J20_RRBS_sig_genotype.csv")
#genotype_stats_FDR <- read.csv("/mnt/data1/isabel/RRBS/J20GenotypeStats.csv", row.names=1, stringsAsFactors=FALSE)

# code colours
ann_colors = list(
	Genotype=c(WT="black", TG=color_J20_TG),
	Age_months=c("6"=age1, "8"=age2, "10"=age3, "12"=age4))

## Genotype
# clustered, methylation
sites_to_plot <- rownames(sig_genotype)
plotdata <- methylation[sites_to_plot,]
plotdata <- plotdata[complete.cases(plotdata), ]

pdf(file = "/mnt/data1/isabel/RRBS/J20_genotype_heatmapclustered_methylation.pdf", width=10, height=15)
pheatmap(
	plotdata[, order(coldata$Genotype, coldata$Age_months)],
	color = blueWhiteRed(100)[0:100],
	cluster_rows=TRUE,
	cluster_cols=TRUE,
	annotation_col=coldata[ order(coldata$Genotype, coldata$Age_months),c("Genotype", "Age_months")],
	annotation_colors=ann_colors,
	show_colnames=FALSE,
	fontsize_row = 5,
	show_rownames = FALSE,
	scale="row")
dev.off()