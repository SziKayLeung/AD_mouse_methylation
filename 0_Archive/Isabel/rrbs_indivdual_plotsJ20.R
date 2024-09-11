#Isabel Castanho I.Castanho@exeter.ac.uk

# RRBS individual plots

# R 3.4.3

setwd("/mnt/data1/isabel/RRBS/")

library("DESeq2")
library("pheatmap")
library("WGCNA")
#library(vioplot)

color_J20_TG <- "#FF5A62"

######################################################## J20 ########################################################
# column data / phenotypes
coldata <- read.csv("/mnt/data1/isabel/RRBS/J20_phenotype_RRBS.csv", row.names=1, stringsAsFactors=FALSE)
colnames(coldata)
rownames(coldata)
coldata$Age_months <- as.factor(coldata$Age_months)
coldata$Genotype <- as.factor(coldata$Genotype)
coldata$Genotype <- relevel(coldata$Genotype, "WT")
levels(coldata$Genotype)
coldata$Group_ID <- as.factor(coldata$Group_ID)
coldata$Group_ID <- factor(coldata$Group_ID, levels =c("WT_06m", "WT_08m", "WT_10m", "WT_12m", "TG_06m", "TG_08m", "TG_10m", "TG_12m"))
levels(coldata$Group_ID)

load(file = "/mnt/data1/isabel/RRBS/J20methylation.RData") # object: methylation

## Top sites
sites_to_plot <- c("chr16_45667595", "chr14_80916197", "chr6_65870804", "chr4_115869229", "chr16_45740128",
	"chr18_75831039", "chr13_54346100", "chr16_43973665", "chr16_20701672", "chr12_105535344")
genes_to_plot <- c("Tmprss7", "Gm9578", "Prdm5", "Mknk1", "Abhd10", 
	"Zbtb7c", "Cplx2", "Zdhhc23", "Fam131a", "Bdkrb2")

## Common top sites


## FAD genes / MAPT/ GWAS genes
# App, Psen1, Psen2, Mapt, Apoe, Clu, Picalm, Crry (mouse: Crry; humans: CR1 + CR2 -> CR1 in GWAS), Bin1, Cd2ap, Cd33, Epha1, Abca7, Trem2, Abi3, Pld3, Adam10, Frmd4a, Sorl1, Ptk2b
# GWAS: Apoe, Clu, Crry, Bin1, Cd2ap, Cd33, Epha1, Abca7, Adam10
# Rare variants: Trem2, Abi3, Pld3
sites_to_plot <- c("chr16_84922199", "chr16_84995182", "chr16_85084104", "chr16_85086498", "chr16_85087931", 
	"chr16_85120676", "chr11_104310400", "chr11_104316397", "chr14_65975708", "chr14_65981137", 
	"chr14_65981630", "chr18_32392554", "chr18_32392717", "chr18_32392818", "chr10_80015045", 
	"chr2_4390183", "chr9_41994042", "chr9_41996600", "chr9_42111450", "chr9_42158913",
	"chr14_66249810")
genes_to_plot <- c("App", "App", "App", "App", "App", 
	"App", "Mapt", "Mapt", "Clu", "Clu", 
	"Clu", "Bin1", "Bin1", "Bin1", "Abca7", 
	"Frmd4a", "Sorl1", "Sorl1", "Sorl1", "Sorl1",
	"Ptk2b")

## EWAS genes
# ANK1, RPL13, CDH23, RHBDF2
# TIMM22–ABR, ACTR3BP2, CLYBL–TM9SF2, SLC2A1–FLJ32224, COQ7–ITPRIPL2, HOXA, FOXK1–AP5Z1–RADIL, DIP2A, SERPINF1, SERPINF2
sites_to_plot <- c("chr8_22974676", "chr8_23034919", "chr8_23116360", "chr8_23123856", "chr8_23123864", 
	"chr8_23123900", "chr10_60303537", "chr10_60334144", "chr10_60346528", "chr10_60411179", 
	"chr10_60436753", "chr10_60495199", "chr10_60495544", "chr10_60540255", "chr10_60595936", 
	"chr10_60596311", "chr10_60638890", "chr10_60638910", "chr10_60649040", "chr11_76473889", 
	"chr11_76623004", "chr14_122179634", "chr14_122181645", "chr14_122200537", "chr14_122222613", 
	"chr14_122222644", "chr14_122222740", "chr14_122406878", "chr14_122407087", "chr6_52175647", 
	"chr6_52176077", "chr6_52190581", "chr6_52205277", "chr6_52208860", "chr6_52226257", 
	"chr6_52226267", "chr6_52233069", "chr6_52244694", "chr6_52248587", "chr6_52252821", 
	"chr6_52260482", "chr5_142435122", "chr5_142435276", "chr5_142466639", "chr5_142474598", 
	"chr5_142506674", "chr10_76295287", "chr11_75436401")
genes_to_plot <- c("Ank1", "Ank1", "Ank1", "Ank1", "Ank1", 
	"Ank1", "Cdh23", "Cdh23", "Cdh23", "Cdh23", 
	"Cdh23", "Cdh23", "Cdh23", "Cdh23", "Cdh23", 
	"Cdh23", "Cdh23", "Cdh23", "Cdh23", "Abr", 
	"Abr", "Clybl", "Clybl", "Clybl", "Clybl", 
	"Clybl", "Clybl", "Clybl", "Clybl", "Hoxaas3", 
	"Hoxaas3", "Hoxaas3", "Hoxaas3", "Hoxaas3", "Hoxa9", 
	"Hoxa9", "Hoxa10", "Hoxa11os", "Hoxa11os", "Hoxa11os", 
	"Hoxa13", "Foxk1", "Foxk1", "Ap5z1", "Ap5z1", 
	"Radil", "Dip2a", "Serpinf2")

## Other genes of interest
# Bace1, Vgf, Fyn, Trpa1, Abca1, Snca, Pld1, Pld2, Bdnf, Pik3c3, Pld4, Kmo (Sarah's paper), Pim3 (Sarah's paper), Rgcc (Sarah's paper)
sites_to_plot <- c("chr5_137031226", "chr10_39473164", "chr10_39473182", "chr10_39552922", "chr4_53104360", 
	"chr4_53126897", "chr4_53177713", "chr11_70554034", "chr18_29731870", "chr1_175661570", 
	"chr15_88862524", "chr14_79285759", "chr14_79301670")
genes_to_plot <- c("Vgf", "Fyn", "Fyn", "Fyn", "Abca1", 
	"Abca1", "Abca1", "Pld2", "Pik3c3", "Kmo", 
	"Pim3", "Rgcc", "Rgcc")

# RNA-seq paper
# Ccdc80, Abca8a, Htr1a, Hspa5
# Cst7, Wdfy1, Grxcr2, Itgax, Ifitm1
# Gfap, Cd68, Itgax, Clec7a, Cst7
sites_to_plot <- c("chr16_45096248", "chr11_110056727", "chr11_110078353", "chr13_105371119", "chr13_105647323",
	"chr1_79763237", "chr7_128103121", "chr7_128136914", "chr7_128205783")
genes_to_plot <- c("Ccdc80", "Abca8a", "Abca8a", "Htr1a", "Htr1a", 
	"Wdfy1", "Itgam", "Itgax", "Itgad")

### INDIVDUAL PLOTS
data <- methylation
for (index in 1:length(sites_to_plot)) {
	site <- sites_to_plot[index]
	gene <- genes_to_plot[index]
	data_plot <- cbind(coldata[,c("Age_months","Genotype", "Group_ID")], t(data[site,]))
	data_plot <- data_plot[complete.cases(data_plot), ]
	pdf(file = paste("/mnt/data1/isabel/RRBS/individual_plots/J20/", gene, "_", site, ".pdf", sep=""), width=5, height=10)
	#dev.new(width=5, height=10)
	par(mar=c(5,5,5,3))
	plot(y=data[site,], x=coldata$Genotype,
		col = c("black", color_J20_TG)[coldata$Genotype],
		pch = 19,
		cex=2,
		xlim=c(0.5,2.5),
		ylim=c(0, 100),
		axes=FALSE, ann=FALSE)
	axis(side=1, at=1:2, labels=c("WT", "TG"), cex.axis=2)
	axis(2, cex.axis=2)
	title(main = paste(gene, " ", site), xlab = "Genotype", ylab = "Methylation (%)", cex.main=2, cex.lab=2)
	boxplot(t(data[site,]) ~ coldata$Genotype, 
		axes=FALSE, ann=FALSE, 
		add = TRUE, at = 1:2)
	dev.off()
}

### Correlations ## Scatter plots methylation ~ gene expression
#load(file = "/mnt/data1/isabel/RRBS/rTg4510methylation.RData") # object: methylation
load(file="/mnt/data1/isabel/RNA-seq/analysis_isabel_new/J20_counts/normalized_counts_J20.RData") # object: normalized_counts

# make sure I have the same samples in coldata, methylation and normalized_counts
identical(rownames(coldata), colnames(methylation)) # TRUE
colnames(normalized_counts)
common_samples <- intersect(rownames(coldata), colnames(normalized_counts))
coldata_common <- coldata[common_samples,c("Age_months","Genotype", "Group_ID")]
methylation_common <- methylation[,common_samples]
normalized_counts_common <- as.data.frame(normalized_counts[,common_samples])
identical(rownames(coldata_common), colnames(normalized_counts_common)) # TRUE
identical(colnames(methylation_common), colnames(normalized_counts_common)) # TRUE

# plots
sites_to_plot <- c("chr16_45096248", "chr11_110056727", "chr11_110078353", "chr13_105371119", "chr13_105647323",
	"chr1_79763237", "chr7_128136914")
genes_to_plot <- c("Ccdc80", "Abca8a", "Abca8a", "Htr1a", "Htr1a", 
	"Wdfy1", "Itgax")

for (index in 1:length(sites_to_plot)) {
	site <- sites_to_plot[index]
	gene <- genes_to_plot[index]
	data_plot <- cbind(coldata_common, t(methylation_common[site,]), t(normalized_counts_common[gene,]))
	data_plot <- data_plot[complete.cases(data_plot),]
	pdf(file = paste("/mnt/data1/isabel/RRBS/individual_plots/J20/", gene, "_", site, "correlation.pdf", sep=""), width=10, height=10.3)
	#dev.new(width=10, height=10)
	par(mar=c(6,6,5,3), las=1,bty="l") # bottom, left, top, and right
	verboseScatterplot(data_plot[,site],
		data_plot[,gene],
		col = "black",
		bg=c("black", color_J20_TG)[data_plot$Genotype],
		pch = 21,
		cex=3,
		axes=FALSE,
		xlab = "Methylation (%)",
		ylab = "Gene expression (normalised counts)",
		#xlim=c(-0.4,0.4),
		#ylim=c(0,40),
		mgp=c(4,1,0), cex.main=2, cex.lab=2)
	axis(side=1, cex.axis=2)
	axis(2, cex.axis=2)
	dev.off()
}