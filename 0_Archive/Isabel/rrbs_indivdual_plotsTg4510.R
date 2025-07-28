#Isabel Castanho I.Castanho@exeter.ac.uk

# RRBS individual plots

# R 3.4.3

setwd("/mnt/data1/isabel/RRBS/")

library("DESeq2")
library("pheatmap")
library("WGCNA")
#library(vioplot)

color_Tg4510_TG <- "#00AEC9"

######################################################## Tg4510 ########################################################
# column data / phenotypes
coldata <- read.csv("/mnt/data1/isabel/RRBS/Tg4510_phenotype_RRBS.csv", row.names=1, stringsAsFactors=FALSE)
colnames(coldata)
rownames(coldata)
coldata$Age_months <- as.factor(coldata$Age_months)
coldata$Genotype <- as.factor(coldata$Genotype)
coldata$Genotype <- relevel(coldata$Genotype, "WT")
levels(coldata$Genotype)
coldata$Group_ID <- as.factor(coldata$Group_ID)
coldata$Group_ID  <- factor(coldata$Group_ID, levels =c("Tg4510_WT_2m", "Tg4510_WT_4m", "Tg4510_WT_6m", "Tg4510_WT_8m", "Tg4510_TG_2m", "Tg4510_TG_4m", "Tg4510_TG_6m", "Tg4510_TG_8m"))
levels(coldata$Group_ID)

load(file = "/mnt/data1/isabel/RRBS/rTg4510methylation.RData") # object: methylation

## Top sites
sites_to_plot <- c("chr5_87231917", "chr18_60917719", "chr18_60917731", "chr18_60917761", "chr12_113190326", 
	"chr12_113266015", "chr7_103488891", "chr12_119325190", "chr12_113199742", "chr1_118279473", 
	"chr12_113190240")
genes_to_plot <- c("Ugt2b37", "Arsi", "Arsi", "Arsi", "Tmem121", 
	"Ighe", "Olfr609", "Gm6768", "Tmem121", "Tsn", 
	"Tmem121")

## FAD genes / MAPT/ GWAS genes
# App, Psen1, Psen2, Mapt, Apoe, Clu, Picalm, Crry (mouse: Crry; humans: CR1 + CR2 -> CR1 in GWAS), Bin1, Cd2ap, Cd33, Epha1, Abca7, Trem2, Abi3, Pld3, Adam10, Frmd4a, Sorl1, Ptk2b
# GWAS: Apoe, Clu, Crry, Bin1, Cd2ap, Cd33, Epha1, Abca7, Adam10
# Rare variants: Trem2, Abi3, Pld3
sites_to_plot <- c("chr16_85088148", "chr16_85091216", "chr16_85091220", "chr16_85091240", "chr16_85091540", 
	"chr16_85269100", "chr11_104247451", "chr11_104299117", "chr10_80013732", "chr10_80013994", 
	"chr10_80003740", "chr2_4551149", "chr2_4025291", "chr9_41934692", "chr9_41956913", 
	"chr9_41981976", "chr9_42055354", "chr9_42000832", "chr14_66265779")
genes_to_plot <- c("App", "App", "App", "App", "App", 
	"App", "Mapt", "Mapt", "Abca7", "Abca7", 
	"Abca7", "Frmd4a", "Frmd4a", "Sorl1", "Sorl1", 
	"Sorl1", "Sorl1", "Sorl1", "Ptk2b")

## EWAS genes
# ANK1, RPL13, CDH23, RHBDF2
# TIMM22–ABR, ACTR3BP2, CLYBL–TM9SF2, SLC2A1–FLJ32224, COQ7–ITPRIPL2, HOXA, FOXK1–AP5Z1–RADIL, DIP2A, SERPINF1, SERPINF2
sites_to_plot <- c("chr8_22978636", "chr8_22978666", "chr8_23035040", "chr8_23035042", "chr8_23045092", 
	"chr8_23058852", "chr8_23143029", "chr8_23143041", "chr8_123103130", "chr10_60321298", 
	"chr10_60346482", "chr10_60455760", "chr10_60523066", "chr10_60634392", "chr10_60634564", 
	"chr10_60686029", "chr11_116625068", "chr11_76508111", "chr11_76508178", "chr11_76542804",
	"chr11_76570619", "chr11_76577662", "chr11_76624326", "chr11_76624371", "chr11_76639885",
	"chr14_122172522", "chr14_122172566", "chr14_122257470", "chr14_122307062", "chr14_122406730",
	"chr7_118501993", "chr6_52199921", "chr6_52209539", "chr6_52217163", "chr6_52224435",
	"chr6_52225800", "chr6_52226299", "chr6_52244601", "chr6_52252901", "chr6_52259995",
	"chr6_52261001", "chr5_142400345", "chr5_142423143", "chr5_142474921", "chr5_142477404", 
	"chr5_142485494", "chr5_142546740", "chr10_76344952", "chr11_75431357", "chr11_75431375")
genes_to_plot <- c("Ank1", "Ank1", "Ank1", "Ank1", "Ank1", 
	"Ank1", "Ank1", "Ank1", "Rpl13", "Cdh23",
	"Cdh23", "Cdh23", "Cdh23", "Cdh23", "Cdh23",
	"Cdh23", "Rhbdf2", "Abr", "Abr", "Abr", 
	"Abr", "Abr", "Abr", "Abr", "Abr",
	"Clybl", "Clybl", "Clybl", "Clybl", "Clybl",
	"Coq7", "Hoxaas3", "Hoxaas3", "Hoxa7", "Hoxa9",
	"Hoxa9", "Hoxa9", "Hoxa11os", "Hoxa11os", "Hoxa13",
	"Hoxa13", "Foxk1", "Foxk1", "Ap5z1", "Ap5z1", 
	"Radil", "Radil", "Dip2a", "Serpinf2", "Serpinf2")

## Other genes of interest
# Bace1, Vgf, Fyn, Trpa1, Abca1, Snca, Pld1, Pld2, Bdnf, Pik3c3, Pld4, Kmo (Sarah's paper), Pim3 (Sarah's paper), Rgcc (Sarah's paper)
sites_to_plot <- c("chr5_137032777", "chr10_39370129", "chr10_39379166", "chr10_39430902", "chr10_39514342", 
	"chr10_39524581", "chr1_14877744", "chr3_27938975", "chr3_27939234", "chr3_28026728", 
	"chr11_70557370", "chr2_109693653", "chr18_29733434", "chr14_79296238", "chr14_79314003")
genes_to_plot <- c("Vgf", "Fyn", "Fyn", "Fyn", "Fyn", 
	"Fyn", "Trpa1", "Pld1", "Pld1", "Pld1", 
	"Pld2", "Bdnf",	"Pik3c3", "Rgcc", "Rgcc")

## RNA-seq paper genes
# Car4, Gpr17, Blnk, Hspa5, Gfap, Cd68, Itgax, Clec7a, Cst7
# C1qb, Mpeg1, Tyrobp
# Atp9a, Faim2, Ppp2r1a 
# Dlgap3, Shank3, Epn1, and Fbxl16
sites_to_plot <- c("chr11_84939460", "chr19_40941266", "chr7_30399773", "chr2_168623002", "chr2_168628800", 
	"chr2_168642111", "chr2_168644715", "chr2_168653245", "chr2_168656080", "chr2_168669802",
	"chr15_99499835", "chr15_99509851", "chr15_99524876", "chr15_99527736", "chr15_99539592",
	"chr4_127164722", "chr4_127178013", "chr4_127195186", "chr4_127198988", "chr4_127236687", 
	"chr4_127237058", "chr15_89517489", "chr15_89524254", "chr15_89548024", "chr15_89548209",
	"chr15_89548683")
genes_to_plot <- c("Car4", "Blnk", "Tyrobp", "Atp9a", "Atp9a",
	"Atp9a", "Atp9a", "Atp9a", "Atp9a", "Atp9a",
	"Faim2", "Faim2", "Faim2", "Faim2", "Faim2",
	"Dlgap3", "Dlgap3", "Dlgap3", "Dlgap3", "Dlgap3",
	"Dlgap3", "Shank3", "Shank3", "Shank3", "Shank3", 
	"Shank3")

### INDIVDUAL PLOTS
data <- methylation
for (index in 1:length(sites_to_plot)) {
	site <- sites_to_plot[index]
	gene <- genes_to_plot[index]
	data_plot <- cbind(coldata[,c("Age_months","Genotype", "Group_ID")], t(data[site,]))
	data_plot <- data_plot[complete.cases(data_plot), ]
	pdf(file = paste("/mnt/data1/isabel/RRBS/individual_plots/Tg4510/", gene, "_", site, ".pdf", sep=""), width=5, height=10)
	#dev.new(width=5, height=10)
	par(mar=c(5,5,5,3))
	plot(y=data[site,], x=coldata$Genotype,
		col = c("black", color_Tg4510_TG)[coldata$Genotype],
		pch = 19,
		cex=2,
		xlim=c(0.5,2.5),
		ylim=c(0, 100),
		axes=FALSE, ann=FALSE)
	axis(side=1, at=1:2, labels=c("WT", "TG"), cex.axis=2)
	axis(2, cex.axis=2)
	title(main = paste(gene, " ", site), ylab = "Methylation (%)", cex.main=2, cex.lab=2)
	boxplot(t(data[site,]) ~ coldata$Genotype, 
		axes=FALSE, ann=FALSE, 
		add = TRUE, at = 1:2)
	dev.off()
}

### Correlations ## Scatter plots methylation ~ gene expression
#load(file = "/mnt/data1/isabel/RRBS/rTg4510methylation.RData") # object: methylation
load(file="/mnt/data1/isabel/RNA-seq/analysis_isabel_new/Tg4510_counts/normalized_counts_Tg4510.RData") # object: normalized_counts

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
sites_to_plot <- c("chr11_84939460", "chr19_40941266")
genes_to_plot <- c("Car4", "Blnk")

sites_to_plot <- c("chr16_85088148", "chr16_85091216", "chr16_85091220", "chr16_85091240", "chr16_85091540", "chr16_85269100")
genes_to_plot <- c("App", "App", "App", "App", "App", "App")

sites_to_plot <- c("chr11_104247451", "chr11_104299117")
genes_to_plot <- c("Mapt", "Mapt")

for (index in 1:length(sites_to_plot)) {
	site <- sites_to_plot[index]
	gene <- genes_to_plot[index]
	data_plot <- cbind(coldata_common, t(methylation_common[site,]), t(normalized_counts_common[gene,]))
	data_plot <- data_plot[complete.cases(data_plot),]
	pdf(file = paste("/mnt/data1/isabel/RRBS/individual_plots/Tg4510/", gene, "_", site, "correlation.pdf", sep=""), width=10, height=10.3)
	#dev.new(width=10, height=10.3)
	par(mar=c(6,6,5,3), las=1,bty="l") # bottom, left, top, and right
	verboseScatterplot(data_plot[,site],
		data_plot[,gene],
		col = "black",
		bg=c("black", color_Tg4510_TG)[data_plot$Genotype],
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

