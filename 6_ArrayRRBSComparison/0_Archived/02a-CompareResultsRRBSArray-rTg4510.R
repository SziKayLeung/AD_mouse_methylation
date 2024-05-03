## Plots and analysis are in here. The data is stroed in a rdat but can also be reran by running the first script
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(ggrepel)
setwd("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/CompareTech/rTg4510/")
load("TechComparison_rTg4510.RData")


pdf("Genotype/pval.correlatiRon.pdf")
plot(GenotypeResults$A.p.val.Genotype, GenotypeResults$R.p.val, xlab = "Array pval", ylab = "RRBS pval")
dev.off()

df <- GenotypeResults[,c("A.p.val.Genotype", "R.p.val")]
colnames(df) <- c("Array", "RRBS")
df.mlt <- melt(df)
pdf("Genotype/pval.distribution.pdf")
ggplot(df.mlt, aes(value, fill = variable)) +
	geom_density(alpha = 0.5)
dev.off()


pdf("Genotype/sd.correlation.pdf")
plot(GenotypeResults$A.std.error.Genotype, GenotypeResults$R.std.error, 
	xlab = "Array std error", ylab = "RRBS std error")
dev.off()

df <- GenotypeResults[,c("A.std.error.Genotype", "R.std.error")]
colnames(df) <- c("Array", "RRBS")
df.mlt <- melt(df)
pdf("Genotype/sd.distribution.pdf")
ggplot(df.mlt, aes(value, fill = variable)) +
	geom_density(alpha = 0.5)
dev.off()


# Compare pseudo R sqrt
# Pseudo R2 is a measure of how well variables of the model explain some phenomenon

cor(GenotypeResults$A.pseudo.R.sqrt, GenotypeResults$R.pseudo.R.sqrt)

pdf("Genotype/R_sqrt_comparison.pdf")
# ggplot(GenotypeResults, aes(x = A.pseudo.R.sqrt, y = R.pseudo.R.sqrt)) +
# 	geom_point()
ggscatter(GenotypeResults, x = 'A.pseudo.R.sqrt', y = 'R.pseudo.R.sqrt',
          add = 'reg.line',
          add.params = list(color = 'blue', fill = 'lightgray'),
          conf.int = T) +
  stat_cor(method = 'pearson') +
  labs(title = "pseudo.R.sqrt",
       x = 'Array', y = 'RRBS') +
   theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Genotype/Log10_pvalCorr.pdf")
# plot(-log10(GenotypeResults$A.p.val.Genotype), -log10(GenotypeResults$R.p.val))
df <- cbind(-log10(GenotypeResults$A.p.val.Genotype), -log10(GenotypeResults$R.p.val))
df <- as.data.frame(df)
colnames(df) <- c("Array", "RRBS")
ggscatter(df, x = 'Array', y = 'RRBS',
          add = 'reg.line',
          add.params = list(color = 'blue', fill = 'lightgray'),
          conf.int = T) +
  stat_cor(method = 'pearson') +
  labs(title = "-log 10 Pvalues",
       x = 'Array', y = 'RRBS') +
   theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Genotype/EstCorr.pdf")
ggscatter(GenotypeResults, x = 'A.estimate.Genotype', y = 'R.estimate',
          add = 'reg.line',
          add.params = list(color = 'blue', fill = 'lightgray'),
          conf.int = T) +
  stat_cor(method = 'pearson') +
  labs(title = "Estimates",
       x = 'Array', y = 'RRBS') +
   theme(plot.title = element_text(hjust = 0.5))
dev.off()



#############################  Tabulate how many are signifcant  ################

A.SignGenotypeResults <- GenotypeResults[which(GenotypeResults$A.FDR_adj_genotype <= 0.05),]
R.SignGenotypeResults <- GenotypeResults[which(GenotypeResults$R.FDR_adj_genotype <= 0.05),]

A.SignInteractionResults <- InteractionResults[which(InteractionResults$A.FDR_adj_Interaction <= 0.05),]
R.SignInteractionResults <- InteractionResults[which(InteractionResults$R.FDR_adj_Interaction <= 0.05),]


A.SignPathologyResults <- PathologyResults[which(PathologyResults$A.FDR_adj_Pathology <= 0.05),]
R.SignPathologyResults <- PathologyResults[which(PathologyResults$R.FDR_adj_Pathology <= 0.05),]

Array <- c(nrow(A.SignGenotypeResults), nrow(A.SignInteractionResults), nrow(A.SignPathologyResults))
RRBS <- c(nrow(R.SignGenotypeResults), nrow(R.SignInteractionResults), nrow(R.SignPathologyResults))

Techtbl <- rbind(Array, RRBS)
Techtbl <- as.data.frame(Techtbl)
colnames(Techtbl) <- c("Genotype", "Interaction", "Pathology")
rownames(Techtbl) <- c("Array", "RRBS")

png("NumberofFDRSignicantSitesof1201.png")
p<-tableGrob(Techtbl)
grid.arrange(p)
dev.off()

range(GenotypeResults$R.FDR_adj_genotype)
range(InteractionResults$R.FDR_adj_Interaction)
range(PathologyResults$R.FDR_adj_Pathology)


Geno <- GenotypeResults[which(GenotypeResults$A.FDR_adj_genotype <= 0.05 & GenotypeResults$R.FDR_adj_genotype <= 0.05),]
Inte <- InteractionResults[which(InteractionResults$A.FDR_adj_Interaction <= 0.05 & InteractionResults$R.FDR_adj_Interaction <= 0.05),]
Path <- PathologyResults[which(PathologyResults$A.FDR_adj_Pathology <= 0.05 & PathologyResults$R.FDR_adj_Pathology <= 0.05),]

Sharedsig <- c(nrow(Geno), nrow(Inte), nrow(Path))
Techtbl <- as.data.frame(t(Sharedsig))
colnames(Techtbl) <- c("Genotype", "Interaction", "Pathology")
rownames(Techtbl) <- c("Shared Significant Sites")

png("NumberofFDRSharedSignicantSitesof1201.png")
p<-tableGrob(Techtbl)
grid.arrange(p)
dev.off()

### GEnotype

if(nrow(A.SignGenotypeResults) > 0){
  p <- A.SignGenotypeResults[,c("A.Gene_Symbol", "A.Chr", "A.Bp","A.FDR_adj_genotype", "A.estimate.Genotype")]
  pdf("Genotype/ArraySignificantSites.pdf", h = nrow(p))
  rownames(p) <- c(1:nrow(p))
  p<-tableGrob(p)
  grid.arrange(p)
  dev.off()
} else {
  print("No sig array sites")
}

if(nrow(R.SignGenotypeResults) > 0){
  p <- R.SignGenotypeResults[,c("R.chr", "R.pos","R.FDR_adj_genotype", "R.estimate")]
  pdf("Genotype/RRBSSignificantSites.pdf", h = nrow(p))
  rownames(p) <- c(1:nrow(p))
  p<-tableGrob(p)
  grid.arrange(p)
  dev.off()
} else {
  print("No sig RRBS sites")
}


###  Interaction 
if(nrow(A.SignInteractionResults) > 0){
  p <- A.SignInteractionResults[,c("A.Gene_Symbol", "A.Chr", "A.Bp","A.FDR_adj_Interaction", "A.estimate.Interaction")]
  pdf("Interaction/ArraySignificantSites.pdf", h = nrow(p))
  rownames(p) <- c(1:nrow(p))
  p<-tableGrob(p)
  grid.arrange(p)
  dev.off()
} else {
  print("No sig array sites")
}

if(nrow(R.SignInteractionResults) > 0){
  p <- R.SignInteractionResults[,c("R.chr", "R.pos","R.FDR_adj_Interaction", "R.estimate")]
  pdf("Interaction/RRBSSignificantSites.pdf", h = nrow(p))
  rownames(p) <- c(1:nrow(p))
  p<-tableGrob(p)
  grid.arrange(p)
  dev.off()
} else {
  print("No sig RRBS sites")
}


### Pathology

if(nrow(A.SignPathologyResults) > 0){
  p <- A.SignPathologyResults[,c("A.Gene_Symbol", "A.Chr", "A.Bp","A.FDR_adj_Pathology", "A.estimate.Pathology")]
  pdf("Pathology/ArraySignificantSites.pdf", h = 25)
  rownames(p) <- c(1:nrow(p))
  p<-tableGrob(p)
  grid.arrange(p)
  dev.off()
} else {
  print("No sig array sites")
}

if(nrow(R.SignPathologyResults) > 0){
  p <- R.SignPathologyResults[,c("R.chr", "R.pos","R.FDR_adj_Pathology", "R.estimate")]
  pdf("Pathology/RRBSSignificantSites.pdf", h = nrow(p))
  rownames(p) <- c(1:nrow(p))
  p<-tableGrob(p)
  grid.arrange(p)
  dev.off()
} else {
  print("No sig RRBS sites")
}


### Shared Significant sites

# Plot effect size with plots having positions

# Pathology only
Path$A.FDR_adj_Pathology <-format(Path$A.FDR_adj_Pathology, scientific=F)
pdf("Pathology/SharedSigSites.pdf")
ggplot(Path, aes( x = Path$A.estimate.Pathology, y = Path$R.estimate)) +
  geom_point() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_classic() +
  labs(x = "Array", y = "RRBS", 
       title = "Effect Sizes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(aes(label =Path$A.Gene_Symbol),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50')

dev.off()
