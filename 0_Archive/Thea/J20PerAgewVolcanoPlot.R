# Analysis of J20 WT vs TG per age group 
# create volcano plots for each

# Dorothea Seiler Vellame 10-01-2019

# load in the filtered data, to split by WT and TG
setwd("/mnt/data1/Thea/IsabelSamples/data")
load("J20AgeMatrixFiltered.Rdata")
RRBSAge<-RRBS[,1:64]
rm(RRBS)

# load the phenotype data 
phenotype<-read.csv("J20_phenotype_RRBS.csv")
pheno<-phenotype[,c(4,5)]
rm(phenotype)


# split the data into each age group and make sub pheno files
ages<-levels(as.factor(pheno[,2]))

RRBSm<-RRBSAge[2:64]

RRBS1<-cbind(RRBSAge[,1],RRBSm[,pheno[,2]==ages[1]])
RRBS2<-cbind(RRBSAge[,1],RRBSm[,pheno[,2]==ages[2]])
RRBS3<-cbind(RRBSAge[,1],RRBSm[,pheno[,2]==ages[3]])
RRBS4<-cbind(RRBSAge[,1],RRBSm[,pheno[,2]==ages[4]])

# refilter each group by standard deviation
SDThreshold<-5

sdMatrix<-apply(RRBS1[,2:ncol(RRBS1)],1,sd,na.rm=TRUE)
k<-which(sdMatrix>SDThreshold)
RRBS1<-RRBS1[k,]

sdMatrix<-apply(RRBS2[,2:ncol(RRBS2)],1,sd,na.rm=TRUE)
k<-which(sdMatrix>SDThreshold)
RRBS2<-RRBS2[k,]

sdMatrix<-apply(RRBS3[,2:ncol(RRBS3)],1,sd,na.rm=TRUE)
k<-which(sdMatrix>SDThreshold)
RRBS3<-RRBS3[k,]

sdMatrix<-apply(RRBS4[,2:ncol(RRBS4)],1,sd,na.rm=TRUE)
k<-which(sdMatrix>SDThreshold)
RRBS4<-RRBS4[k,]

pheno1<-pheno[pheno[,2]==ages[1],1]
pheno2<-pheno[pheno[,2]==ages[2],1]
pheno3<-pheno[pheno[,2]==ages[3],1]
pheno4<-pheno[pheno[,2]==ages[4],1]


setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("statisticalTestTwoGroups.R")

J20stats1<-statisticalTestTwoGroups(RRBS1,pheno1)
J20stats2<-statisticalTestTwoGroups(RRBS2,pheno2)
J20stats3<-statisticalTestTwoGroups(RRBS3,pheno3)
J20stats4<-statisticalTestTwoGroups(RRBS4,pheno4)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(J20stats1,J20stats2,J20stats3,J20stats4,file="J20GenotypePerAge.Rdata")


# plot(log2(stats1[,3]),log10(stats1[,4]),xlab = "log2 fold change", ylab = "log10 p value")
# 
# 
# # create volcano plots
# library(ggplot2)
# library(limma)
# 
# setwd("/mnt/data1/Thea/IsabelSamples/data")
# pdf(file="J20VolcanoPlot.pdf")
# 
# fit<-lmFit(lapply(RRBS1[,2:ncol(RRBS1)],as.numeric),as.factor(pheno1))
# eb<-ebayes(fit)
# fc = fit$coef
# sig = -log10(eb$p.value)
# df <- data.frame(fc, sig)
# df$thre <- as.factor(abs(fc) < 0.4 & sig < -log10(0.05))
# 
# ggplot(data=df, aes(x=fc, y = sig, color=thre)) +
#   geom_point(alpha=.6, size=1.2) +
#   theme(legend.position="none") +
#   xlab("Effect size") +
#   ylab("-log10 p value")
# 
# 
# fit<-lmFit(RRBS2[,2:ncol(RRBS2)],as.factor(pheno2))
# eb<-ebayes(fit)
# fc = fit$coef
# sig = -log10(eb$p.value)
# df <- data.frame(fc, sig)
# df$thre <- as.factor(abs(fc) < 0.4 & sig < -log10(0.05))
# 
# ggplot(data=df, aes(x=fc, y = sig, color=thre)) +
#   geom_point(alpha=.6, size=1.2) +
#   theme(legend.position="none") +
#   xlab("Effect size") +
#   ylab("-log10 p value")
# 
# 
# 
# fit<-lmFit(RRBS3[,2:ncol(RRBS3)],as.factor(pheno3))
# eb<-ebayes(fit)
# fc = fit$coef
# sig = -log10(eb$p.value)
# df <- data.frame(fc, sig)
# df$thre <- as.factor(abs(fc) < 0.4 & sig < -log10(0.05))
# 
# ggplot(data=df, aes(x=fc, y = sig, color=thre)) +
#   geom_point(alpha=.6, size=1.2) +
#   theme(legend.position="none") +
#   xlab("Effect size") +
#   ylab("-log10 p value")
# 
# 
# 
# fit<-lmFit(RRBS4[,2:ncol(RRBS4)],as.factor(pheno4))
# eb<-ebayes(fit)
# fc = fit$coef
# sig = -log10(eb$p.value)
# df <- data.frame(fc, sig)
# df$thre <- as.factor(abs(fc) < 0.4 & sig < -log10(0.05))
# 
# ggplot(data=df, aes(x=fc, y = sig, color=thre)) +
#   geom_point(alpha=.6, size=1.2) +
#   theme(legend.position="none") +
#   xlab("Effect size") +
#   ylab("-log10 p value")
# 
# 
# dev.off()
