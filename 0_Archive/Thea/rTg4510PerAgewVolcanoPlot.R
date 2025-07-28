# Analysis of rTg4510 WT vs TG per age group 
# create volcano plots for each

# Dorothea Seiler Vellame 10-01-2019

# load in the filtered data, to split by WT and TG
setwd("/mnt/data1/Thea/IsabelSamples/data")
load("rTg4510MatrixFiltered.Rdata")
RRBSAge<-RRBS[,1:63]
rm(RRBS)

# load the phenotype data 
phenotype<-read.csv("rTg4510_phenotype_RRBS.csv")
pheno<-phenotype[,c(4,5)]
rm(phenotype)


# split the data into each age group and make sub pheno files
ages<-levels(as.factor(pheno[,2]))

RRBSm<-RRBSAge[,2:63]

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

rTg4510stats1<-statisticalTestTwoGroups(RRBS1,pheno1)
rTg4510stats2<-statisticalTestTwoGroups(RRBS2,pheno2)
rTg4510stats3<-statisticalTestTwoGroups(RRBS3,pheno3)
rTg4510stats4<-statisticalTestTwoGroups(RRBS4,pheno4)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(rTg4510stats1,rTg4510stats2,rTg4510stats3,rTg4510stats4,file="rTg4510GenotypePerAge.Rdata")

# # create volcano plots
# library(ggplot2)
# library(limma)
#
# plot(rTg4510stats1[,3],rTg4510stats1[,5])
# 
# setwd("/mnt/data1/Thea/IsabelSamples/data")
# pdf(file="rTg4510VolcanoPlot.pdf")
# 
# fit<-lmFit(RRBS1[,2:ncol(RRBS1)],as.factor(pheno1))
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
