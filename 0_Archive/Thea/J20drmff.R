## Carry out dmrff on J20 samples per age group

# Dorothea Seiler Vellame 23-01-2019

# source the needed functions
setwd("/mnt/data1/Thea/IsabelSamples/dmrff/tests")
source("functions.r")

# load in the filtered data, to split by WT and TG
setwd("/mnt/data1/Thea/IsabelSamples/data")
load("J20MatrixFiltered.Rdata")
RRBSAge<-RRBS[,1:64]
rm(RRBS)

# remove sex and mitochondrial chromosomes
RRBSAge<-RRBSAge[-which(substring(RRBSAge[,1],1,4)== "chrM" | substring(RRBSAge[,1],1,4)== "chrX" | substring(RRBSAge[,1],1,4)== "chrY"),]

# load the phenotype data 
phenotype<-read.csv("J20_phenotype_RRBS.csv") # keep the full phenotype file

variable.name="Age_months"

# split the data into each age group and make sub pheno files
ages<-levels(as.factor(phenotype[,5]))

RRBSm<-RRBSAge[2:64]

# create an annotation file
annot<-strsplit(as.character(RRBSAge[,1]),"_")
annotation<-unlist(annot)
annot<-matrix(ncol=2,nrow=length(RRBSAge[,1]))
annot[,1]=annotation[seq(1,length(annotation),2)]
annot[,2]=as.numeric(annotation[seq(2,length(annotation),2)])
# keep only chromosome number and make numeric
annot[,1]<-as.numeric(substring(annot[,1],4))
colnames(annot)<-c("chr","pos")
annot<-as.data.frame(annot)



# create each methylation matrix
RRBS1<-RRBSm[,phenotype[,5]==ages[1]]
RRBS2<-RRBSm[,phenotype[,5]==ages[2]]
RRBS3<-RRBSm[,phenotype[,5]==ages[3]]
RRBS4<-RRBSm[,phenotype[,5]==ages[4]]

pheno1<-phenotype[phenotype[,5]==ages[1],]
pheno2<-phenotype[phenotype[,5]==ages[2],]
pheno3<-phenotype[phenotype[,5]==ages[3],]
pheno4<-phenotype[phenotype[,5]==ages[4],]


# load libraries
library(limma)
setwd("/mnt/data1/Thea/IsabelSamples/dmrff/tests")
source("functions.r")

setwd("/mnt/data1/Thea/IsabelSamples/dmrff/R")
source("dmrff.r")
source("candidates.r")
source("bumphunter.r")
source("msg.r")
source("stats.r")
source("shrink.r")
source("region-members.r")
source("pre.r")
source("meta.r")
source("knit-report.r")
source("ivwfe.r")
source("impute-matrix.r")
source("annotate.r")

library(parallel)

idx<-order(annot$chr,as.numeric(annot$pos))
methylation <- RRBS1[idx,]
annot <- annot[idx,]

design1<-model.matrix(~ pheno1[,4])

fit <- lmFit(RRBS1, design1)

fit <- eBayes(fit)


dmrs <- dmrff(fit$coef[,2],
              
              sqrt(fit$s2.post) * fit$stdev.unscaled[,2],
              
              fit$p.value[,2],
              
              as.matrix(RRBS1),
              
              annot$chr, 
              
              as.numeric(annot$pos))

