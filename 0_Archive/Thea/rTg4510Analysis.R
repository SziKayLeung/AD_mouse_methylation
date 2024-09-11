# rTg4510 analysis (same structure as J20)

# Dorothea Seiler Vellame 23-10-2018

# load the functions I have created using
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("RRBSMatrix.R")
source("RRBSFiltering.R")

# Combine rTg4510 samples into one matrix
#RRBSmatrix<-RRBSMatrix("/mnt/data1/Thea/IsabelSamples/data/covFolder")

# save intermediate matrix
#setwd("/mnt/data1/Thea/IsabelSamples/data/")
#save(RRBSmatrix,file="rTg4510MatrixUnfiltered.Rdata")

# load the saved matrix
setwd("/mnt/data1/Thea/IsabelSamples/data/")
load("rTg4510MatrixUnfiltered.Rdata")


# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("rTg4510_phenotype_RRBS.csv")
groupPheno<-phenotype[,6]


# Filter the data 
RRBS<-RRBSFiltering(RRBSmatrix, RDThreshold = 10,  groupPheno, sampleMin = 6, 
                    densityPlots=TRUE,
                    pdfFileName = "rTg4510DensityPlots.pdf",
                    SDThreshold = 5,
                    covCols=TRUE)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBS, file="rTg4510MatrixFiltered.Rdata")


## carry out statistical analysis

# load the function
setwd("/mnt/data1/Thea/IsabelSamples/RScripts/")
source("staticticalTest.R")

# load filtered data
setwd("/mnt/data1/Thea/IsabelSamples/data/")
load("rTg4510MatrixFiltered.Rdata")

# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("rTg4510_phenotype_RRBS.csv")
groupPheno<-phenotype[,c(4,5)]

# run analysis
rTg4510stats<-statisticalTest(RRBS,groupPheno)
# save results
setwd("/mnt/data1/Thea/IsabelSamples/data")
save(rTg4510stats,file="rTg4510stats.Rdata")


# plot the significant sites
setwd("/mnt/data1/Thea/IsabelSamples/data")
load("rTg4510stats.Rdata")
load("rTg4510MatrixFiltered.Rdata")
phenotype<-read.csv("rTg4510_phenotype_RRBS.csv")
groupPheno<-phenotype[,c(4,5)]

k<-which(RRBS[,1]==rownames(rTg4510stats[[1]])[1:3])
k<-c(40994,40995,40996)

# plot
par(mfrow=c(1,3))

plot(t(groupPheno[,2]),t(RRBS[k[1],2:63]),
     col=ifelse(groupPheno[,1]=="WT","blue","red"),
     xlab="Age",ylab="% methylation",pch=16,
     main=RRBS[k[1],1])

plot(t(groupPheno[,2]),t(RRBS[k[2],2:63]),
     col=ifelse(groupPheno[,1]=="WT","blue","red"),
     xlab="Age",ylab="% methylation",pch=16,
     main=RRBS[k[2],1])

plot(t(groupPheno[,2]),t(RRBS[k[3],2:63]),
     col=ifelse(groupPheno[,1]=="WT","blue","red"),
     xlab="Age",ylab="% methylation",pch=16,
     main=RRBS[k[3],1])

legend("topright",levels(groupPheno[,1]),col=c("blue","red"),pch=16,cex = 1.5)



