# J20 Analysis with more relaxed NA filtering

# same structure as before

# Dorothea Seiler Vellame 08-11-2018

# load the functions I have created using
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("RRBSFiltering.R")

setwd("/mnt/data1/Thea/IsabelSamples/data")
load("J20MatrixUnfiltered.Rdata")


# Create phenotype that is all 1 because I don't want to seperate by group
groupPheno<-as.factor(rep(1,(ncol(RRBSmatrix)-1)/2))


# Filter the data 
RRBS<-RRBSFiltering(RRBSmatrix, RDThreshold = 10,  groupPheno, sampleMin = 6, 
                    densityPlots=TRUE,
                    pdfFileName = "J20RelaxedDensityPlots.pdf",
                    SDThreshold = 5,
                    covCols=FALSE)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBS, file="J20MatrixRelaxedFiltered.Rdata")


# load the function
setwd("/mnt/data1/Thea/IsabelSamples/RScripts/")
source("staticticalTest.R")

# load filtered data
setwd("/mnt/data1/Thea/IsabelSamples/data/")
load("J20MatrixRelaxedFiltered.Rdata")

# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("J20_phenotype_RRBS.csv")
groupPheno<-phenotype[,c(4,5)]

# run analysis
J20stats<-statisticalTest(RRBS,groupPheno)
# save results
setwd("/mnt/data1/Thea/IsabelSamples/data")
save(J20stats,file="J20Relaxedstats.Rdata")
