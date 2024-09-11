# J20 analysis

# Dorothea Seiler Vellame 22-10-2018

# load the functions I have created using
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("RRBSMatrix.R")
source("RRBSFiltering.R")

# Combine J20 samples into one matrix
#RRBSmatrix<-RRBSMatrix("/mnt/data1/Thea/IsabelSamples/data/covFolder")

# save intermediate matrix
#setwd("/mnt/data1/Thea/IsabelSamples/data/")
#save(RRBSmatrix,file="J20MatrixUnfiltered.Rdata")

# load the saved matrix
setwd("/mnt/data1/Thea/IsabelSamples/data/")
load("J20MatrixUnfiltered.Rdata")


# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("J20_phenotype_RRBS.csv")
groupPheno<-phenotype[,6]


# Filter the data 
RRBS<-RRBSFiltering(RRBSmatrix, RDThreshold = 10,  groupPheno, sampleMin = 6, 
                    densityPlots=TRUE,
                    pdfFileName = "J20DensityPlots.pdf",
                    SDThreshold = 5,
                    covCols=TRUE)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBS, file="J20MatrixFiltered.Rdata")

## carry out statistical analysis

# load the function
setwd("/mnt/data1/Thea/IsabelSamples/RScripts/")
source("staticticalTest.R")

# load filtered data
setwd("/mnt/data1/Thea/IsabelSamples/data/")
load("J20MatrixFiltered.Rdata")

# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("J20_phenotype_RRBS.csv")
groupPheno<-phenotype[,c(4,5)]

# run analysis
J20stats<-statisticalTest(RRBS,groupPheno)
# save results
save(J20stats,file="J20stats.Rdata")