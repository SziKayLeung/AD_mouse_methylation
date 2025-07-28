# rTg4510 Analysis with more relaxed NA filtering - filtering at sampleMin=1

# same structure as before

# Dorothea Seiler Vellame 08-11-2018

# load the functions I have created using
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("RRBSFiltering.R")

setwd("/mnt/data1/Thea/IsabelSamples/data")
load("rTg4510MatrixUnfiltered.Rdata")


# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("rTg4510_phenotype_RRBS.csv")
groupPheno<-phenotype[,6]


# Filter the data, setting min number of samples per group to 1
RRBS<-RRBSFiltering(RRBSmatrix, RDThreshold = 10,  groupPheno, sampleMin = 1, 
                    densityPlots=TRUE,
                    pdfFileName = "rTg4510Relaxed1DensityPlots.pdf",
                    SDThreshold = 5,
                    covCols=FALSE)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBS, file="rTg4510MatrixRelaxed1Filtered.Rdata")


## carry out statistical analysis

# load the function
setwd("/mnt/data1/Thea/IsabelSamples/RScripts/")
source("staticticalTest.R")

# load filtered data
setwd("/mnt/data1/Thea/IsabelSamples/data/")
load("rTg4510MatrixRelaxed1Filtered.Rdata")

# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("rTg4510_phenotype_RRBS.csv")
groupPheno<-phenotype[,c(4,5)]

# run analysis
rTg4510stats<-statisticalTest(RRBS,groupPheno)
# save results
setwd("/mnt/data1/Thea/IsabelSamples/data")
save(rTg4510stats,file="rTg4510Relaxed1stats.Rdata")
