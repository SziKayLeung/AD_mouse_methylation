# rTg4510 analysis by age and phenotype seperately

# Dorothea Seiler Vellame 08-01-2019

# load the functions I have created using
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("RRBSFiltering.R")

# load the saved matrix
setwd("/mnt/data1/Thea/IsabelSamples/data/")
load("rTg4510MatrixUnfiltered.Rdata")

### ANALYSE BY AGE ONLY ####################################################

# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("rTg4510_phenotype_RRBS.csv")
groupPheno<-as.factor(phenotype[,5])   # age column of phenotype


# Filter the data 
RRBS<-RRBSFiltering(RRBSmatrix, RDThreshold = 10,  groupPheno, sampleMin = 5, 
                    densityPlots=TRUE,
                    pdfFileName = "rTg4510AgeDensityPlots.pdf",
                    SDThreshold = 5,
                    covCols=TRUE)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBS, file="rTg4510AgeMatrixFiltered.Rdata")

# load the data
load("rTg4510AgeMatrixFiltered.Rdata")
RRBSAge<-RRBS[,1:63]

## load statistical test function
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("statisticalTestPerPheno.R")

RRBSAgeStat<-statisticalTestPerPheno(RRBSAge,groupPheno)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBSAgeStat, file="rTg4510AgeStats.Rdata")

### ANALYSE BY GENOTYPE ####################################################
## load statistical test function
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("statisticalTestTwoGroups.R")

# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("rTg4510_phenotype_RRBS.csv")
groupPheno<-as.factor(phenotype[,4])   # genotype column of phenotype


# Filter the data 
RRBS<-RRBSFiltering(RRBSmatrix, RDThreshold = 10,  groupPheno, sampleMin = 5, 
                    densityPlots=TRUE,
                    pdfFileName = "rTg4510GenotypeDensityPlots.pdf",
                    SDThreshold = 5,
                    covCols=TRUE)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBS, file="rTg4510GenotypeMatrixFiltered.Rdata")


# load the data
load("rTg4510GenotypeMatrixFiltered.Rdata")
RRBSGenotype<-RRBS[,1:63]

## load statistical test function
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("statisticalTestPerPheno.R")

RRBSGenotypeStat<-statisticalTestTwoGroups(RRBSGenotype,groupPheno)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBSGenotypeStat, file="rTg4510GenotypeStats.Rdata")