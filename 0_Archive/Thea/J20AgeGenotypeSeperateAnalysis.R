# J20 analysis by age and phenotype seperately

# Dorothea Seiler Vellame 07-01-2019

# load the functions I have created using
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("RRBSFiltering.R")

# load the saved matrix
setwd("/mnt/data1/Thea/IsabelSamples/data/")
load("J20MatrixUnfiltered.Rdata")

### ANALYSE BY AGE ONLY ####################################################

# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("J20_phenotype_RRBS.csv")
groupPheno<-as.factor(phenotype[,5])   # age column of phenotype


# Filter the data 
RRBS<-RRBSFiltering(RRBSmatrix, RDThreshold = 10,  groupPheno, sampleMin = 5, 
                    densityPlots=TRUE,
                    pdfFileName = "J20AgeDensityPlots.pdf",
                    SDThreshold = 5,
                    covCols=TRUE)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBS, file="J20AgeMatrixFiltered.Rdata")

# load the data
load("J20AgeMatrixFiltered.Rdata")
RRBSAge<-RRBS[,1:64]

## load statistical test function
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("statisticalTestPerPheno.R")

RRBSAgeStat<-statisticalTestPerPheno(RRBSAge,groupPheno)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBSAgeStat, file="J20AgeStats.Rdata")

### ANALYSE BY GENOTYPE ####################################################
## load statistical test function
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("statisticalTestTwoGroups.R")


# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("J20_phenotype_RRBS.csv")
groupPheno<-as.factor(phenotype[,4])   # genotype column of phenotype


# Filter the data 
RRBS<-RRBSFiltering(RRBSmatrix, RDThreshold = 10,  groupPheno, sampleMin = 5, 
                    densityPlots=TRUE,
                    pdfFileName = "J20GenotypeDensityPlots.pdf",
                    SDThreshold = 5,
                    covCols=TRUE)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBS, file="J20GenotypeMatrixFiltered.Rdata")


# load the data
load("J20GenotypeMatrixFiltered.Rdata")
RRBSGenotype<-RRBS[,1:64]

## load statistical test function
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("statisticalTestPerPheno.R")

RRBSGenotypeStat<-statisticalTestTwoGroups(RRBSGenotype,groupPheno)

setwd("/mnt/data1/Thea/IsabelSamples/data")
save(RRBSGenotypeStat, file="J20GenotypeStats.Rdata")