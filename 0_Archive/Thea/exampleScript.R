## Example script for data filtering for Isabel

## This script will explain how to use RRBSMatrix.R and RRBSFiltering.R

# First you need to upload and unzip the .cov data files. They are currently 
# on my hard drive but I will upload them to knight as they are finished. 
# They can be unzipped using gunzip *.cov.gz but I will do this too if I 
# get the chance. 

# The file path is "/mnt/data1/Thea/IsabelSamples/data/covFolder"


# load the functions I have created using
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
source("RRBSMatrix.R")
source("RRBSFiltering.R")

# You will need to run the following:
RRBSmatrix<-RRBSMatrix("/mnt/data1/Thea/EpigeneticClock/MouseEpigeneticClock/bismarkFiles/")

# load the phenotype data and extract the groups column
setwd("/mnt/data1/Thea/IsabelSamples/data/")
phenotype<-read.csv("Tg4510_samples_table.csv")
groupPheno<-phenotype[,6]

# you will need to enter your read depth threshold and the maximum 
# number of NA's that you will allow per site. I have suggested fairly 
# arbitrary values in the comments at the start of the RRBSFiltering.R
# script
RRBS<-RRBSFiltering(RRBSmatrix, RDThreshold = 10,  groupPheno, sampleMin = 6, 
                    densityPlots=TRUE,
                    SDThreshold = 5,
                    covCols=TRUE)
setwd("/mnt/data1/Thea/IsabelSamples/RScripts")
save(RRBS, file="RRBSrTg4510.Rdata")

# Now RRBS will be a matrix containing n+1 columns where the first is 
# chromosome number and start position and the rest are methylation
# values per sample. This should be a useable format for analysis.

# If you want to understand what the functions do I have included 
# descriptors at the beggining of each. Any questions feel free to 
# email me! Good luck with analysis :) 
