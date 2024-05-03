# title: "RRBS - BiSeq rTg4510 - Script to merge smoothed DNA methylation (Biseq) and raw methylation, so we don't lose sites that were not asigned to clusters"
# author: "Isabel Castanho"
# date: 25/08/2020

# load packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("BiSeq")

library(BiSeq)

setwd("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/")

# Get raw methylation matrix
load(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/biseq_rTg4510.RData")
# rrbs

# Filter raw methylation data by coverage
dim(rrbs) # 4,824,627 sites
rrbs.reduced <- filterBySharedRegions(object=rrbs, groups=colData(rrbs)$Genotype, perc.samples=0.5, minCov=10)

rrbs.rel <- rawToRel(rrbs.reduced) # convert to methylation values
RRBS_rawbetas <- as.data.frame(methLevel(rrbs.rel))
nrow(RRBS_rawbetas) # 1,229,351 sites
coordinates <- as.data.frame(rrbs.reduced@rowRanges)
rownames(RRBS_rawbetas) <- paste(coordinates$seqnames, ":", coordinates$start, sep="")

# Get smooth methylation matrix
RRBS_smoothbetas <- as.data.frame(methLevel(predictedMeth))
nrow(RRBS_smoothbetas) # 1,761,997 sites
coordinates <- as.data.frame(predictedMeth@rowRanges)
rownames(RRBS_smoothbetas) <- paste(coordinates$seqnames, ":", coordinates$start, sep="")

# Replace raw methylation by smoothed methylation where possible
## rows that I want to replace
rowsToBeReplaced <- rownames(RRBS_smoothbetas)

## delete rows from raw data frame
RRBS_rawbetas_reduced <- RRBS_rawbetas[-which(rownames(RRBS_rawbetas) %in% rowsToBeReplaced),]
nrow(RRBS_rawbetas_reduced) # 180,470 sites

## add rows from smooth data frame
RRBS_completebetas <- rbind(RRBS_rawbetas_reduced, RRBS_smoothbetas)
## order by genomic position (row names)
RRBS_completebetas <- RRBS_completebetas[order(rownames(RRBS_completebetas)),]
nrow(RRBS_completebetas) # 1,942,467 sites
# 713,116 sites added to clean raw data
# 180,470 sites added to smooth data

save(RRBS_completebetas,
     file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/rTg4510_RRBSbetasComplete.RData")

# Packages & versions
sessionInfo()
# R version 3.6.0 (2019-04-26)
# BiSeq_1.26.0