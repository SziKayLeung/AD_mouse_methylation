# title: RRBS - BiSeq rTg4510 - Methylation difference between mean for TG in latest time point (T4) and all WT
## Pathology results
# author: Isabel Castanho (I.S.Castanho@exeter.ac.uk)
# date: 02/12/2020

# Setup
### set up parameters to run parallel
library(foreach)
library(doParallel) # set up parameters to run parallel
num.cores <- detectCores()
cl<-makeCluster(num.cores)
registerDoParallel(cl)
clusterEvalQ(cl, library(foreach))

setwd("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/")

# Phenotypic data
coldata <- read.csv("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/Tg4510_coldata_RRBS.csv", row.names=1, stringsAsFactors=FALSE)
coldata$Age_months <- as.numeric(coldata$Age_months)
coldata$Genotype <- as.factor(coldata$Genotype)
coldata$Genotype <- relevel(coldata$Genotype, "WT")
levels(coldata$Genotype)

# ENTORHINAL CORTEX
class(coldata$ECX)

# Remove samples that have NA for pathology
coldata_ECX <- coldata[,c("Genotype", "Age_months", "Histology_no", "ECX")]
coldata_ECX_clean <- na.omit(coldata_ECX)
coldata_pathology <- coldata_ECX_clean

# Load the data
load(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/rTg4510_RRBSbetasComplete.RData")
# RRBS_completebetas
RRBS_completebetas_ECX <- RRBS_completebetas[,rownames(coldata_pathology)]

# get methylation data
dat <- RRBS_completebetas_ECX

identical(rownames(coldata_pathology), colnames(dat))

age <- "8"

# Create function which performs analysis for each row
testRow <- function(row){
  ### Initialize dataframe where means and mean difference will be stored
  newDat <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(newDat) <- c("mean_allWT", "n WT", "mean_T4TG", "n T4 TG", "meanDiffT4")

  ### Calculate difference between mean for WT and mean for TG in latest time point (T4)
  site <- rownames(dat)[i]
  dat1<- cbind(coldata_pathology[,c("Age_months","Genotype")], t(dat[site,]))
  mean_allWT <- mean(subset(dat1, subset=(dat1$Genotype=="WT"))[,site], na.rm=TRUE)
  dat1_NAomit <- na.omit(dat1)
  n_WT <- length(subset(dat1_NAomit, subset=(dat1_NAomit$Genotype=="WT"))[,site])
  mean_T4TG <- mean(subset(dat1, subset=(dat1$Age_months==age & dat1$Genotype=="TG"))[,site], na.rm=TRUE)
  n_T4_Tg <- length(subset(dat1_NAomit, subset=(dat1_NAomit$Age_months==age & dat1_NAomit$Genotype=="TG"))[,site])
  meanDiffT4 <- mean_allWT - mean_T4TG

  ### save row in final dataframe
  newDat[1,1] <- mean_allWT
  newDat[1,2] <- n_WT
  newDat[1,3] <- mean_T4TG
  newDat[1,4] <- n_T4_Tg
  newDat[1,5] <- meanDiffT4

  newList <- list(newDat)
  return(newList)
}

# Run function to calculate means and mean difference
res <- foreach(i=1:nrow(dat)) %dopar% {
  testRow(row = dat[i,])
}

# Save results as a dataframe because res is saved as list (from foreach dopar)
DatList <- lapply(res, `[[`, 1)
DatFrame <-  as.data.frame(do.call(rbind, DatList))
rownames(DatFrame) <- rownames(dat)

save(DatFrame,
     file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/rTg4510_DMPsPathology_meanDiffTG-WT.RData")

stopCluster(cl)


# session info
sessionInfo()