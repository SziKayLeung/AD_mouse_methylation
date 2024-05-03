# title: RRBS - BiSeq J20 - Differentially Methylated Positions (DMPs) - ~ Pathology
# author: Isabel Castanho (I.S.Castanho@exeter.ac.uk)
# date: 26/08/2020

# Setup
setwd("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/")

library(betareg)
# library(cgwtools) # package that allows to resave objects in R using resave(..., list = character(), file)

# parallel computing
library(parallel)
numcores <- detectCores()

color_J20_TG <- "#FF5A62"

source("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/code/RRBS-GitLab/BetaRegressionDMPs-Pathology.R")

# Phenotypic data
coldata <- read.csv("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/J20_coldata_RRBS.csv", row.names=1, stringsAsFactors=FALSE)
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
load(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/J20_RRBSbetasComplete.RData")
# RRBS_completebetas

# Beta regression ~Pathology
RRBS_completebetas_ECX <- RRBS_completebetas[,rownames(coldata_pathology)]

betaResults <- BetaRegressionRRBSpathology(betas = RRBS_completebetas_ECX,
                                           pheno = coldata_pathology,
                                           formula = ~ECX,
                                           link = "probit",
                                           num.cores = numcores)


betaResultsPathology <- betaResults
save(betaResultsPathology,
     file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/J20_betaResultsDMPsPathology.RData")
write.csv(betaResultsPathology,
          file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/J20_betaResultsDMPsPathology.csv")

# Change columns from factor to numeric
betaResults[1:3] <- lapply(betaResults[1:3], as.character)
betaResults[1:3] <- lapply(betaResults[1:3], as.numeric)

# Check and remove NAs
## check how many NAs there are in the dataframe
na_df <- as.data.frame(betaResults[rowSums(is.na(betaResults)) > 0,])
nrow(na_df) # number of sites that have no data (NA)

# remove NAs
betaResultsClean <- betaResults[complete.cases(betaResults),]

# Add FDR correction
FDR_adj_pathology <- p.adjust(betaResultsClean[,"p.val.Pathology"], method = "fdr")

stats_table <- cbind(betaResultsClean[,"Position"],
                     FDR_adj_pathology,
                     betaResultsClean[,c("p.val.Pathology", "estimate.Pathology", "std.error.Pathology")])

write.csv(stats_table,
          file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMPs/J20/DMPsPathology_stats_table_J20.csv")

sig_pathology <- stats_table[which(stats_table[,"FDR_adj_pathology"] < 0.05),]
nrow(sig_pathology)
write.csv(sig_pathology,
          file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMPs/J20/DMPsPathology_J20_sig_pathology.csv")

# session info
sessionInfo()

# R version 3.6.0 (2019-04-26)
# Platform: x86_64-redhat-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] doParallel_1.0.16 iterators_1.0.13  foreach_1.5.1     lmtest_0.9-38     zoo_1.8-8         betareg_3.1-3    
# 
# loaded via a namespace (and not attached):
#   [1] lattice_0.20-38   codetools_0.2-16  grid_3.6.0        stats4_3.6.0      flexmix_2.3-17    sandwich_3.0-0   
# [7] Formula_1.2-4     tools_3.6.0       compiler_3.6.0    nnet_7.3-12       modeltools_0.2-23