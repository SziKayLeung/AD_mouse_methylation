# title: RRBS - BiSeq J20 - Methylation difference between mean for TG in latest time point (T4) and all WT
## Interaction (Genotype*Age) results
# author: Isabel Castanho (I.S.Castanho@exeter.ac.uk)
# date: 04/10/2021

# Setup
### set up parameters to run parallel
library(foreach)
library(doParallel) # set up parameters to run parallel
num.cores <- detectCores()
cl<-makeCluster(num.cores)
registerDoParallel(cl)
clusterEvalQ(cl, library(foreach))

setwd("/lustre/projects/Research_Project-191406/isabel/RRBS_new/")

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
# RRBS_completebetas_ECX <- RRBS_completebetas[,rownames(coldata_pathology)]

# get DMPs to make the process faster
DMPs <- read.csv("/lustre/projects/Research_Project-191406/isabel/RRBS_new/DMPs/J20/DMPsInteractionModel_J20_sig_interaction.csv",
                    row.names=1, stringsAsFactors=FALSE)

# get methylation data
dat <- RRBS_completebetas[rownames(DMPs),]

identical(rownames(coldata), colnames(dat))

age <- "12"

# Create function which performs analysis for each row
testRow <- function(row){
  ### Initialize dataframe where means and mean difference will be stored
  newDat <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(newDat) <- c("mean_allWT", "n WT", "mean_T4TG", "n T4 TG", "meanDiffT4")

  ### Calculate difference between mean for WT and mean for TG in latest time point (T4)
  site <- rownames(dat)[i]
  dat1<- cbind(coldata[,c("Age_months","Genotype")], t(dat[site,]))
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
     file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/J20_DMPsInteraction_meanDiffTG-WT.RData")

stopCluster(cl)


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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] doParallel_1.0.16           iterators_1.0.13            foreach_1.5.1               BiSeq_1.26.0               
# [5] Formula_1.2-4               SummarizedExperiment_1.16.1 DelayedArray_0.12.3         BiocParallel_1.20.1        
# [9] matrixStats_0.57.0          Biobase_2.46.0              GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        
# [13] IRanges_2.20.2              S4Vectors_0.24.4            BiocGenerics_0.32.0        
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.5               compiler_3.6.0           XVector_0.26.0           bitops_1.0-6            
# [5] tools_3.6.0              zlibbioc_1.32.0          digest_0.6.27            bit_4.0.4               
# [9] memoise_1.1.0            RSQLite_2.2.1            annotate_1.64.0          lattice_0.20-38         
# [13] rlang_0.4.8              Matrix_1.2-17            DBI_1.1.0                GenomeInfoDbData_1.2.2  
# [17] rtracklayer_1.46.0       vctrs_0.3.5              Biostrings_2.54.0        bit64_4.0.5             
# [21] lmtest_0.9-38            grid_3.6.0               nnet_7.3-12              globaltest_5.40.0       
# [25] flexmix_2.3-17           AnnotationDbi_1.48.0     survival_3.2-7           XML_3.99-0.3            
# [29] lokern_1.1-8.1           blob_1.2.1               codetools_0.2-16         splines_3.6.0           
# [33] Rsamtools_2.2.3          modeltools_0.2-23        sfsmisc_1.1-7            GenomicAlignments_1.22.1
# [37] xtable_1.8-4             betareg_3.1-3            sandwich_3.0-0           RCurl_1.98-1.2          
# [41] zoo_1.8-8  