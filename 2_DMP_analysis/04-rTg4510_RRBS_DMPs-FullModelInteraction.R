# title: RRBS - BiSeq rTg4510 - Differentially Methylated Positions (DMPs) - ~ Genotype + Age + Genotype*Age
# author: Isabel Castanho (I.S.Castanho@exeter.ac.uk)
# date: 26/08/2020

# Setup
setwd("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/")

library(betareg)
# library(cgwtools) # package that allows to resave objects in R using resave(..., list = character(), file)

# parallel computing
library(parallel)
numcores <- detectCores()

color_Tg4510_TG <- "#00AEC9"

source("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/code/RRBS-GitLab/BetaRegressionDMPs.R")

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


# Beta regression ~Genotype + Age + Genotype*Age
betaResults <- BetaRegressionRRBSinteraction(betas = RRBS_completebetas,
                                             pheno = coldata,
                                             formula = ~Genotype + Age_months + Genotype*Age_months,
                                             formulaNull = ~Genotype + Age_months,
                                             link = "probit",
                                             num.cores = numcores)

betaResultsInteraction <- betaResults
save(betaResultsInteraction,
     file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/rTg4510_betaResultsDMPsInteraction.RData")
write.csv(betaResultsInteraction,
          file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/rTg4510_betaResultsDMPsInteraction.csv")

# Change columns from factor to numeric
betaResults[1:14] <- lapply(betaResults[1:14], as.character)
betaResults[1:14] <- lapply(betaResults[1:14], as.numeric)

# Check and remove NAs
## check how many NAs there are in the dataframe
na_df <- as.data.frame(betaResults[rowSums(is.na(betaResults)) > 0,]) 
nrow(na_df) # number of sites that have no data (NA)

# remove NAs
betaResultsClean <- betaResults[complete.cases(betaResults),]

# Add FDR correction
FDR_adj_genotype <- p.adjust(betaResultsClean[,"p.val.Genotype"], method = "fdr")
FDR_adj_age <- p.adjust(betaResultsClean[,"p.val.Age"], method = "fdr")
FDR_adj_interaction <- p.adjust(betaResultsClean[,"p.val.Interaction"], method = "fdr")

stats_table <- cbind(betaResultsClean[,"Position"],
                     FDR_adj_genotype,
                     betaResultsClean[,c("p.val.Genotype", "meth.group1.WT", "meth.group2.TG",
                                         "meth.diff.Genotype", "estimate.Genotype", "std.error.Genotype", "pseudo.R.sqrt")],
                     FDR_adj_age,
                     betaResultsClean[,c("p.val.Age", "estimate.Age", "std.error.Age")],
                     FDR_adj_interaction,
                     betaResultsClean[,c("p.val.Interaction", "estimate.Interaction", "std.error.Interaction",
                                         "p.val.modelLRT")])

write.csv(stats_table,
          file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMPs/rTg4510/DMPsInteractionModel_stats_table_rTg4510.csv")

sig_genotype <- stats_table[which(stats_table[,"FDR_adj_genotype"] < 0.05),]
nrow(sig_genotype)
write.csv(sig_genotype,
          file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMPs/rTg4510/DMPsInteractionModel_rTg4510_sig_genotype.csv")

sig_age <- stats_table[which(stats_table[,"FDR_adj_age"] < 0.05),]
nrow(sig_age)
write.csv(sig_age,
          file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMPs/rTg4510/DMPsInteractionModel_rTg4510_sig_age.csv")

sig_interaction <- stats_table[which(stats_table[,"FDR_adj_interaction"] < 0.05),]
nrow(sig_interaction)
write.csv(sig_interaction,
          file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMPs/rTg4510/DMPsInteractionModel_rTg4510_sig_interaction.csv")

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
#   [1] doParallel_1.0.16           iterators_1.0.13            foreach_1.5.1               lmtest_0.9-38              
# [5] zoo_1.8-8                   betareg_3.1-3               BiSeq_1.26.0                Formula_1.2-4              
# [9] SummarizedExperiment_1.16.1 DelayedArray_0.12.3         BiocParallel_1.20.1         matrixStats_0.57.0         
# [13] Biobase_2.46.0              GenomicRanges_1.38.0        GenomeInfoDb_1.22.1         IRanges_2.20.2             
# [17] S4Vectors_0.24.4            BiocGenerics_0.32.0        
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.5               compiler_3.6.0           XVector_0.26.0           bitops_1.0-6             tools_3.6.0             
# [6] zlibbioc_1.32.0          digest_0.6.27            bit_4.0.4                memoise_1.1.0            RSQLite_2.2.1           
# [11] annotate_1.64.0          lattice_0.20-38          rlang_0.4.8              Matrix_1.2-17            DBI_1.1.0               
# [16] GenomeInfoDbData_1.2.2   rtracklayer_1.46.0       vctrs_0.3.5              Biostrings_2.54.0        bit64_4.0.5             
# [21] grid_3.6.0               nnet_7.3-12              globaltest_5.40.0        flexmix_2.3-17           AnnotationDbi_1.48.0    
# [26] survival_3.2-7           XML_3.99-0.3             lokern_1.1-8.1           blob_1.2.1               codetools_0.2-16        
# [31] splines_3.6.0            Rsamtools_2.2.3          modeltools_0.2-23        sfsmisc_1.1-7            GenomicAlignments_1.22.1
# [36] xtable_1.8-4             sandwich_3.0-0           RCurl_1.98-1.2 