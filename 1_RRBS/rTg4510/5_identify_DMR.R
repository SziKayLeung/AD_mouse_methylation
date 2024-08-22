#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## 17.06.2024: run DMR analysis on genotype results
## model: - ~ Genotype + Age + Genotype*Age
## --------------------------------

##--------------- packages ----------------

suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(BiSeq))
suppressMessages(library(betareg))


## ---------- import data -----------------

# directory input
rootDir = "/lustre/projects/Research_Project-191406"
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/"

source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "1_RRBS/functions/identifyDMR.R"))

## results from running BiSeq clustering (needed for cluster id)
load(file = paste0(dirnames$processed, "/rrbs/biseq_rTg4510.RData"))

## smoothed BiSeq RRBS dataset
rTg4510_rrbs_sig <- get(load(file = paste0(dirnames$annotated, "/rrbs/rTg4510_rrbs_annoSigResultsDMPs.RData")))


## ---------- datawrangle format for input -----------------

# extract cluster.id 
originalClusterID <- as.data.frame(GRangesobject)
originalClusterID <- originalClusterID %>% mutate(Position = paste0(seqnames,":",start))
retainedClusterID <- lapply(rTg4510_rrbs_sig, function(x) originalClusterID %>% filter(Position %in% x[["Position"]]))

# select column for downstream and retain original column names
betaColNames <- c("chr", "pos", "p.val", "meth.group1", "meth.group2", "meth.diff", "estimate", "std.error", "pseudo.R.sqrt", "cluster.id")  
rTg4510_rrbs_DMR_input <- list()
rTg4510_rrbs_DMR_input$Genotype <- 
  rTg4510_rrbs_sig$Genotype %>% select(Position, p.val.Genotype, 
                                       meth.group1.WT, meth.group2.TG, meth.diff.Genotype, estimate.Genotype, std.error.Genotype, pseudo.R.sqrt) %>% 
  left_join(., retainedClusterID$Genotype[,c("Position","cluster.id")], by = "Position") %>%
  tidyr::separate_wider_delim(., cols = Position, delim = ":", names = c("chr", "pos")) %>%
  `colnames<-`(betaColNames)


## ---------- null hypothesis -----------------

## Take sample beginning with letter "S" rows of the ColData predictedMeth: full set of samples of each time point and genotype
S_samples <- row.names(phenotype$rTg4510 %>% filter(grepl("S",rownames(.))))
M_samples <- row.names(phenotype$rTg4510 %>% filter(grepl("M",rownames(.))))
null_samples <- c(S_samples, M_samples)
predictedMethNull <- predictedMeth[, rownames(colData(predictedMeth)) %in% null_samples]
colData(predictedMethNull)$group.null <- rep(c(1,2), each = 8)

# To ensure that the P values are roughly uniformly distributed to get a variance of the Z scores that is Gaussian with variance 1,
# we recommend to estimate the variogram (and hence the correlation of Z scores) under the null hypothesis

# calculate betaResultsNull
betaResultsNull <- betaRegression(formula = ~group.null,
                                  link = "probit",
                                  object = predictedMethNull,
                                  type="BR",
                                  mc.cores = 16)


# Save all ojects from betaResultsNull
save(betaResultsNull, file = paste0(dirnames$processed, "/rrbs/biseqBetaResultsNull_rTg4510.RData"))
