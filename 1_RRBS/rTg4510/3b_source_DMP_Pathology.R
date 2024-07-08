#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## 13.06.2024: run DMP analysis using BiSeq
## model: - ~ Pathology
## --------------------------------

## ---------- packages -----------------

suppressMessages(library(betareg))
suppressMessages(library(parallel))
suppressMessages(library(dplyr))


## ---------- input data -----------------

message("import files")
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "1_RRBS/functions/BetaRegressionDMPs-Pathology.R"))

# smoothed BiSeq RRBS dataset
load(file = paste0(dirnames$processed, "/rrbs/rTg4510_RRBS_SmoothBetas.RData"))

# Remove samples that have NA for pathology 
message("Number of samples: ", nrow(phenotype$rTg4510))
phenotype$rTg4510 <- phenotype$rTg4510[!is.na(phenotype$rTg4510$ECX),]
message("Number of samples after removing samples without pathology data: ", nrow(phenotype$rTg4510))


## ---------- Beta regression -----------------

# common samples 
common_samples <- intersect(row.names(phenotype$rTg4510), colnames(RRBS_smoothbetas))
RRBS_smoothbetas <- RRBS_smoothbetas %>% select(all_of(common_samples))
phenotype$rTg4510 <- phenotype$rTg4510 %>% filter(row.names(.) %in% common_samples) 

# Beta regression ~ Pathology
message("Running beta regression: ~ Pathology")
betaResults <- BetaRegressionRRBSpathology(betas = RRBS_smoothbetas,
                                           pheno = phenotype$rTg4510,
                                           formula = ~ ECX,
                                           link = "probit",
                                           num.cores = detectCores())


# Change columns from factor to numeric
betaResults[1:3] <- lapply(betaResults[1:3], as.character)
betaResults[1:3] <- lapply(betaResults[1:3], as.numeric)

# Check and remove NAs
## check how many NAs there are in the dataframe
na_df <- as.data.frame(betaResults[rowSums(is.na(betaResults)) > 0,]) 
nrow(na_df) # number of sites that have no data (NA)

# remove NAs
betaResultsClean <- betaResults[complete.cases(betaResults),]


## ---------- Output -----------------

# Add FDR correction
betaResultsClean_Pathology <- betaResultsClean %>% 
  mutate(FDR_adj_pathology = p.adjust(betaResultsClean[,"p.val.Pathology"], method = "fdr")) %>% 
  select(Position, FDR_adj_pathology, p.val.Pathology, estimate.Pathology, std.error.Pathology)
save(betaResultsClean_Pathology, file = paste0(dirnames$differential, "/rrbs/rrTg4510_betaResultsDMP_Pathology.RData"))

# Filter by significance (genotype and interaction effect)
rTg4510_rrbs_DMP_Pathology <- betaResultsClean_Pathology %>% filter(FDR_adj_pathology < 0.05)
save(rTg4510_rrbs_DMP_Pathology, file = paste0(dirnames$differential, "/rrbs/rTg4510_sigResultsDMP_Pathology.RData"))

# merge with genotype results
# significant results
load(paste0(dirnames$differential, "/rrbs/rTg4510_sigResultsDMP_GenotypeInteraction.RData"))
rTg4510_rrbs_DMP$Pathology <- rTg4510_rrbs_DMP_Pathology
save(rTg4510_rrbs_DMP, file = paste0(dirnames$differential, "/rrbs/rTg4510_sigResultsDMPs.RData"))

load(paste0(dirnames$differential, "/rrbs/rTg4510_betaResultsDMP_GenotypeInteraction.RData"))
betaResultsCleanMerged <- list(
  Genotype = betaResultsClean,
  Pathology = betaResultsClean_Pathology
)
betaResultsClean <- betaResultsCleanMerged
save(betaResultsClean, file = paste0(dirnames$differential, "/rrbs/rTg4510_betaResultsDMPs.RData"))
