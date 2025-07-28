#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## 23.05.2024: run DMP analysis using BiSeq
## model: - ~ Genotype + Age + Genotype*Age
## --------------------------------

## ---------- packages -----------------

suppressMessages(library(betareg))
suppressMessages(library(parallel))
suppressMessages(library(dplyr))


## ---------- input data -----------------

message("import files")
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "1_RRBS/functions/BetaRegressionDMPs.R"))

# smoothed BiSeq RRBS dataset
load(file = paste0(dirnames$processed, "/rrbs/rTg4510_RRBS_SmoothBetas.RData"))

# include only the samples that are in the smoothed beta regression matrix
message("Number of samples:", nrow(phenotype$rTg4510))

## ---------- Beta regression -----------------

# Beta regression ~Genotype + Age + Genotype*Age
message("Running beta regression: ~ Genotype + Age + Genotype*Age")
betaResults <- BetaRegressionRRBSinteraction(betas = RRBS_smoothbetas,
                                             pheno = phenotype$rTg4510,
                                             formula = ~Genotype + Age_months + Genotype*Age_months,
                                             formulaNull = ~Genotype + Age_months,
                                             link = "probit",
                                             num.cores = detectCores())

# Change columns from factor to numeric
betaResults[1:14] <- lapply(betaResults[1:14], as.character)
betaResults[1:14] <- lapply(betaResults[1:14], as.numeric)

# Check and remove NAs
## check how many NAs there are in the dataframe
na_df <- as.data.frame(betaResults[rowSums(is.na(betaResults)) > 0,]) 
nrow(na_df) # number of sites that have no data (NA)

# remove NAs
betaResultsClean <- betaResults[complete.cases(betaResults),]


## ---------- Output -----------------

# Add FDR correction
betaResultsClean <- betaResultsClean %>% mutate(FDR_adj_genotype = p.adjust(betaResultsClean[,"p.val.Genotype"], method = "fdr"),
                            FDR_adj_age = p.adjust(betaResultsClean[,"p.val.Age"], method = "fdr"),
                            FDR_adj_interaction = p.adjust(betaResultsClean[,"p.val.Interaction"], method = "fdr")) %>% 
  select(Position, FDR_adj_genotype, p.val.Genotype, meth.group1.WT, meth.group2.TG,
         meth.diff.Genotype, estimate.Genotype, std.error.Genotype, pseudo.R.sqrt,
         FDR_adj_age, p.val.Age, estimate.Age, std.error.Age,
         FDR_adj_interaction,p.val.Interaction, estimate.Interaction, std.error.Interaction,p.val.modelLRT)


save(betaResultsClean, file = paste0(dirnames$differential, "/rrbs/rTg4510_betaResultsDMP_GenotypeInteraction.RData"))

# Filter by significance (genotype and interaction effect)
rTg4510_rrbs_DMP <- list(
  Genotype = betaResultsClean %>% filter(FDR_adj_genotype < 0.05),
  Interaction = betaResultsClean %>% filter(FDR_adj_interaction < 0.05)
)
save(rTg4510_rrbs_DMP, file = paste0(dirnames$differential, "/rrbs/rTg4510_sigResultsDMP_GenotypeInteraction.RData"))
