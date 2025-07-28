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

load(paste0(dirnames$processed, "/rrbs/biseqBetaResultsNull_rTg4510.RData"))

## ---------- calculate vario -----------------

vario <- makeVariogram(betaResultsNull)
vario.sm <- smoothVariogram(vario, sill = 0.9) 

## ---------- output -----------------

save(vario, vario.sm, file = paste0(dirnames$processed, "/rrbs/varioNull_rTg4510.RData"))
