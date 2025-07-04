#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: convert human AD genes to mouse AD genes  
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## Date: 13/11/2023
##
## ---------- Notes -----------------
## 
## Using E.Walker gene conversion file 
## human and mouse homologues
## human_Thal file with single BCAR1 gene cannot be read, therefore 
## concatenated to the end of the Shireby2022_human_braak.txt fie
## and manually created a converted_Shireby2022_human_Thal.txt with the converted mouse gene
## 
## with exception of Shireby2022_human_metaS8.txt, all input files provided by I.Castanho from S2
## human_meta_S8.txt takes all annotated genes in column R "UCSC Nearest Gene" including those with delimited ";"

## ------------ input ------------

# directory input
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config"))

## from import.config
# read in gene conversion file --> mouse2human
# read in gene lists from G.Shireby and R.Smith lists (EWAS) --> inputHumanLists


## ------------ merge ------------

# using the mouse2human file with the same "human" column, keep only the mouse homologues genes commonly detected
convertedLists <- lapply(inputHumanLists, function(x) x %>% left_join(., mouse2human, by = "human") %>% dplyr::select(mouse))
# remove NA
convertedLists <- lapply(convertedLists, function(x) x[!is.na(x)])
# remove duplicates
convertedLists <- lapply(convertedLists, function(x) unique(x))

names(convertedLists) <- c("Shireby2022Braak", "Shireby2022Thal", "Shireby2022Meta", "Smith2021")

## ------------ output ------------

# write output to input directory
for(n in 1:length(convertedLists)){
  write.table(convertedLists[[n]], paste0(dirnames$humanAnnot, names(convertedLists)[[n]], "_geneList.txt"),quote=F,col.names=F,row.names = F)
}
