#!/usr/bin/env Rscript
## ----------Script-----------------  
##
## Purpose: Annotate differentially-methylated positions using chipseeker in array data
##         
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------

##--------------- packages 

suppressMessages(library(stringr))
suppressMessages(library(dplyr))


# directory input
rootDir = "/lustre/projects/Research_Project-191406"
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/"

source(paste0(scriptDir, "import.config"))
# source functions to run script
source(paste0(scriptDir, "1_RRBS/functions/chipSeekerAnnotation.R"))

##-------------- input -------

# read J20 significant results
sigResults <- list.files(path = paste0(dirnames$differential,"/array"), pattern = "^J20.*sigResultsDMPs.RData", full.names = T)
for(i in 1:length(sigResults)){load(sigResults[i])}

##-------------- annotate peaks using ChipSeeker, rename column names

anno_chipseeker <- list(
  ECX = lapply(J20_array_ECX_DMP, function(x) ChipRenameCol(ChipAnnotatePeaks(x, merge = "yes"))),
  HIP = lapply(J20_array_HIP_DMP, function(x) ChipRenameCol(ChipAnnotatePeaks(x, merge = "yes"))),
  tissue = lapply(J20_array_tissue_DMP, function(x) ChipRenameCol(ChipAnnotatePeaks(x, merge = "yes")))
)


##-------------- Write output to xlsx file 

for(i in 1:3){
  write.csv(anno_chipseeker$ECX[[i]], paste0(dirnames$annotated, "/array/J20_array_ECX_", names(anno_chipseeker$ECX)[[i]],".csv"))
  write.csv(anno_chipseeker$HIP[[i]], paste0(dirnames$annotated, "/array/J20_array_HIP_", names(anno_chipseeker$HIP)[[i]],".csv"))
}
write.csv(anno_chipseeker$tissue[[1]], paste0(dirnames$annotated, "/array/J20_array_tissue_", names(anno_chipseeker$tissue)[[1]],".csv"))


##-------------- save output to Rdata (all files)

J20_array_anno <- anno_chipseeker
save(J20_array_anno, file = paste0(dirnames$annotated, "/array/J20_array_annoSigResultsDMPs.RData"))
