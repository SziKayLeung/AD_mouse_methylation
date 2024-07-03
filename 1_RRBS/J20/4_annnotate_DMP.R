#!/usr/bin/env Rscript
## ----------Script-----------------  
##
## Purpose: annotate DMP in J20 ECX from RRBS dataset
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
J20_rrbs_DMP <- get(load(paste0(dirnames$differential, "/rrbs/J20_sigResultsDMPs.RData")))

mouseRef_gtf <- as.data.frame(rtracklayer::import(paste0(dirnames$utils, "gencode.vM22.annotation.gtf"))) 


##-------------- annotate peaks using ChipSeeker

anno_chipseeker <- lapply(J20_rrbs_DMP, function(x) ChipAnnotatePeaks(x))


##-------------- classify and plot annotations

# classify peaks into simplified categories 
class_anno_chipseeker <- lapply(anno_chipseeker, function(x) ChipClassifyAnno_df(x))


##-------------- further annotate peaks in intron and exon (defined from ChipSeeker) using reference gtf

anno_chipseeker_mergedref <- lapply(anno_chipseeker, function(x) ChipMergeRefGtf(x, mouseRef_gtf))


##-------------- finalise annotations from ChipSeeker and merging with reference gtf

anno_final <- list()
for(i in 1:3){anno_final[[i]] <- ChipFinaliseAnno(anno_chipseeker[[i]], anno_chipseeker_mergedref[[i]], J20_rrbs_DMP[[i]])}
names(anno_final) <- names(anno_chipseeker)[1:3]


##-------------- Write output to xlsx file 

for(i in 1:3){
  write.csv(anno_final[[i]], paste0(dirnames$annotated, "/rrbs/J20_rrbs_ECX_", names(anno_final)[[i]],".csv"))
}

##-------------- save output to Rdata (all files)

J20_rrbs_anno <- anno_final
save(J20_rrbs_anno, file = paste0(dirnames$annotated, "/rrbs/J20_rrbs_annoSigResultsDMPs.RData"))
