#!/usr/bin/env Rscript
## ----------Script-----------------  
##
## Purpose: Re-annotate differentially-methylated positions using chipseeker and alignment of gtf
##          Barplot of chipseeker annotations in rTg4510 and J20 mouse model (location)
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------
## 
## Motivation
## 1. E.Walker used chipseeker for annotation of DMP
## 2. I.Castanho noticed misannotation of top DMP (in exon) to Mir6388 rather than Gm30948.
## 3. S.Leung to validate and correct annotations (see ChipSeeker_wrongannotations.html)
## 4. S.Leung noticed another separate method for annotation, by annotating the EnsemblID in the the “annotation” column from ChipSeeker output with mouse reference gtf
##    Such approach, while limited to DMPs located in exon and intron, further annotated DMP to Gm30948     
## 5. After discussion with J.Mill, E.Hannon..., to include column for Chipseeker annotations and annotation from merging with reference gtf
## 6. Use Chipseeker column to manually output annotation locations in order to include results across all and significant sites, and in both mouse models
## 
## Function
## ChipAnnotatePeaks - original function called for chipseeker annotation used by E.Walker
## 

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

# read rTg4510 significant results
rTg4510_DMP <- get(load(paste0(dirnames$differential, "/rrbs/rTg4510_sigResultsDMPs.RData")))

# beta values
rTg4510_rrbs_beta <- get(load(file = paste0(dirnames$processed, "/rrbs/rTg4510_RRBS_SmoothBetas.RData")))
J20_rrbs_beta <- get(load(file = paste0(dirnames$processed, "/rrbs/J20_RRBS_SmoothBetas.RData")))

# mouse reference
mouseRef_gtf <- as.data.frame(rtracklayer::import(paste0(dirnames$utils, "gencode.vM22.annotation.gtf"))) 


##-------------- annotate peaks using ChipSeeker

anno_chipseeker <- lapply(rTg4510_DMP, function(x) ChipAnnotatePeaks(x))

# annotate all RRBS peaks
anno_rrbs_all <- list( 
  rTg4510 = ChipAnnotatePeaks(rTg4510_rrbs_beta %>% tibble::rownames_to_column(., var = "Position")),
  J20 = ChipAnnotatePeaks(J20_rrbs_beta %>% tibble::rownames_to_column(., var = "Position"))
)

##-------------- classify and plot annotations

# classify peaks into simplified categories 
class_anno_chipseeker <- lapply(anno_chipseeker, function(x) ChipClassifyAnno_df(x))


##-------------- further annotate peaks in intron and exon (defined from ChipSeeker) using reference gtf

anno_chipseeker_mergedref <- lapply(anno_chipseeker, function(x) ChipMergeRefGtf(x, mouseRef_gtf))


##-------------- finalise annotations from ChipSeeker and merging with reference gtf

anno_final <- list()
for(i in 1:3){anno_final[[i]] <- ChipFinaliseAnno(anno_chipseeker[[i]], anno_chipseeker_mergedref[[i]], rTg4510_DMP[[i]])}
names(anno_final) <- names(anno_chipseeker)[1:3]

                   
##-------------- Write output to xlsx file 

for(i in 1:3){
  write.csv(anno_final[[i]], paste0(dirnames$annotated, "/rrbs/rTg4510_rrbs_ECX_", names(anno_final)[[i]],".csv"))
}

##-------------- save output to Rdata (all files)

rTg4510_rrbs_anno <- anno_final
save(rTg4510_rrbs_anno, file = paste0(dirnames$annotated, "/rrbs/rTg4510_rrbs_annoSigResultsDMPs.RData"))
save(anno_rrbs_all, file = paste0(dirnames$annotated, "/rrbs/rrbs_annoAllPositions.RData"))
