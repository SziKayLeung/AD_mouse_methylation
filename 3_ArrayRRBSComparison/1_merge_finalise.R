#!/usr/bin/env Rscript
## ----------Script-----------------  
##
## Purpose: Process and merge RRBS and array data from rTg4510 and J20
##          Common probes between the array and RRBS
## 
## Date:    03 July 2024        
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------

# directory input
rootDir = "/lustre/projects/Research_Project-191406"
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/"
source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "3_ArrayRRBSComparison/functions/summaryStatsDMP.R"))

#-------------- input -------------

## smoothed BiSeq RRBS dataset
rTg4510_rrbs_beta <- get(load(file = paste0(dirnames$processed, "/rrbs/rTg4510_RRBS_SmoothBetas.RData")))
rTg4510_array_beta <- get(load(file = paste0(dirnames$processed, "/array/Normalised_rTg4510_array_ECX.RData")))
rTg4510_array_HIP_beta <- get(load(file = paste0(dirnames$processed, "/array/Normalised_rTg4510_array_HIP.RData")))
J20_rrbs_beta <- get(load(file = paste0(dirnames$processed, "/rrbs/J20_RRBS_SmoothBetas.RData")))
J20_array_beta <- get(load(file = paste0(dirnames$processed, "/array/Normalised_J20_array_ECX.RData")))
J20_array_HIP_beta <- get(load(file = paste0(dirnames$processed, "/array/Normalised_J20_array_HIP.RData")))

# remove probes with no chromosomal position (due to patch)
# duplicated probes per position in array; take average

datawrangle_beta <- function(inputBeta){
  inputBeta <- inputBeta[!is.na(inputBeta$position),] %>% dplyr::select(-Gene_Symbol, -annotation, -cpg) 
  inputBeta <- as.data.frame(inputBeta %>% group_by(position) %>% summarise_all("mean")) %>% tibble::column_to_rownames(., var = "position")
  return(inputBeta)
}

rTg4510_array_beta <- datawrangle_beta(rTg4510_array_beta)
J20_array_beta <- datawrangle_beta(J20_array_beta)
rTg4510_array_HIP_beta <- datawrangle_beta(rTg4510_array_HIP_beta)
J20_array_HIP_beta <- datawrangle_beta(J20_array_HIP_beta)

# save datawrangled results
save(rTg4510_array_beta, file = paste0(dirnames$processed, "/array/Normalised_rTg4510_array_ECX_Final.RData"))
save(J20_array_beta, file = paste0(dirnames$processed, "/array/Normalised_J20_array_ECX_Final.RData"))
save(rTg4510_array_HIP_beta, file = paste0(dirnames$processed, "/array/Normalised_rTg4510_array_HIP_Final.RData"))
save(J20_array_HIP_beta, file = paste0(dirnames$processed, "/array/Normalised_J20_array_HIP_Final.RData"))

## differential all results
rTg4510_rrbs_results <- get(load( file = paste0(dirnames$differential, "/rrbs/rTg4510_betaResultsDMPs.RData")))
J20_rrbs_results <- get(load( file = paste0(dirnames$differential, "/rrbs/J20_betaResultsDMPs.RData")))

## annotated significant results
rTg4510_rrbs_sig <- get(load(file = paste0(dirnames$annotated, "/rrbs/rTg4510_rrbs_annoSigResultsDMPs.RData")))
rTg4510_array_sig <- get(load(file = paste0(dirnames$annotated, "/array/rTg4510_array_annoSigResultsDMPs.RData")))

J20_rrbs_sig <- get(load(file = paste0(dirnames$annotated, "/rrbs/J20_rrbs_annoSigResultsDMPs.RData")))
J20_array_sig <- get(load(file = paste0(dirnames$annotated, "/array/J20_array_annoSigResultsDMPs.RData")))

# merge RRBS and Array results
rTg4510_merge_sig <- MergeArrayRRBSSigResults(rrbsBeta=rTg4510_rrbs_beta, 
                                              arrayBeta=rTg4510_array_beta, 
                                              rrbsSigResults=rTg4510_rrbs_sig, 
                                              arraySigResults=rTg4510_array_sig, 
                                              phenotypeInput=phenotype$rTg4510)
J20_merge_sig <- MergeArrayRRBSSigResults(J20_rrbs_beta, J20_array_beta, J20_rrbs_sig, J20_array_sig, phenotype$J20)

# save signfiicant results
rTg4510_ECX_sigResultsFull <- rTg4510_merge_sig$sig_ECX_Full
J20_ECX_sigResultsFull <- J20_merge_sig$sig_ECX_Full
save(rTg4510_ECX_sigResultsFull, file = paste0(dirnames$annotated,"/final/rTg4510_ECX_sigResultsDMPs.RData"))
save(J20_ECX_sigResultsFull, file = paste0(dirnames$annotated,"/final/J20_ECX_sigResultsDMPs.RData"))

# save significant beta 
rTg4510_sig_betaFull <- rTg4510_merge_sig$sig_betaECX_Full
J20_sig_betaFull  <- J20_merge_sig$sig_betaECX_Full
save(rTg4510_sig_betaFull, file = paste0(dirnames$annotated,"/final/rTg4510_ECX_sigBeta.RData"))
save(J20_sig_betaFull, file = paste0(dirnames$annotated,"/final/J20_ECX_sigBeta.RData"))

# array tissue
rTg4510_tissue_sigResults <- datawrangle_tissue_array(rTg4510_array_sig$tissue$Genotype)
J20_tissue_sigResults <- datawrangle_tissue_array(J20_array_sig$tissue$Genotype)
save(rTg4510_tissue_sigResults, file = paste0(dirnames$annotated,"/final/rTg4510_tissueArray_sigResultsDMPs.RData"))
save(J20_tissue_sigResults, file = paste0(dirnames$annotated,"/final/J20_tissueArray_sigResultsDMPs.RData"))

## entorhinal cortex
sigResArrayECX <- list(
  rTg4510 = sigArrayResults(rTg4510_array_beta, rTg4510_array_sig$ECX, phenotype$rTg4510),
  J20 = sigArrayResults(J20_array_beta, J20_array_sig$ECX, phenotype$J20)
)
save(sigResArrayECX, file = paste0(dirnames$annotated,"/final/entorhinalcortexArray_sigResultsDMPs.RData"))

## hippocampus 
sigResArrayHIP  <- list(
  rTg4510 = sigArrayResults(rTg4510_array_HIP_beta, rTg4510_array_sig$HIP, phenotype$rTg4510_HIP),
  J20 = sigArrayResults(J20_array_HIP_beta, J20_array_sig$HIP, phenotype$J20_HIP)
)
save(sigResArrayHIP, file = paste0(dirnames$annotated,"/final/hippocampusArray_sigResultsDMPs.RData"))
