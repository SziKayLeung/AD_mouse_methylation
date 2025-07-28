#!/usr/bin/env Rscript
## ----------Script-----------------  
##
## Purpose: Annotate differentially-methylated positions using chipseeker in array data
##         
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------

library("stringr")
library("cowplot")
library("ggplot2")

# directory input
rootDir = "/lustre/projects/Research_Project-191406"
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/"
source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "3_ArrayRRBSComparison/functions/summaryStatsDMP.R"))

LOGEN_ROOT="/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen"
source(paste0(LOGEN_ROOT, "/aesthetics_basics_plots/draw_density.R"))
source(paste0(LOGEN_ROOT, "/aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "/aesthetics_basics_plots/draw_venn.R"))


## smoothed BiSeq RRBS dataset
rTg4510_rrbs_beta <- get(load(file = paste0(dirnames$processed, "/rrbs/rTg4510_RRBS_SmoothBetas.RData"))); rm(RRBS_smoothbetas)
rTg4510_array_beta <- get(load(file = paste0(dirnames$processed, "/array/Normalised_rTg4510_array_ECX_Final.RData")))
J20_rrbs_beta <- get(load(file = paste0(dirnames$processed, "/rrbs/J20_RRBS_SmoothBetas.RData"))); rm(RRBS_smoothbetas)
J20_array_beta <- get(load(file = paste0(dirnames$processed, "/array/Normalised_J20_array_ECX_Final.RData")))

## annotated significant results
rTg4510_rrbs_sig <- get(load(file = paste0(dirnames$annotated, "/rrbs/rTg4510_rrbs_annoSigResultsDMPs.RData"))); rm(rTg4510_rrbs_anno)
rTg4510_array_sig <- get(load(file = paste0(dirnames$annotated, "/array/rTg4510_array_annoSigResultsDMPs.RData"))); rm(rTg4510_array_anno)
J20_rrbs_sig <- get(load(file = paste0(dirnames$annotated, "/rrbs/J20_rrbs_annoSigResultsDMPs.RData"))); rm(J20_rrbs_anno)
J20_array_sig <- get(load(file = paste0(dirnames$annotated, "/array/J20_array_annoSigResultsDMPs.RData"))); rm(J20_array_anno)

# common probes between rrbs and array
rTg4510_common_probes <- commonStatsDescription(rTg4510_rrbs_beta, rTg4510_array_beta)
J20_common_probes <- commonStatsDescription(J20_rrbs_beta, J20_array_beta)


#-------------- correlation of all common probes detected -------------

# rTg4510
corrPlotCommonProbes(rTg4510_rrbs_beta, rTg4510_array_beta, rTg4510_common_probes)
# J20
corrPlotCommonProbes(J20_rrbs_beta, J20_array_beta, J20_common_probes)

#-------------- DMP genotype: direction of effect -------------

plot_grid(twovenndiagrams(rTg4510_rrbs_sig$Genotype$Position,rTg4510_rrbs_sig$Pathology$Position,"Genotype","Pathology"))
message("Number of overlapping probes between array and rrbs (sig): ", 
        length(intersect(rTg4510_rrbs_sig$Genotype$Position,rTg4510_array_sig$ECX$Genotype$position)))

## --------------------- rrbs significant, validated by array 

plot_platform_commonProbes(rTg4510_rrbs_beta, rTg4510_array_beta, rTg4510_rrbs_sig$Genotype, 
                           phenotype$rTg4510, rTg4510_common_probes)


# individual probes plot
plot_grid(
  plot_DMP(betaMatrix=rTg4510_rrbs_beta, phenotypeFile=phenotype$rTg4510, position=intersect(rTg4510_rrbs_sig$Genotype$Position, rTg4510_common_probes)),
  plot_DMP(betaMatrix=rTg4510_array_beta, phenotypeFile=phenotype$rTg4510, position=intersect(rTg4510_rrbs_sig$Genotype$Position, rTg4510_common_probes)),
  ncol = 1
)


## --------------------- array significant, validated by RRBS

plot_platform_commonProbes(rTg4510_rrbs_beta, rTg4510_array_beta, rTg4510_array_sig$ECX$Genotype, 
                           phenotype$rTg4510, rTg4510_common_probes)

plot_grid(
  plot_DMP(betaMatrix=rTg4510_rrbs_beta, phenotypeFile=phenotype$rTg4510, position=intersect(rTg4510_array_sig$ECX$Genotype$Position, rTg4510_common_probes)) + ylim(0,1),
  plot_DMP(betaMatrix=rTg4510_array_beta, phenotypeFile=phenotype$rTg4510, position=intersect(rTg4510_array_sig$ECX$Genotype$Position, rTg4510_common_probes)) + ylim(0,1),
  ncol = 1
)

## -------------------- array tissue difference ---------------

cols2Keep <- c("ChIPseeker_GeneEnsembl","ChIPseeker_GeneSymbol", "ChIPseeker_GeneName","ChIPseeker_Annotation","ChIPseeker_TransEnsembl","distanceToTSS")
sig_tissue_Genotype <- rTg4510_array_sig$tissue$Genotype %>%
  mutate(position = Position, platform = "Array") %>%
  select(platform, position, cpg, "PrZ.Tissue", "PrZ.Tissue", "PrZ.Tissue.GenotypeTG", "FDR_adj_TissueGenotype", all_of(cols2Keep)) %>%
  setNames(c("Platform", "Position", "cpg", "P.value_Tissue", "P.value_TissueGenotype", "FDR_adj_TissueGenotype", cols2Keep))  %>%
  arrange(FDR_adj_TissueGenotype)
save(sig_tissue_Genotype, file = paste0(dirnames$annotated,"/final/rTg4510_tissueArray_sigResultsDMPs.RData"))
