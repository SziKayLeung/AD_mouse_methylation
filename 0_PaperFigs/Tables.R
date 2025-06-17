# directory input
rootDir = "/lustre/projects/Research_Project-191406"
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "0_PaperFigs/paper_import.config.R"))

## ------------------- load ------------------

rTg4510_ECX_sigResultsDMPs <- get(load(paste0(dirnames$annotated,"/final/rTg4510_ECX_sigResultsDMPs.RData")))
rTg4510_tissueArray_sigResultsDMPs <- get(load(paste0(dirnames$annotated,"/final/rTg4510_tissueArray_sigResultsDMPs.RData")))

J20_ECX_sigResultsDMPs <- get(load(paste0(dirnames$annotated,"/final/J20_ECX_sigResultsDMPs.RData")))
J20_tissueArray_sigResultsDMPs <- get(load(paste0(dirnames$annotated,"/final/J20_tissueArray_sigResultsDMPs.RData")))

# Pathology tables, insert column if also identified in genotype analysis, identified in interaction analysis
rTg4510_ECX_sigResultsDMPs$Pathology <- rTg4510_ECX_sigResultsDMPs$Pathology %>% mutate(
  InteractionEffect = ifelse(Position %in% rTg4510_ECX_sigResultsDMPs$GenotypeAge$Position, TRUE, FALSE)
)
J20_ECX_sigResultsDMPs$Pathology <- J20_ECX_sigResultsDMPs$Pathology %>% mutate(
  InteractionEffect = ifelse(Position %in% J20_ECX_sigResultsDMPs$GenotypeAge$Position, TRUE, FALSE)
)


sigResArrayHIP$rTg4510$Genotype <- sigResArrayHIP$rTg4510$Genotype %>% 
  mutate(ECX = ifelse(Position %in% c(sigResArrayECX$rTg4510$Genotype$Position,intersect(sigResArrayECX$rTg4510$GenotypeAge$Position,sigResArrayECX$rTg4510$Pathology$Position)), TRUE, FALSE))

sigResArrayHIP$rTg4510$Pathology <- sigResArrayHIP$rTg4510$Pathology %>% 
  mutate(ECX = ifelse(Position %in% c(sigResArrayECX$rTg4510$Genotype$Position,intersect(sigResArrayECX$rTg4510$GenotypeAge$Position,sigResArrayECX$rTg4510$Pathology$Position)), TRUE, FALSE))

sigResArrayHIP$J20$Genotype <- sigResArrayHIP$J20$Genotype %>% 
  mutate(ECX = ifelse(Position %in% c(sigResArrayECX$J20$Genotype$Position,intersect(sigResArrayECX$J20$GenotypeAge$Position,sigResArrayECX$J20$Pathology$Position)), TRUE, FALSE))

sigResArrayHIP$J20$Pathology <- sigResArrayHIP$J20$Pathology %>% 
  mutate(ECX = ifelse(Position %in% c(sigResArrayECX$J20$Genotype$Position,intersect(sigResArrayECX$J20$GenotypeAge$Position,sigResArrayECX$J20$Pathology$Position)), TRUE, FALSE)) 


# unique DMPs in ECX array but not in hippocampus array
ECXrTg4510Unique2HIP <- setdiff(
  # ECX
  c(rTg4510_array_sig$ECX$Genotype$position,
    rTg4510_array_sig$ECX$Interaction$position,
   rTg4510_array_sig$ECX$Pathology$position),
  # HIP
  c(rTg4510_array_sig$HIP$Genotype$position,
    rTg4510_array_sig$HIP$Interaction$position,
    rTg4510_array_sig$HIP$Pathology$position))

ECXrTg4510Unique2HIPTable <- rbind(
  rTg4510_array_sig$ECX$Genotype[rTg4510_array_sig$ECX$Genotype$position %in% ECXrTg4510Unique2HIP,] %>% 
    arrange(FDR_adj_Genotype) %>% 
    dplyr::select(position, cpg, FDR_adj_Genotype, ChIPseeker_GeneSymbol, distanceToTSS) %>% 
    rename("FDR_adj_Genotype" = "FDR_adj") %>% mutate(test = "genotype"),
  rTg4510_array_sig$ECX$Interaction[rTg4510_array_sig$ECX$Interactionposition %in% ECXrTg4510Unique2HIP,] %>% 
    arrange(FDR_adj_GenotypeAge) %>% 
    dplyr::select(position, cpg, FDR_adj_GenotypeAge, ChIPseeker_GeneSymbol, distanceToTSS) %>% 
    rename("FDR_adj_GenotypeAge" = "FDR_adj") %>% mutate(test = "interaction"),
  rTg4510_array_sig$ECX$Pathology[rTg4510_array_sig$ECX$Pathology$position %in% ECXrTg4510Unique2HIP,] %>% 
    arrange(FDR_adj_Pathology) %>% 
    dplyr::select(position, cpg, FDR_adj_Pathology, ChIPseeker_GeneSymbol, distanceToTSS) %>%
    rename("FDR_adj_Pathology" = "FDR_adj") %>% mutate(test = "pathology")
)

plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                  ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr11:97685376", pathology = TRUE)
plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                  ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr5:52157628", pathology = TRUE)


# combine pathology + genotype
combine_sig_genotype_pathology <- function(sigGenotype, sigPathology){
  
  commonAnnotationCols <- c("Position","meanWT", "meanTG", "delta", "directionTG", "ChIPseeker_GeneEnsembl", "ChIPseeker_GeneSymbol", "ChIPseeker_GeneName","ChIPseeker_Annotation","ChIPseeker_TransEnsembl","distanceToTSS")
  
  sigPathology <- sigPathology %>% mutate(PathologyEffect = TRUE)
  dat <- merge(sigGenotype,
               sigPathology,
               by = "Position", all = T) %>%
    mutate(cpg = coalesce(cpg.x, cpg.y),
           Platform = coalesce(Platform.x, Platform.y),
           meanWT = coalesce(meanWT.x, meanWT.y),
           meanTG = coalesce(meanTG.x, meanTG.y),
           delta  = coalesce(delta.x, delta.y),
           directionTG = coalesce(directionTG.x, directionTG.y),
           ChIPseeker_GeneSymbol = coalesce(ChIPseeker_GeneSymbol.x, ChIPseeker_GeneSymbol.y),
           ChIPseeker_GeneName = coalesce(ChIPseeker_GeneName.x, ChIPseeker_GeneName.y),
           ChIPseeker_GeneEnsembl = coalesce(ChIPseeker_GeneEnsembl.x, ChIPseeker_GeneEnsembl.y),
           ChIPseeker_Annotation = coalesce(ChIPseeker_Annotation.x, ChIPseeker_Annotation.y),
           ChIPseeker_TransEnsembl = coalesce(ChIPseeker_TransEnsembl.x, ChIPseeker_TransEnsembl.y),
           distanceToTSS = coalesce(distanceToTSS.x, distanceToTSS.y)) %>% 
    dplyr::select(-contains(".x"),-contains(".y")) %>%
    mutate(across(contains("Effect"), ~ tidyr::replace_na(., FALSE))) %>%
    dplyr::select(Position, Platform, cpg, "GenotypeEffect","PathologyEffect","InteractionEffect", contains("_Genotype"), contains("_Pathology"), all_of(commonAnnotationCols))
  
  return(dat)
}
combinedSigResults <- list(
  rTg4510 = combine_sig_genotype_pathology(rTg4510_ECX_sigResultsDMPs$Genotype, rTg4510_ECX_sigResultsDMPs$Pathology),
  J20 = combine_sig_genotype_pathology(J20_ECX_sigResultsDMPs$Genotype, J20_ECX_sigResultsDMPs$Pathology)
)

## ------------------- table output ------------------

write.csv(rTg4510_ECX_sigResultsDMPs$Genotype, paste0(dirnames$paper,"/tables/rTg4510_ECX_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(rTg4510_ECX_sigResultsDMPs$GenotypeAge, paste0(dirnames$paper,"/tables/rTg4510_ECX_sigResultsDMPs_GenotypeAge.csv"), quote = T)
write.csv(rTg4510_ECX_sigResultsDMPs$Pathology, paste0(dirnames$paper,"/tables/rTg4510_ECX_sigResultsDMPs_Pathology.csv"), quote = T)
write.csv(rTg4510_tissueArray_sigResultsDMPs, paste0(dirnames$paper,"/tables/rTg4510_tissueArray_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResArrayHIP$rTg4510$Genotype, paste0(dirnames$paper,"/tables/rTg4510_HIP_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResArrayHIP$rTg4510$Pathology, paste0(dirnames$paper,"/tables/rTg4510_HIP_sigResultsDMPs_Pathology.csv"), quote = T)
write.csv(combinedSigResults$rTg4510, paste0(dirnames$paper,"/tables/rTg4510_ECX_sigResultsDMPs_Combined.csv"), quote = T)

# J20
write.csv(J20_ECX_sigResultsDMPs$Genotype, paste0(dirnames$paper,"/tables/J20_ECX_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(J20_ECX_sigResultsDMPs$GenotypeAge, paste0(dirnames$paper,"/tables/J20_ECX_sigResultsDMPs_GenotypeAge.csv"), quote = T)
write.csv(J20_ECX_sigResultsDMPs$Pathology, paste0(dirnames$paper,"/tables/J20_ECX_sigResultsDMPs_Pathology.csv"), quote = T)
write.csv(J20_tissueArray_sigResultsDMPs, paste0(dirnames$paper,"/tables/J20_tissueArray_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResArrayHIP$J20$Genotype, paste0(dirnames$paper,"/tables/J20_HIP_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResArrayHIP$J20$Pathology, paste0(dirnames$paper,"/tables/J20_HIP_sigResultsDMPs_Pathology.csv"), quote = T)
write.csv(combinedSigResults$J20, paste0(dirnames$paper,"/tables/J20_ECX_sigResultsDMPs_Combined.csv"), quote = T)

# GO analysis
extract_postions_as_bed(sigRes$rTg4510$Genotype$Position, paste0(dirnames$GO,"rTg4510_sig_ECX_Genotype_coordinates.bed"))
extract_postions_as_bed(sigRes$J20$Genotype$Position, paste0(dirnames$GO,"J20_sig_ECX_Genotype_coordinates.bed"))
extract_postions_as_bed(unique(sigRes$J20$PathologyCommonInteraction$Position), paste0(dirnames$GO,"J20_sig_ECX_PathologyCommonInteraction_coordinates.bed"))
extract_postions_as_bed(unique(sigRes$rTg4510$PathologyCommonInteraction$Position), 
                        paste0(dirnames$GO,"rTg4510_sig_ECX_PathologyCommonInteraction_coordinates.bed"))
extract_postions_as_bed(unique(c(sigRes$rTg4510$Genotype$Position,sigRes$rTg4510$PathologyCommonInteraction$Position)),
                        paste0(dirnames$GO,"rTg4510_sig_ECX_GenotypePathologyCommonInteraction_coordinates.bed"))
# note over 1M sites and GREAT cannot handle 
extract_postions_as_bed((unique(c(row.names(rTg4510_rrbs_beta), row.names(rTg4510_array_beta)))), paste0(dirnames$GO,"rTg4510_RRBSArray_ECX_coordinates.bed"))
# subset to promoter sites in RRBS only and all Array sites
# background needs to include the significant sites tested
GO_background_sites <- c(unique(c(anno_rrbs_all$rTg4510[anno_rrbs_all$rTg4510$Promoter == "TRUE","position"], 
                                 row.names(rTg4510_array_beta),
                                 sigRes$rTg4510$Genotype$Position,
                                 sigRes$rTg4510$PathologyCommonInteraction$Position,
                                 sigRes$J20$Genotype$Position,
                                 sigRes$J20$PathologyCommonInteraction$Position
)))
setdiff(sigRes$J20$PathologyCommonInteraction$Position, GO_background_sites)

extract_postions_as_bed(GO_background_sites, paste0(dirnames$GO,"rTg4510_PromoterRRBSArray_ECX_coordinates.bed"))



