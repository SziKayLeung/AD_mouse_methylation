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

## ------------------- table output ------------------

write.csv(rTg4510_ECX_sigResultsDMPs$Genotype, paste0(dirnames$paper,"/tables/rTg4510_ECX_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(rTg4510_ECX_sigResultsDMPs$GenotypeAge, paste0(dirnames$paper,"/tables/rTg4510_ECX_sigResultsDMPs_GenotypeAge.csv"), quote = T)
write.csv(rTg4510_ECX_sigResultsDMPs$Pathology, paste0(dirnames$paper,"/tables/rTg4510_ECX_sigResultsDMPs_Pathology.csv"), quote = T)
write.csv(rTg4510_tissueArray_sigResultsDMPs, paste0(dirnames$paper,"/tables/rTg4510_tissueArray_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResArrayHIP$rTg4510$Genotype, paste0(dirnames$paper,"/tables/rTg4510_HIP_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResArrayHIP$rTg4510$Pathology, paste0(dirnames$paper,"/tables/rTg4510_HIP_sigResultsDMPs_Pathology.csv"), quote = T)


# J20
write.csv(J20_ECX_sigResultsDMPs$Genotype, paste0(dirnames$paper,"/tables/J20_ECX_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(J20_ECX_sigResultsDMPs$GenotypeAge, paste0(dirnames$paper,"/tables/J20_ECX_sigResultsDMPs_GenotypeAge.csv"), quote = T)
write.csv(J20_ECX_sigResultsDMPs$Pathology, paste0(dirnames$paper,"/tables/J20_ECX_sigResultsDMPs_Pathology.csv"), quote = T)
write.csv(J20_tissueArray_sigResultsDMPs, paste0(dirnames$paper,"/tables/J20_tissueArray_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResArrayHIP$J20$Genotype, paste0(dirnames$paper,"/tables/J20_HIP_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResArrayHIP$J20$Pathology, paste0(dirnames$paper,"/tables/J20_HIP_sigResultsDMPs_Pathology.csv"), quote = T)


# GO analysis
extract_postions_as_bed(sigRes$rTg4510$Genotype$Position, paste0(dirnames$GO,"rTg4510_sig_ECX_Genotype_coordinates.bed"))
extract_postions_as_bed((unique(c(row.names(rTg4510_rrbs_beta), row.names(rTg4510_array_beta)))), paste0(dirnames$GO,"rTg4510_RRBSArray_ECX_coordinates.bed"))