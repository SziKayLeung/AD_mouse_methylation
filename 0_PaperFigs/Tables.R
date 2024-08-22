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

## ------------------- table output ------------------

write.csv(rTg4510_ECX_sigResultsDMPs$Genotype, paste0(dirnames$paper,"/tables/rTg4510_ECX_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(rTg4510_ECX_sigResultsDMPs$GenotypeAge, paste0(dirnames$paper,"/tables/rTg4510_ECX_sigResultsDMPs_GenotypeAge.csv"), quote = T)
write.csv(rTg4510_ECX_sigResultsDMPs$Pathology, paste0(dirnames$paper,"/tables/rTg4510_ECX_sigResultsDMPs_Pathology.csv"), quote = T)
write.csv(rTg4510_tissueArray_sigResultsDMPs, paste0(dirnames$paper,"/tables/rTg4510_tissueArray_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResHip$rTg4510$Genotype, paste0(dirnames$paper,"/tables/rTg4510_HIP_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResHip$rTg4510$Pathology, paste0(dirnames$paper,"/tables/rTg4510_HIP_sigResultsDMPs_Pathology.csv"), quote = T)


# J20
write.csv(J20_ECX_sigResultsDMPs$Genotype, paste0(dirnames$paper,"/tables/J20_ECX_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(J20_ECX_sigResultsDMPs$GenotypeAge, paste0(dirnames$paper,"/tables/J20_ECX_sigResultsDMPs_GenotypeAge.csv"), quote = T)
write.csv(J20_ECX_sigResultsDMPs$Pathology, paste0(dirnames$paper,"/tables/J20_ECX_sigResultsDMPs_Pathology.csv"), quote = T)
write.csv(J20_tissueArray_sigResultsDMPs, paste0(dirnames$paper,"/tables/J20_tissueArray_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResHip$J20$Genotype, paste0(dirnames$paper,"/tables/J20_HIP_sigResultsDMPs_Genotype.csv"), quote = T)
write.csv(sigResHip$J20$Pathology, paste0(dirnames$paper,"/tables/J20_HIP_sigResultsDMPs_Pathology.csv"), quote = T)
