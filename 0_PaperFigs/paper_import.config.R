rootDir = "/lustre/projects/Research_Project-191406"
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "0_PaperFigs/Functions.R"))

sigBeta <- list(
  rTg4510 = get(load(paste0(dirnames$annotated,"/final/rTg4510_ECX_sigBeta.RData"))),
  J20 = get(load(paste0(dirnames$annotated,"/final/J20_ECX_sigBeta.RData")))
);  rm(rTg4510_sig_betaFull); rm(J20_sig_betaFull)

sigRes <- list(
  rTg4510 = get(load(file = paste0(dirnames$annotated,"/final/rTg4510_ECX_sigResultsDMPs.RData"))),
  J20 = get(load(paste0(dirnames$annotated,"/final/J20_ECX_sigResultsDMPs.RData")))
);  rm(rTg4510_ECX_sigResultsFull); rm(J20_ECX_sigResultsFull)

transgeneTg4510 <- c("Wdr60","Esyt2","Ncapg2","Ptprn2","Prnp","Prn","Mapt")
transgeneJ20 <- c("Zbt20")

sigRes_transgeneRem <- list(
  rTg4510 = lapply(sigRes$rTg4510, function(x) x %>% filter(!ChIPseeker_GeneSymbol %in% transgeneTg4510)),
  J20 = lapply(sigRes$J20, function(x) x %>% filter(!ChIPseeker_GeneSymbol %in% transgeneJ20))
)


## smoothed BiSeq RRBS dataset
rTg4510_rrbs_beta <- get(load(file = paste0(dirnames$processed, "/rrbs/rTg4510_RRBS_SmoothBetas.RData")))
rTg4510_array_beta <- get(load(file = paste0(dirnames$processed, "/array/Normalised_rTg4510_array_ECX_Final.RData")))
rTg4510_array_HIP_beta <- get(load(file = paste0(dirnames$processed, "/array/Normalised_rTg4510_array_HIP_Final.RData")))
J20_rrbs_beta <- get(load(file = paste0(dirnames$processed, "/rrbs/J20_RRBS_SmoothBetas.RData")))
J20_array_beta <- get(load(file = paste0(dirnames$processed, "/array/Normalised_J20_array_ECX_Final.RData")))
J20_array_HIP_beta <- get(load(file = paste0(dirnames$processed, "/array/Normalised_J20_array_HIP_Final.RData")))

# annotated RRBS positions
anno_rrbs_all <- get(load(file = paste0(dirnames$annotated, "/rrbs/rrbs_annoAllPositions.RData")))

# all results
rTg4510_rrbs_results <- get(load( file = paste0(dirnames$differential, "/rrbs/rTg4510_betaResultsDMPs.RData")))
J20_rrbs_results <- get(load( file = paste0(dirnames$differential, "/rrbs/J20_betaResultsDMPs.RData")))

# significant results
rTg4510_rrbs_sig <- get(load(file = paste0(dirnames$annotated, "/rrbs/rTg4510_rrbs_annoSigResultsDMPs.RData")))
rTg4510_array_sig <- get(load(file = paste0(dirnames$annotated, "/array/rTg4510_array_annoSigResultsDMPs.RData")))
J20_rrbs_sig <- get(load(file = paste0(dirnames$annotated, "/rrbs/J20_rrbs_annoSigResultsDMPs.RData")))
J20_array_sig <- get(load(file = paste0(dirnames$annotated, "/array/J20_array_annoSigResultsDMPs.RData")))

# brown results
rTg4510_brown_sig <- read.csv(paste0(dirnames$annotated,"/final/rTg4510_sigBrown.csv"))

# clock
load(file = paste0(dirnames$clock, "/rTg4510Clock.RData"))
load(file = paste0(dirnames$clock, "/J20Clock.RData"))

# significant array results in rTg4510 and J20
sigResArrayECX <- get(load(file = paste0(dirnames$annotated, "/final/entorhinalcortexArray_sigResultsDMPs.RData")))
sigResArrayHIP <- get(load(file = paste0(dirnames$annotated, "/final/hippocampusArray_sigResultsDMPs.RData")))

# common significant DMPs in Pathology and Interaction model
progressiveSites <- list(
  rTg4510 = intersect(sigRes$rTg4510$Pathology$Position, sigRes$rTg4510$GenotypeAge$Position),
  J20 = intersect(sigRes$J20$Pathology$Position, sigRes$J20$GenotypeAge$Position)
)
sigRes$rTg4510$PathologyCommonInteraction <- sigRes$rTg4510$Pathology[sigRes$rTg4510$Pathology$Position %in% progressiveSites$rTg4510, ]
sigRes$J20$PathologyCommonInteraction <- sigRes$J20$Pathology[sigRes$J20$Pathology$Position %in% progressiveSites$J20, ]

# common significant DMPs between ECX and HIP
commonECXHIPSites <- list(
  rTg4510_genotype = intersect(sigResArrayECX$rTg4510$Genotype$Position,sigResArrayHIP$rTg4510$Genotype$Position),
  rTg4510_interaction_pathology_ECX =  intersect(sigResArrayECX$rTg4510$GenotypeAge$Position, sigResArrayECX$rTg4510$Pathology$Position),
  rTg4510_interaction_pathology_HIP = intersect(sigResArrayHIP$rTg4510$GenotypeAge$Position, sigResArrayHIP$rTg4510$Pathology$Position),
  J20_interaction_pathology_HIP = intersect(sigResArrayHIP$J20$GenotypeAge$Position, sigResArrayHIP$J20$Pathology$Position)
)
hip_specific <- unique(setdiff(c(sigResArrayHIP$rTg4510$Genotype$Position, sigResArrayHIP$rTg4510$Pathology$Position),
                        c(sigResArrayECX$rTg4510$Genotype$Position, sigResArrayECX$rTg4510$Pathology$Position)))
sigResArrayHIP$rTg4510$PathologyCommonInteraction <- sigResArrayHIP$rTg4510$Pathology[sigResArrayHIP$rTg4510$Pathology$Position %in% commonECXHIPSites$rTg4510_interaction_pathology_HIP,]
sigResArrayHIP$J20$PathologyCommonInteraction <- sigResArrayHIP$J20$Pathology[sigResArrayHIP$J20$Pathology$Position %in% commonECXHIPSites$J20_interaction_pathology_HIP,]

## ------ human comparisons -----

# human significant DMPs
humanAllGeneList <- read.csv(paste0(dirnames$humanAnnot,"Shireby2022Meta_geneList.txt"), stringsAsFactors = F, header = F)[[1]]
humanDMPres <- read.csv(paste0(dirnames$humanAnnot, "GShireby2022_S8.csv"), header = T) %>% 
  mutate(UCSC.Nearest.Gene = gsub("ANK1;MIR486", "ANK1", UCSC.Nearest.Gene),
         UCSC.Nearest.Gene = gsub("CDH23;C10orf54", "CDH23", UCSC.Nearest.Gene))

# overlap genes
rTg4510HumanGenes <- intersect(c(sigRes$rTg4510$Genotype$ChIPseeker_GeneSymbol,sigRes$rTg4510$PathologyCommonInteraction$ChIPseeker_GeneSymbol),humanAllGeneList)
J20HumanGenes <- intersect(c(sigRes$J20$Genotype$ChIPseeker_GeneSymbol,sigRes$J20$PathologyCommonInteraction$ChIPseeker_GeneSymbol),humanAllGeneList)

# ECX rTg4510 and J20 overlap with human
sigRes$rTg4510_Human <- dplyr::bind_rows(sigRes$rTg4510$Genotype,sigRes$rTg4510$PathologyCommonInteraction) %>% 
  filter(ChIPseeker_GeneSymbol %in% rTg4510HumanGenes) 
sigRes$rTg4510_Human <- merge(sigRes$rTg4510_Human, mouse2human, all.x = T, by.x = "ChIPseeker_GeneSymbol", by.y = "mouse")

sigRes$J20_Human <- dplyr::bind_rows(sigRes$J20$Genotype,sigRes$J20$PathologyCommonInteraction) %>% 
  filter(ChIPseeker_GeneSymbol %in% J20HumanGenes) 
sigRes$J20_Human <- merge(sigRes$J20_Human, mouse2human, all.x = T, by.x = "ChIPseeker_GeneSymbol", by.y = "mouse")
