
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "0_PaperFigs/Functions.R"))
source(paste0(scriptDir, "0_PaperFigs/paper_import.config.R"))


message("Number of sites in RRBS in rTg4510: ", nrow(rTg4510_rrbs_beta))
message("Number of sites in RRBS in J20: ", nrow(J20_rrbs_beta))

length(unique(as.character(row.names(rTg4510_rrbs_beta))))
length(unique(as.character(row.names(rTg4510_array_beta))))
length(intersect(unique(as.character(row.names(rTg4510_array_beta))),unique(as.character(row.names(rTg4510_rrbs_beta)))))
message("Number of sites (RRBS and array) in rTg4510:")
(1347524 + 23588) - (3541 * 2)

length(unique(as.character(row.names(J20_rrbs_beta))))
length(unique(as.character(row.names(J20_array_beta))))
length(intersect(unique(as.character(row.names(J20_array_beta))),unique(as.character(row.names(J20_rrbs_beta)))))
message("Number of sites (RRBS and array) in rTg4510:")
(1372261 + 23588) - (3557 * 2)

length(unique(as.character(row.names(J20_array_HIP_beta))))
length(unique(as.character(row.names(rTg4510_array_HIP_beta))))

message("Mean number of reads in J20: ", mean(rawRRBSReads$J20$No..raw.reads))
message("Mean number of reads in rTg4510: ", mean(rawRRBSReads$rTg4510$No..raw.reads))
message("SD number of reads in J20: ", sd(rawRRBSReads$J20$No..raw.reads))
message("SD number of reads in rTg4510: ", sd(rawRRBSReads$rTg4510$No..raw.reads))

message("t-test of read difference between WT and TG in J20")
t.test(No..raw.reads ~ Genotype, rawRRBSReads$J20)
t.test(No..raw.reads ~ Genotype, rawRRBSReads$rTg4510)


## ------ epigenetic clock ------ 

# correlation in rTg4510
cor.test(rTg4510Clocks$ECX$Age_months, rTg4510Clocks$ECX$DNAmAgeClockCortex)
cor.test(J20Clocks$ECX$Age_months, J20Clocks$ECX$DNAmAgeClockCortex)

## ------ changes in ECX: Genotype + Interaction effect 

message("Number of significant DMPs (FDR < 0.05) genotype: ", nrow(sigRes$rTg4510$Genotype))
message("Number of genes associated with significant DMPs genotype: ", length(unique(sigRes$rTg4510$Genotype$ChIPseeker_GeneSymbol)))
message("Number of significant DMPs (FDR < 0.05) genotype, hypomethylated: ", nrow(sigRes$rTg4510$Genotype %>% filter(directionTG == "down")))
message("Number of significant DMPs (FDR < 0.05) genotype, hypermethylated: ", nrow(sigRes$rTg4510$Genotype %>% filter(directionTG == "up")))
binom.test(686, 1143, p = 0.5, alternative = "two.sided")
message("Number of significant DMPs (FDR < 0.05) pathology: ", nrow(sigRes$rTg4510$Pathology))
message("Number of significant DMPs (FDR < 0.05) unique to pathology: ", length(setdiff(sigRes$rTg4510$Pathology$Position, sigRes$rTg4510$Genotype$Position)))
message("Number of common significant DMPs pathology and interaction: ", length(progressiveSites$rTg4510))
message("Number of additional sites in DMPs pathology and interactions vs genotype: ", length(intersect(sigRes$rTg4510$Genotype$Position, progressiveSites$rTg4510)))
message("Number of genes associated with significant DMPs interaction and pathology: ", length(unique(sigRes$rTg4510$PathologyCommonInteraction$ChIPseeker_GeneSymbol)))


message("Number of significant DMPs (FDR < 0.05) genotype: ", nrow(sigRes$J20$Genotype))
message("Number of significant DMPs (FDR < 0.05) genotype, hypomethylated: ", nrow(sigRes$J20$Genotype %>% filter(directionTG == "down")))
message("Number of significant DMPs (FDR < 0.05) genotype, hypermethylated: ", nrow(sigRes$J20$Genotype %>% filter(directionTG == "up")))
res <- binom.test(1673, 2521, p = 0.5, alternative = "two.sided")
res$p.value

message("Number of significant DMPs (FDR < 0.05) pathology: ", nrow(sigRes$J20$Pathology))
message("Number of significant DMPs (FDR < 0.05) unique to pathology: ", length(setdiff(sigRes$J20$Pathology$Position, sigRes$J20$Genotype$Position)))
message("Number of significant DMPs (FDR < 0.05) common to interaction and pathology: ", length(progressiveSites$J20))
message("Number of significant DMPs (FDR < 0.05) common to interaction and pathology but not genotype: ", length(setdiff(progressiveSites$J20, sigRes$J20$Genotype$Position)))


message("Common number of sites between genotype in J20 and rTg4510: ", length(intersect(sigRes$rTg4510$Genotype$Position, sigRes$J20$Genotype$Position)))
sigRes$rTg4510$Genotype[sigRes$rTg4510$Genotype$Position %in% intersect(sigRes$rTg4510$Genotype$Position, sigRes$J20$Genotype$Position),]

message("Common number of sites between pathology in J20 and rTg4510: ", length(intersect(progressiveSites$rTg4510, progressiveSites$J20)))

# effect size after removing transgene effects
message("mean delta difference in rTg4510 - hypo,", 
        median(sigRes_transgeneRem$rTg4510$Genotype[sigRes_transgeneRem$rTg4510$Genotype$directionTG == "down", "delta"]))
message("mean delta difference in J20 - hypo,", median(sigRes_transgeneRem$J20$Genotype[sigRes_transgeneRem$J20$Genotype$directionTG == "down", "delta"]))

message("mean delta difference in rTg4510 - hyper,", 
        median(sigRes_transgeneRem$rTg4510$Genotype[sigRes_transgeneRem$rTg4510$Genotype$directionTG == "up", "delta"]))
message("mean delta difference in J20 - hyper,",median(sigRes_transgeneRem$J20$Genotype[sigRes_transgeneRem$J20$Genotype$directionTG == "up", "delta"]))

# top-ranked significant by genotype
sigRes$rTg4510$Genotype[1,]
Prn <- sigRes$rTg4510$Genotype %>% filter(ChIPseeker_GeneSymbol == "Prn")
mean(Prn$delta)
rTg4510_brown_sig[rTg4510_brown_sig$Gene == "Prn",]
sigRes$rTg4510$Genotype %>% filter(ChIPseeker_GeneSymbol == "Prnp")

## ------ changes in Hippocampus vs ECX (array)

message("Total number of signficant DMPs associated with tau in rTg4510 TG vs WT mice (genotype) in hippocampus: ", length(rTg4510_array_sig$HIP$Genotype$position))
message("Number of signficant genotype in hippocampus that were significant in ECX rTg4510: ", length(commonECXHIPSites$rTg4510_genotype))
message("Signficiant genotype-associated DMPs in ECX and HIP")
rTg4510_array_sig$HIP$Genotype[rTg4510_array_sig$HIP$Genotype$position %in% commonECXHIPSites$rTg4510_genotype,]
message("Number of unique genotype DMPs in hippocampus:", nrow(sigResArrayHIP$rTg4510$Genotype[sigResArrayHIP$rTg4510$Genotype$Position %in% hip_specific,]))
message("Number of genes associated with unique genotype DMPs in hippocampus: ", 
        length(unique(sigResArrayHIP$rTg4510$Genotype[sigResArrayHIP$rTg4510$Genotype$Position %in% hip_specific,"ChIPseeker_GeneSymbol"])))
message("Number of unique pathology DMPs in hippocampus:", 
        nrow(sigResArrayHIP$rTg4510$PathologyCommonInteraction[sigResArrayHIP$rTg4510$PathologyCommonInteraction$Position %in% hip_specific,]))
message("Number of genes associated with unique genotype DMPs in hippocampus: ", 
        length(unique(sigResArrayHIP$rTg4510$PathologyCommonInteraction[sigResArrayHIP$rTg4510$PathologyCommonInteraction$Position %in% hip_specific,"ChIPseeker_GeneSymbol"])))

message("Number of common sites between pathology in hippocampus and genotye+pathology in entorhainl cortex: ", 
length(intersect(sigResArrayHIP$rTg4510$PathologyCommonInteraction$Position, 
                 c(sigResArrayECX$rTg4510$Genotype$Position,                                                                       
                   intersect(sigResArrayECX$rTg4510$GenotypeAge$Position, sigResArrayECX$rTg4510$Pathology$Position)))))

message("Total number of signficant DMPs associated with amyloid in J20 TG vs WT mice (genotype) in hippocampus: ", length(J20_array_sig$HIP$Genotype$position))
message("Total number of signficant DMPs associated with amyloid in J20 TG vs WT mice (pathology) in hippocampus: ", length(sigResArrayHIP$J20$PathologyCommonInteraction$Position))

message("Total number of DMPs from array in rTg4510, ECX & HIP: ", length(unique(c(sigResArrayECX$rTg4510$Genotype$Position,commonECXHIPSites$rTg4510_interaction_pathology_ECX,  
  sigResArrayHIP$rTg4510$Genotype$Position, commonECXHIPSites$rTg4510_interaction_pathology_HIP))))

message("Total number of DMPs from array in J20, ECX & HIP: ", length(unique(c(J20_array_sig$ECX$Genotype$position,intersect(J20_array_sig$ECX$Interaction$position, J20_array_sig$ECX$Pathology$position), J20_array_sig$HIP$Genotype$position, intersect(J20_array_sig$HIP$Interaction$position, J20_array_sig$HIP$Pathology$position)))))


## ------ gene ontology results

write.table(sigRes$J20$Genotype$ChIPseeker_GeneSymbol, paste0(dirnames$differential,"/merged/J20_Genotype_sigGenes.txt"), quote = F, row.names = F, col.names = F)
write.table(sigRes$J20$Pathology$ChIPseeker_GeneSymbol, paste0(dirnames$differential,"/merged/J20_Pathology_sigGenes.txt"), quote = F, row.names = F, col.names = F)
write.table(sigRes$rTg4510$Genotype$ChIPseeker_GeneSymbol, paste0(dirnames$differential,"/merged/rTg4510_Genotype_sigGenes.txt"), quote = F, row.names = F, col.names = F)
write.table(sigRes$rTg4510$Pathology$ChIPseeker_GeneSymbol, paste0(dirnames$differential,"/merged/rTg4510_Pathology_sigGenes.txt"), quote = F, row.names = F, col.names = F)

commonDMPGenesPathologyBetweenModels <- unique(intersect(sigRes$J20$Pathology$ChIPseeker_GeneSymbol, sigRes$rTg4510$Pathology$ChIPseeker_GeneSymbol))
length(commonDMPGenesPathologyBetweenModels)
write.table(commonDMPGenesPathologyBetweenModels, paste0(dirnames$differential,"/merged/ComparisonModel_Pathology_sigGenes.txt"), quote = F, row.names = F, col.names = F)
