
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
        median(sigRes_transgeneRem$rTg4510$Genotype[sigRes_transgeneRem$rTg4510$Genotype$directionTG == "down", "BetaSize_Genotype"]))
message("mean delta difference in J20 - hypo,", median(sigRes_transgeneRem$J20$Genotype[sigRes_transgeneRem$J20$Genotype$directionTG == "down", "BetaSize_Genotype"]))

message("mean delta difference in rTg4510 - hyper,", 
        median(sigRes_transgeneRem$rTg4510$Genotype[sigRes_transgeneRem$rTg4510$Genotype$directionTG == "up", "BetaSize_Genotype"]))
message("mean delta difference in J20 - hyper,",median(sigRes_transgeneRem$J20$Genotype[sigRes_transgeneRem$J20$Genotype$directionTG == "up", "BetaSize_Genotype"]))

comparison_effect_size <- rbind(sigRes$rTg4510$Genotype[,c("Platform","BetaSize_Genotype")] %>% mutate(model = "rTg4510"),
                                sigRes$J20$Genotype[,c("Platform","BetaSize_Genotype")] %>% mutate(model = "J20")) 

# p-value of effect sizes
down_effect_size <- comparison_effect_size[comparison_effect_size$BetaSize_Genotype < 0,]
with(down_effect_size, shapiro.test(BetaSize_Genotype[model == "rTg4510"])) # not normal
with(down_effect_size, shapiro.test(BetaSize_Genotype[model == "J20"])) # not normal
res <- wilcox.test(down_effect_size[down_effect_size$model == "rTg4510","BetaSize_Genotype"], 
                   down_effect_size[down_effect_size$model == "J20","BetaSize_Genotype"])
res$p.value

up_effect_size <- comparison_effect_size[comparison_effect_size$BetaSize_Genotype >= 0,]
with(up_effect_size, shapiro.test(BetaSize_Genotype[model == "rTg4510"])) # not normal
with(up_effect_size, shapiro.test(BetaSize_Genotype[model == "J20"])) # not normal
res <- wilcox.test(up_effect_size[up_effect_size$model == "rTg4510","BetaSize_Genotype"], 
                   up_effect_size[up_effect_size$model == "J20","BetaSize_Genotype"])
res$p.value
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

# FDR values for all ECX array  
rTg4510_array_results$Genotype$FDR_adj_Genotype <- p.adjust(rTg4510_array_results$Genotype[,"PrZ.GenotypeTG"], method = "fdr")
rTg4510_array_results$Pathology$FDR_adj_Pathology <- p.adjust(rTg4510_array_results$Pathology[,"PrZ.Pathology"], method = "fdr")
J20_array_results$Genotype$FDR_adj_Genotype <- p.adjust(J20_array_results$Genotype[,"PrZ.GenotypeTG"], method = "fdr")
J20_array_results$Pathology$FDR_adj_Pathology <- p.adjust(J20_array_results$Pathology[,"PrZ.Pathology"], method = "fdr")

# Pxk
sigResArrayHIP$rTg4510$Genotype[sigResArrayHIP$rTg4510$Genotype$Position == "chr14:8146212",]
rTg4510_array_results$Genotype[rTg4510_array_results$Genotype$X == "cg04016300",]
# Mef2c
sigResArrayHIP$rTg4510$Genotype[sigResArrayHIP$rTg4510$Genotype$Position == "chr13:83504232",]
rTg4510_array_results$Genotype[rTg4510_array_results$Genotype$X == "cg05512265",]
# Faf1
sigResArrayHIP$rTg4510$Pathology[sigResArrayHIP$rTg4510$Pathology$Position == "chr5:30890202",]
rTg4510_array_results$Pathology[rTg4510_array_results$Pathology$X == "cg17903645",]
# Meis2
sigResArrayHIP$rTg4510$Pathology[sigResArrayHIP$rTg4510$Pathology$Position == "chr2:116018971",]
rTg4510_array_results$Pathology[rTg4510_array_results$Pathology$X == "cg05692692",]
# Mir568 
sigResArrayHIP$J20$Genotype[sigResArrayHIP$J20$Genotype$Position == "chr16:43609394",]
J20_array_results$Genotype[J20_array_results$Genotype$X == "cg06690459",]
# Mcpt2
sigResArrayHIP$J20$Genotype[sigResArrayHIP$J20$Genotype$Position == "chr13:76810803",]
J20_array_results$Genotype[J20_array_results$Genotype$X == "cg02533269",]
# Sox4
sigResArrayHIP$J20$Pathology[sigResArrayHIP$J20$Pathology$Position == "chr13:28949481",]
J20_array_results$Pathology[J20_array_results$Pathology$X == "cg22103215",]
# Cetn3
sigResArrayHIP$J20$Pathology[sigResArrayHIP$J20$Pathology$Position == "chr13:81828611",]
J20_array_results$Pathology[J20_array_results$Pathology$X == "cg22103215",]


## ----- human ----

message("common rTg4510, J20 and human gene: ", intersect(rTg4510HumanGenes,J20HumanGenes))

## ------ gene ontology results

write.table(sigRes$J20$Genotype$ChIPseeker_GeneSymbol, paste0(dirnames$differential,"/merged/J20_Genotype_sigGenes.txt"), quote = F, row.names = F, col.names = F)
write.table(sigRes$J20$Pathology$ChIPseeker_GeneSymbol, paste0(dirnames$differential,"/merged/J20_Pathology_sigGenes.txt"), quote = F, row.names = F, col.names = F)
write.table(sigRes$rTg4510$Genotype$ChIPseeker_GeneSymbol, paste0(dirnames$differential,"/merged/rTg4510_Genotype_sigGenes.txt"), quote = F, row.names = F, col.names = F)
write.table(sigRes$rTg4510$Pathology$ChIPseeker_GeneSymbol, paste0(dirnames$differential,"/merged/rTg4510_Pathology_sigGenes.txt"), quote = F, row.names = F, col.names = F)

commonDMPGenesPathologyBetweenModels <- unique(intersect(sigRes$J20$Pathology$ChIPseeker_GeneSymbol, sigRes$rTg4510$Pathology$ChIPseeker_GeneSymbol))
length(commonDMPGenesPathologyBetweenModels)
write.table(commonDMPGenesPathologyBetweenModels, paste0(dirnames$differential,"/merged/ComparisonModel_Pathology_sigGenes.txt"), quote = F, row.names = F, col.names = F)
## ----- effect sizes ----

filter_for_effectsize <- function(res, position = NULL, gene = NULL){
  
  if(!is.null(position)){
    dat <-res[res$Position == position,] 
    print(dat) 
    message("Beta effect size for ", position, " : ",  round(as.numeric(as.character(dat %>% dplyr::select(contains("BetaSize")))),2))
  }else if(!is.null(gene)){
    dat <- res %>% filter(ChIPseeker_GeneSymbol == gene)
    print(dat)
    message("Number of DMPs associated to gene: ", nrow(dat))
    message("Mean beta effect size: ", round(mean(dat %>% dplyr::pull(contains("Beta"))),2))
  }
}

# Mapt
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Genotype, position = "chr11:104318231")
# Prn
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Genotype, gene = "Prn")
# rTg4510 genotype
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Genotype, position = "chr12:80436248")
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Genotype, gene = "Ncapg2")
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Genotype, gene = "Fgf14")
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Genotype, gene = "Arsi")
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Genotype, gene = "Creb3l4")
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Genotype, gene = "As3mt")
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Genotype, gene = "Ank1")
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Genotype, gene = "Prdm16")
filter_for_effectsize(J20_ECX_sigResultsDMPs$Genotype, gene = "Nutf2")
filter_for_effectsize(J20_ECX_sigResultsDMPs$Genotype, gene = "Tenm2")
# rTg4510 pathology
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Pathology, position = "chr11:97685376")
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Pathology, position = "chr8:87750175")
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Pathology, position = "chr14:21369243")
filter_for_effectsize(rTg4510_ECX_sigResultsDMPs$Pathology, position = "chr11:34314785")
# J20 pathology
filter_for_effectsize(J20_ECX_sigResultsDMPs$Pathology, position = "chr16:81200030")
filter_for_effectsize(J20_ECX_sigResultsDMPs$Pathology, position = "chr5:38684760")
filter_for_effectsize(J20_ECX_sigResultsDMPs$Pathology, gene = "Prmt8")
filter_for_effectsize(J20_ECX_sigResultsDMPs$Pathology, gene = "Zmiz1")
filter_for_effectsize(J20_ECX_sigResultsDMPs$Genotype, position = "chr14:40966816")
filter_for_effectsize(J20_ECX_sigResultsDMPs$Genotype, position = "chr4:154519364")
