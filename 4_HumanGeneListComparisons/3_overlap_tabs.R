library("stringr")
library("knitr")
library(ggrepel)
output_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/Figures/"

#read in results file from meta analysis
rootDir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/"
outputHumanListDir <- "/lustre/projects/Research_Project-MRC148213/lsl693/rrbs_ad_mice/0_ZenOutput/4_humanAnnotations/"
humanTauGeneList <- read.csv(paste0(outputHumanListDir,"Shireby2022Braak_geneList.txt"), stringsAsFactors = F, header = F)[[1]]
humanAmyloidGeneList <- read.csv(paste0(outputHumanListDir,"Shireby2022Thal_geneList.txt"), stringsAsFactors = F, header = F)[[1]]
humanAllGeneList <- read.csv(paste0(outputHumanListDir,"Shireby2022Meta_geneList.txt"), stringsAsFactors = F, header = F)[[1]]
humanDMPres <- read.csv(paste0(outputHumanListDir, "GShireby2022_S8.csv"), header = T)

#load in files
rTg4510DMPRes <- list.files(path = paste0(dirnames$paper,"/tables"), pattern="rTg4510",full.names = T)
J20DMPRes <- list.files(path = paste0(dirnames$paper,"/tables"), pattern="J20",full.names = T)

#create gene lists

identifyHumanOverlap <- function(temp,converted){
  overlap_list <- list()
  overlap_df <- list()
  
  #create summary table to hold information from for loop
  tab <- data.frame(matrix(ncol = 6, nrow = 6))
  colnames(tab) <- c("Model.Term", "N.mouse.probes", "N.mouse.genes","N.mouse.probes.overlaps", "N.mouse.gene.overlaps","N.human.gene.overlaps")
  tab$Model.Term <- c("ECX_Genotype", "ECX_Interaction", "ECX_Pathology", "HIP_Genotype", "HIP_Pathology", "TissueArray")
  
  for(i in 1:length(temp)){
    
    # read in results file
    file <- read.csv(temp[i], stringsAsFactors = F)
    
    # create list of genes
    RRBS_genes <- unique(file$ChIPseeker_GeneSymbol)
    
    chipseeker_overlap <- intersect(RRBS_genes, converted)
    overlap_list[[i]] <- sort(chipseeker_overlap)
    overlap_df[[i]] <- file[file$ChIPseeker_GeneSymbol%in% chipseeker_overlap,]
    
    #message("Number of mouse probes per common gene")
    #print(overlap_df[[i]] %>% group_by(ChIPseeker_GeneSymbol) %>% tally())
    
    #populate table
    tab[i, "N.mouse.probes"] <- nrow(file)
    tab[i, "N.mouse.genes"] <- length(RRBS_genes)
    tab[i, "N.mouse.probes.overlaps"] <- nrow(file[file$ChIPseeker_GeneSymbol %in% chipseeker_overlap,])
    tab[i,"N.mouse.gene.overlaps"] <- paste0(length(chipseeker_overlap)," (", round(length(chipseeker_overlap)/length(RRBS_genes) * 100,2), "%)")
    tab[i,"N.human.gene.overlaps"] <- paste0(length(chipseeker_overlap)," (", round(length(chipseeker_overlap)/length(converted) * 100,2), "%)")
    
  }
  
  print(kable(tab, caption = "RRBS"))
  
  names(overlap_df) <- c("ECX_Genotype", "ECX_Interaction", "ECX_Pathology", "HIP_Genotype", "HIP_Pathology", "TissueArray")
  
  # meth.diff.Genotype = meth.group2.TG - meth.group1.WT 
  # therefore negative refers to higher methylation in TG mouse
  #overlap_df$genotype <- overlap_df$genotype %>% mutate(direction = ifelse(meth.diff.Genotype > 0, "DOWN","UP"))
  #overlap_df$interaction <- overlap_df$interaction %>% mutate(direction = ifelse(meth.diff.Genotype > 0, "DOWN","UP"))
  
  output <- list(tab, overlap_df)
  names(output) <- c("tab","overlap_df")
  return(output)
}


identifyConsistentES <- function(overlapRes, humanDMPres){
  ES <- list()
  human_overlap_df <- list()

  overlapRes$human_overlap_df <- lapply(overlapRes$overlap_df, function(x) 
    dplyr::filter(humanDMPres, grepl(paste(toupper(unique(x$ChIPseeker_GeneSymbol)), collapse = "|"),(UCSC.Nearest.Gene))) 
    %>% mutate(UCSC.Nearest.Gene = word(UCSC.Nearest.Gene,c(1),sep=fixed(";"))))
  
  
  for(i in names(overlapRes$human_overlap_df)[1:5]){
    mouse_mod <- overlapRes$overlap_df[[i]] %>% dplyr::select(Position,contains("FDR_adj"),delta, ChIPseeker_GeneSymbol) %>% 
      mutate(toUpperGeneSymbol = toupper(ChIPseeker_GeneSymbol))
    human_mod <- overlapRes$human_overlap_df[[i]] %>% dplyr::select(DNAm.site,Effect_fixed....,P_fixed,UCSC.Nearest.Gene)
    
    ES[[i]] <-  mouse_mod %>% 
      full_join(., human_mod, by = c("toUpperGeneSymbol" = "UCSC.Nearest.Gene"),relationship = "many-to-many") %>%
      `colnames<-`(c("mousePosition","mouseFDRGenotype","mouseESGenotype","mouseGene","humanGene","humanProbe","humanES","humanPvalue")) %>%
      mutate(mouseDir = ifelse(mouseESGenotype < 0, "DOWN","UP"), humanDir = ifelse(humanES < 0, "DOWN","UP")) %>% 
      mutate(consistentDir = ifelse(mouseDir == humanDir, TRUE,FALSE))
  }
  #write.csv(ES$genotype,paste0(output_dir,"human_mouse_genotypeEffectSizeComparison.csv"), quote=F, row.names=F)
  
  for(i in 1:length(ES)){
    
    consistentGene <- unique(ES[[i]] %>% filter(consistentDir == TRUE) %>% .[,c("mouseGene")])
    consistentProbes <- unique(ES[[i]] %>% filter(consistentDir == TRUE) %>% .[,c("mousePosition")])
    print(consistentGene)
    overlapRes$tab[i, "N.mouse.gene.overlaps.consistentES"] <- length(consistentGene)
    overlapRes$tab[i, "N.mouse.probe.overlaps.consistentES"] <- length(consistentProbes)
    
  }
  print(overlapRes$tab)
  
  output <- list(overlapRes, ES)
  names(output) <- c("overlapRes", "ES")

  return(output)
}

rTg4510Overlap <- identifyHumanOverlap(rTg4510DMPRes, humanAllGeneList)
J20Overlap <- identifyHumanOverlap(J20DMPRes, humanAllGeneList)

rTg4510OverlapConsistent <- identifyConsistentES(rTg4510Overlap, humanDMPres)

# plots
ggplot(rTg4510OverlapConsistent$ES$ECX_Genotype, aes(x = mouseESGenotype, y = humanES, colour = humanGene, label = mousePosition)) + 
  geom_point(size = 2) + 
  geom_text_repel(show.legend = FALSE) +
  geom_hline(yintercept=0, linetype = "dotted") +
  geom_vline(xintercept=0, linetype = "dotted") +
  theme_classic() +
  labs(x = "Effect size in rTg4510 TG vs WT", y = "Effect size in human AD studies") +
  scale_colour_discrete(name = "Gene")

ggplot(rTg4510OverlapConsistent$ES$ECX_Interaction, aes(x = mouseESGenotype, y = humanES, colour = humanGene, label = mousePosition)) + 
  geom_point(size = 2) + 
  geom_text_repel(show.legend = FALSE) +
  geom_hline(yintercept=0, linetype = "dotted") +
  geom_vline(xintercept=0, linetype = "dotted") +
  theme_classic() +
  labs(x = "Effect size in rTg4510 TG vs WT", y = "Effect size in human AD studies") +
  scale_colour_discrete(name = "Gene")

ggplot(rTg4510OverlapConsistent$ES$ECX_Pathology, aes(x = mouseESGenotype, y = humanES, colour = humanGene, label = mousePosition)) + 
  geom_point(size = 2) + 
  geom_text_repel(show.legend = FALSE) +
  geom_hline(yintercept=0, linetype = "dotted") +
  geom_vline(xintercept=0, linetype = "dotted") +
  theme_classic() +
  labs(x = "Effect size in rTg4510 TG vs WT", y = "Effect size in human AD studies") +
  scale_colour_discrete(name = "Gene")

unique(rTg4510OverlapConsistent$ES$ECX_Genotype %>% filter(consistentDir == TRUE) %>% .[["humanGene"]])
cor.test(rTg4510OverlapConsistent$ES$ECX_Genotype$mouseESGenotype,rTg4510OverlapConsistent$ES$ECX_Genotype$humanES)


