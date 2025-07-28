#!/usr/bin/env Rscript
## ----------Script-----------------  
##
## Purpose: Supplementary Figures
##         
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------

#-------------- Input -------------

scriptDir = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "0_PaperFigs/paper_import.config.R"))

## ------------ Figure 1: Overview -------
# generated via ppt

## ------------ Figure 2: rTg4510 -------

# Figure 2A: Manhattan plot of rtg4510 genotype (export plots as png (width = 1200, height = 260))
pManhattan <- list()
pManhattan$rTg4510_Genotype <- plot_manhattan(rTg4510_rrbs_results$Genotype, rTg4510_array_results$Genotype, mode = "Genotype")
pManhattan$rTg4510_Genotype

# Figure 2B: top-ranked DMPs in rTg4510 genotype 
Dcaf5 <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Dcaf5", "ENSMUST00000054145.7", boxplot = TRUE, colour = "rTg4510")
Arsi <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Arsi", "ENSMUST00000040359.5", colour = "rTg4510")
Creb3l4 <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Creb3l4", "ENSMUST00000029547.9", boxplot = TRUE, colour = "rTg4510")
As3mt <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "As3mt", "ENSMUST00000003655.8", colour = "rTg4510")

# Figure 2C: Manhattan plot of rTg4510 pathology (export plots as png (width = 1200, height = 260))
pManhattan$rTg4510_Pathology <- plot_manhattan(rTg4510_rrbs_results$Pathology, rTg4510_array_results$Pathology, mode = "Pathology")
pManhattan$rTg4510_Pathology

# Figure 2D: top-ranked DMPs in rTg4510 pathology 
Insyn2b <- plotGeneTrackDMP(sigRes$rTg4510$Pathology, sigBeta$rTg4510$Pathology, phenotype$rTg4510, "Insyn2b", "ENSMUST00000165963.8", colour = "rTg4510", boxplot = TRUE, pathology = TRUE)
Zfp423 <- plotGeneTrackDMP(sigRes$rTg4510$Pathology, sigBeta$rTg4510$Pathology, phenotype$rTg4510, "Zfp423", "ENSMUST00000109655.8", colour = "rTg4510", boxplot = TRUE, pathology = TRUE, position = "chr8:87750175")
Ankrd52 <- plotGeneTrackDMP(sigRes$rTg4510$Pathology, sigBeta$rTg4510$Pathology, phenotype$rTg4510, "Ankrd52", "ENSMUST00000014642.9", colour = "rTg4510", boxplot = TRUE, pathology = TRUE)
Adk <- plotGeneTrackDMP(sigRes$rTg4510$Pathology, sigBeta$rTg4510$Pathology, phenotype$rTg4510, "Adk", "ENSMUST00000045376.10", colour = "rTg4510", boxplot = TRUE, pathology = TRUE)
Cisd3 <- plotGeneTrackDMP(sigRes$rTg4510$Pathology, sigBeta$rTg4510$Pathology, phenotype$rTg4510, "Cisd3", "ENSMUST00000107584.7", colour = "rTg4510", boxplot = TRUE, pathology = TRUE)


## ------------ Figure 3: J20 -------

# Figure 3A: Manhattan plot of J20 genotype (export plots as png (width = 1200, height = 260))
pManhattan$J20_Genotype <- plot_manhattan(J20_rrbs_results$Genotype, J20_array_results$Genotype, mode = "Genotype")
pManhattan$J20_Genotype 

# Figure 3B: top-ranked DMPs in J20 genotype  
Nutf2 <- plotGeneTrackDMP(sigRes$J20$Genotype, sigBeta$J20$Genotype, phenotype$J20, "Nutf2", "ENSMUST00000008594.8", colour = "J20", boxplot = TRUE)
Tenm2 <- plotGeneTrackDMP(sigRes$J20$Genotype, sigBeta$J20$Genotype, phenotype$J20, "Tenm2", "ENSMUST00000102801.7", colour = "J20", boxplot = TRUE)

# Figure 3C: Manhattan plot of J20 pathology (export plots as png (width = 1200, height = 260))
pManhattan$J20_Pathology <- plot_manhattan(J20_rrbs_results$Pathology, J20_array_results$Pathology, mode = "Pathology")
pManhattan$J20_Pathology

# Figure 3D: top-ranked DMPs in J20 pathology 
Grk2 <- plotGeneTrackDMP(sigRes$J20$Pathology, sigBeta$J20$Pathology, phenotype$J20, "Grk2", "ENSMUST00000167511.2", colour = "J20", boxplot = TRUE, pathology = TRUE)
Fgfr2 <- plotGeneTrackDMP(sigRes$J20$Pathology, sigBeta$J20$Pathology, phenotype$J20, "Fgfr2", "ENSMUST00000117073.1", colour = "J20", boxplot = TRUE, pathology = TRUE)
Ncam2 <- plotGeneTrackDMP(sigRes$J20$PathologyCommonInteraction, sigBeta$J20$Pathology, phenotype$J20, "Ncam2", "ENSMUST00000037785.13", colour = "J20", pathology =  TRUE, boxplot = TRUE)
Zmiz1 <- plotGeneTrackDMP(sigRes$J20$PathologyCommonInteraction, sigBeta$J20$Pathology, phenotype$J20, "Zmiz1", "ENSMUST00000162645.7", colour = "J20", boxplot = TRUE, pathology = TRUE)


## ------------ Figure 4: ECX vs HIP -------

# Figure 4A: Venn diagram of hippocampus vs entorhinal cortex
HipECXVennrTg4510 <- plot_grid(venn.diagram(
  x = list(rTg4510_array_sig$ECX$Genotype$position, rTg4510_array_sig$ECX$Pathology$position,  
           rTg4510_array_sig$HIP$Genotype$position, rTg4510_array_sig$HIP$Pathology$position),
  category.names = c("ECX_Genotype" , "ECX_Pathology", "HIP_Genotype", "HIP_Pathology"),
  fill = pastelColours,
  cex = 0.9,
  cat.cex = 0.9, 
  filename = NULL
))

HipECXVennJ20 <- plot_grid(venn.diagram(
  x = list(J20_array_sig$ECX$Genotype$position, J20_array_sig$ECX$Pathology$position,  
           J20_array_sig$HIP$Genotype$position, J20_array_sig$HIP$Pathology$position),
  category.names = c("ECX_Genotype" , "ECX_Pathology", "HIP_Genotype", "HIP_Pathology"),
  fill = pastelColours,
  cex = 0.9,
  cat.cex = 0.9, 
  filename = NULL
))


# Figure 4B: Top-ranked DMP across rTg4510 ECX and HIP
Dennd1a = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                              ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr2:37946161", 
                              pathology = TRUE, gene = "Dennd1a")

Rapgefl1 = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                               ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr11:98838683", 
                               pathology = TRUE, gene = "Rapgefl1")


# Figure 4C: Top-ranked DMP across rTg4510 HIP but not ECX
HIPrTg4510plots <- list(
  Pxk = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                          ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr14:8146212", 
                          gene = "Pxk"),
  Mef2c = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                            ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr13:83504232", 
                            gene = "Mef2c"),
  Agbl5 = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                            ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr5:30890202", pathology = TRUE, 
                            gene = "Agbl5"),
  Meis2 = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                            ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr2:116018971", pathology = TRUE, 
                            gene = "Meis2")
)

# Figure 4D: Top-ranked DMP across J20 HIP but not ECX
HIPJ20plots <- list(
  Mir568  = plot_DMP_byTissue(ECXbetaMatrix=J20_array_beta, HIPbetaMatrix=J20_array_HIP_beta, 
                              ECXphenotypeFile=phenotype$J20, HIPphenotypeFile=phenotype$J20_HIP, position ="chr16:43609394", 
                              gene = "Mir568", model = "J20"),
  Mctp1 = plot_DMP_byTissue(ECXbetaMatrix=J20_array_beta, HIPbetaMatrix=J20_array_HIP_beta, 
                            ECXphenotypeFile=phenotype$J20, HIPphenotypeFile=phenotype$J20_HIP, position ="chr13:76810803", 
                            gene = "Mctp1", model = "J20"),
  Sox4 = plot_DMP_byTissue(ECXbetaMatrix=J20_array_beta, HIPbetaMatrix=J20_array_HIP_beta, 
                           ECXphenotypeFile=phenotype$J20, HIPphenotypeFile=phenotype$J20_HIP, position ="chr13:28949481", 
                           gene = "Sox4", model = "J20", pathology = TRUE),
  Cetn3 = plot_DMP_byTissue(ECXbetaMatrix=J20_array_beta, HIPbetaMatrix=J20_array_HIP_beta, 
                            ECXphenotypeFile=phenotype$J20, HIPphenotypeFile=phenotype$J20_HIP, position ="chr13:81828611", pathology = TRUE, 
                            gene = "Cetn3", model = "J20")
)


## ------------ Figure 5: Human comparison -------

# Figure 5A: Venn diagram of rTg4510 vs J20 vs Human  
sigRes$rTg4510 <- lapply(sigRes$rTg4510, function(x) x %>% filter(ChIPseeker_GeneSymbol != "NA"))
sigRes$J20 <- lapply(sigRes$J20, function(x) x %>% filter(ChIPseeker_GeneSymbol != "NA"))
pHuman1 <- venn.diagram(
  x = list(c(sigRes$rTg4510$Genotype$ChIPseeker_GeneSymbol,sigRes$rTg4510$Pathology$ChIPseeker_GeneSymbol),
           c(sigRes$J20$Genotype$ChIPseeker_GeneSymbol, sigRes$J20$Pathology$ChIPseeker_GeneSymbol), 
           humanAllGeneList),
  category.names = c("rTg4510","J20","Human"),
  fill = c(label_colour("rTg4510"), label_colour("J20"),"yellow"),
  filename = NULL
)


# Figure 5B: Ank1 RRBS and pyrosequencing
pAnk1DMP <- plot_gene_track(betaMatrix=sigBeta$rTg4510$Genotype, phenotypeFile=phenotype$rTg4510, 
                            position = "chr8:23023192", gene="Ank1", transcript="ENSMUST00000110688.8", colour = "rTg4510")

tAnk1DMP <- plot_DMP(betaMatrix=sigBeta$rTg4510$Genotype, phenotypeFile=phenotype$rTg4510, 
                     position = c("chr8:23023240","chr8:23023210","chr8:23023192"), table = TRUE) %>% mutate(method = "RRBS")

ank1PyroPos <- c(
  `Pos1Meth` = "chr8:23023192",
  `Pos3Meth` = "chr8:23023240"
)
ank1PyroPosdf <- reshape2::melt(ank1PyroPos, value.name = "Position") %>% tibble::rownames_to_column(., var = "prnpPosition")

tAnk1Pyro <- input_pyro$ank1 %>% dplyr::select(SAMPLE, Age, Group.ID, Pos1Meth, Pos3Meth) %>% 
  reshape2::melt(id = c("Age","Group.ID","SAMPLE"), variable.name = "Position", value.name = "methylation") %>% 
  mutate(Age = as.factor(stringr::str_remove(Age,"m"))) %>% 
  mutate(Group.ID = factor(Group.ID, levels = c("WT","TG"))) %>%
  merge(., phenotype$rTg4510, by.x = "SAMPLE", by.y = 0)%>% 
  merge(., reshape2::melt(ank1PyroPos, value.name = "position"), by.x = "Position", by.y = 0) %>% 
  dplyr::rename("sample"= "Position") %>% mutate(method = "Pyrosequencing") %>%
  mutate(methylation = methylation/100) 

pAnk1PyroRRBS <- rbind(tAnk1DMP, tAnk1Pyro %>% dplyr::select(colnames(tAnk1DMP))) %>% 
  mutate(method = factor(method, levels = c("RRBS","Pyrosequencing"))) %>%
  ggplot(., aes(x = Genotype, y = methylation, fill = Genotype)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = Genotype),width = 0.25, size = 2) +
  scale_fill_manual(values = c(alpha("black",0.2), color_Tg4510_TG),guide="none") +
  scale_colour_manual(values = c("black", color_Tg4510_TG),guide="none") +
  labs(x = "Genotype", y = "Methylation") +
  facet_nested(~ position + method) +
  mytheme + 
  theme(panel.border = element_rect(fill = NA, color = "grey", linetype = "dotted"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) 



## ------------ pdf outputs -------

pdf(paste0(output, "Figures/rTg4510_DMPs.pdf"), width = 21, height = 12)
plot_grid(Dcaf5, Arsi, Creb3l4, As3mt, scale = 0.95)
plot_grid(Cisd3, Zfp423, Adk, Insyn2b, scale = 0.95)
dev.off()

pdf(paste0(output, "Figures/J20_DMPs_2B.pdf"), width = 21, height = 8)
plot_grid(Nutf2, Tenm2, scale = 0.95)
dev.off()

pdf(paste0(output, "Figures/J20_DMPs_2D.pdf"), width = 21, height = 12)
plot_grid(Grk2, Fgfr2, Ncam2, Zmiz1, scale = 0.95)
dev.off()

pdf(paste0(output, "Figures/venn_HIP_ECX.pdf"), width = 10, height = 5)
plot_grid(HipECXVennrTg4510,HipECXVennJ20, scale = 0.85)
dev.off()

pdf(paste0(output, "Figures/rTg4510_ECX_HIP.pdf"), width = 21, height = 4)
plot_grid(Dennd1a, Rapgefl1, scale = 0.95, labels = c("i","ii"))
dev.off()

pdf(paste0(output, "Figures/rTg4510_HIP_notECX.pdf"), width = 21, height = 8)
plot_grid(plotlist = HIPrTg4510plots, scale = 0.95, labels = c("i","ii","iii","iv"), label_size = 18)
dev.off()

pdf(paste0(output, "Figures/J20_HIP_notECX.pdf"), width = 21, height = 8)
plot_grid(plotlist = HIPJ20plots, scale = 0.95, labels = c("i","ii","iii","iv"))
dev.off()

pdf(paste0(output, "Figures/Venn_human_comp.pdf"), width = 5, height = 5)
plot_grid(pHuman1, scale = 0.95)
dev.off()

pdf(paste0(output, "Figures/Ank1.pdf"), width = 12, height = 8)
plot_grid(pAnk1DMP, pAnk1PyroRRBS, ncol = 1, rel_heights = c(0.3,0.7))
dev.off()
