#!/usr/bin/env Rscript
## ----------Script-----------------  
##
## Purpose: Supplementary Figures
##         
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------


#-------------- input -------------

rootDir = "/lustre/projects/Research_Project-191406"
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "0_PaperFigs/Functions.R"))
source(paste0(scriptDir, "0_PaperFigs/paper_import.config.R"))
source(paste0(scriptDir, "3_ArrayRRBSComparison/functions/summaryStatsDMP.R"))
source(paste0(LOGEN_ROOT, "/aesthetics_basics_plots/draw_venn.R"))

# old results
#res <- read.csv("/lustre/projects/Research_Project-191406/EmmaW/RRBSAnnotatedResults/rTg4510/DMPs/DMPsPathology_rTg4510_sig_pathology_1500bptssAn#no.csv")
#res <- res %>% arrange(FDR_adj_pathology)
#load(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/rTg4510_RRBSbetasComplete.RData")
#plot_DMP(betaMatrix=RRBS_completebetas, phenotypeFile=phenotype$rTg4510, position = "chr11:6228592", pathology = TRUE) 


#-------------- Figures -------------

# plot RRBS CpG sites by annotation
plot_annotate_sites()

# correlation of probes for RRBS vs Array
# rTg4510
rTg4510_common_probes <- commonStatsDescription(rTg4510_rrbs_beta, rTg4510_array_beta)
rTg4510_corr_arrayRRBS <- corrPlotCommonProbes(rTg4510_rrbs_beta, rTg4510_array_beta, rTg4510_common_probes)
# J20
J20_common_probes <- commonStatsDescription(J20_rrbs_beta, J20_array_beta)
J20_corr_arrayRRBS <- corrPlotCommonProbes(J20_rrbs_beta, J20_array_beta, J20_common_probes)
plot_grid(rTg4510_corr_arrayRRBS, J20_corr_arrayRRBS)

# heatmap of top 1000 genotype and pathology
pCluster = list(
  rTg4510Genotype = cluster_DMP(model="rTg4510", arrayBeta = rTg4510_array_beta, rrbsBeta = rTg4510_rrbs_beta, 
                                phenotypeInput=phenotype$rTg4510, lstPositions=sigRes$rTg4510$Genotype$Position),
  rTg4510Pathology = cluster_DMP(model="rTg4510", arrayBeta = rTg4510_array_beta, rrbsBeta = rTg4510_rrbs_beta, 
                                 phenotypeInput=phenotype$rTg4510, lstPositions=sigRes$rTg4510$Pathology$Position),
  J20Genotype = cluster_DMP(model="J20", arrayBeta = J20_array_beta, rrbsBeta = J20_rrbs_beta, 
                            phenotypeInput=phenotype$J20, lstPositions=sigRes$J20$Genotype$Position),
  J20Pathology = cluster_DMP(model="J20", arrayBeta = J20_array_beta, rrbsBeta = J20_rrbs_beta, 
                             phenotypeInput=phenotype$J20, lstPositions=sigRes$J20$Pathology$Position)
)
plot_grid(pCluster$rTg4510Genotype$gtable, pCluster$rTg4510Pathology$gtable)
plot_grid(pCluster$J20Genotype$gtable, pCluster$J20Pathology$gtable)

# magniture of effect size in J20 vs rTg4510  
ggplot(comparison_effect_size, aes(x = BetaSize_Genotype, fill = model)) + geom_density(alpha = 0.3) +
  theme_classic() +
  labs(x = "Effect size (Genotype)", y = "Density", fill = "Mouse model")

# manhattan plots 
pManhattan <- list(
  rTg4510Genotype = plot_manhattan_final(sigRes$rTg4510, "Genotype"),
  rTg4510Pathology = plot_manhattan_final(sigRes$rTg4510, "Pathology"),
  J20Genotype = plot_manhattan_final(sigRes$J20, "Genotype"),
  J20Pathology = plot_manhattan_final(sigRes$J20, "Pathology")
)
plot_grid(pManhattan$rTg4510Genotype,pManhattan$rTg4510Pathology)
plot_grid(pManhattan$J20Genotype,pManhattan$J20Pathology)

# venn diagram of positions
plot_grid(twovenndiagrams(sigRes$rTg4510$Genotype$Position, sigRes$rTg4510$Pathology$Position, "Genotype","Pathology"))
plot_grid(twovenndiagrams(sigRes$J20$Genotype$Position, sigRes$J20$Pathology$Position, "Genotype","Pathology"))

# top-ranked DMPs in rTg4510 genotype due to transgene
Mapt <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Mapt", "ENSMUST00000100347.10", boxplot = TRUE, colour = "rTg4510")
Prnp <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, c("Prn","Prnp"), "ENSMUST00000091288.12", colour = "rTg4510")
Fgf14 <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Fgf14", "ENSMUST00000095529.9", boxplot = TRUE, colour = "rTg4510")
Ncapg2 <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Ncapg2", "ENSMUST00000084828.4", boxplot = TRUE, colour = "rTg4510")

# top-ranked DMPs in rTg4510 genotype 
Dcaf5 <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Dcaf5", "ENSMUST00000054145.7", boxplot = TRUE, colour = "rTg4510")
Arsi <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Arsi", "ENSMUST00000040359.5", colour = "rTg4510")
Ugt2b37 <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Ugt2b37", "ENSMUST00000075858.3", colour = "rTg4510")
Creb3l4 <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Creb3l4", "ENSMUST00000029547.9", boxplot = TRUE, colour = "rTg4510")
As3mt <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "As3mt", "ENSMUST00000003655.8", colour = "rTg4510")

# top-ranked DMPs in rTg4510 pathology 
Insyn2b <- plotGeneTrackDMP(sigRes$rTg4510$Pathology, sigBeta$rTg4510$Pathology, phenotype$rTg4510, "Insyn2b", "ENSMUST00000165963.8", colour = "rTg4510", boxplot = TRUE, pathology = TRUE)
Zfp423 <- plotGeneTrackDMP(sigRes$rTg4510$Pathology, sigBeta$rTg4510$Pathology, phenotype$rTg4510, "Zfp423", "ENSMUST00000109655.8", colour = "rTg4510", boxplot = TRUE, pathology = TRUE, position = "chr8:87750175")
Ankrd52 <- plotGeneTrackDMP(sigRes$rTg4510$Pathology, sigBeta$rTg4510$Pathology, phenotype$rTg4510, "Ankrd52", "ENSMUST00000014642.9", colour = "rTg4510", boxplot = TRUE, pathology = TRUE)
Adk <- plotGeneTrackDMP(sigRes$rTg4510$Pathology, sigBeta$rTg4510$Pathology, phenotype$rTg4510, "Adk", "ENSMUST00000045376.10", colour = "rTg4510", boxplot = TRUE, pathology = TRUE)

# plot DMP and Track for Prnp
pPrnPrnPDMP <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, c("Prn","Prnp"), "ENSMUST00000091288.12")

# input pyrosequencing results
prnpPyroPos <- c(
  `Pos1Meth` = "chr2:131910162",
  `Pos2Meth` = "chr2:131910164",
  `Pos3Meth` = "chr2:131910180",
  `Pos4Meth` = "chr2:131910201"
)
prnpPyroPosdf <- reshape2::melt(prnpPyroPos, value.name = "Position") %>% tibble::rownames_to_column(., var = "prnpPosition")

pPrnPrnpPyro <- input_pyro$prnp %>% 
  # keep only the samples that were in the final dataset
  filter(SAMPLE %in% row.names(phenotype$rTg4510)) %>% 
  mutate(Sample.group = Group.ID) %>% dplyr::select(Age, Sample.group, contains("Pos")) %>% 
  reshape2::melt(id = c("Age","Sample.group"), variable.name = "Position", value.name = "methylation") %>% 
  mutate(Sample.group = factor(Sample.group, levels = c("WT","TG"))) %>%
  ggplot(., aes(x = Sample.group, y = methylation, fill = Sample.group)) + geom_boxplot(aes(fill = Sample.group), outlier.shape = NA) + 
  facet_grid(~Position, labeller = as_labeller(prnpPyroPos)) +
  geom_jitter(aes(colour = Sample.group),width = 0.25, size = 2) +
  #geom_point(aes(fill = Sample.group), size = 2, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values = c(alpha("black",0.2), color_Tg4510_TG),guide="none") +
  scale_colour_manual(values = c("black", color_Tg4510_TG),guide="none") +
  labs(x = "Genotype", y = "Methylation (%)") + theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "grey", linetype = "dotted"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) 



plot_pyro_rrbs_corr(input_pyro$prnp, prnpPyroPosdf, rTg4510_rrbs_beta, phenotype$rTg4510)

# plot DMP and Track for Ank1
pAnk1DMP <- plotGeneTrackDMP(sigResults=sigRes$rTg4510$Genotype, betaMatrix=sigBeta$rTg4510$Genotype, phenotypeFile=phenotype$rTg4510, 
                             gene="Ank1", transcript="ENSMUST00000110688.8", colour = "rTg4510", boxplot = TRUE)

tAnk1DMP <- plot_DMP(betaMatrix=sigBeta$rTg4510$Genotype, phenotypeFile=phenotype$rTg4510, 
         position = c("chr8:23023240","chr8:23023210","chr8:23023192"), table = TRUE) %>% 
  mutate(method = "RRBS")


ank1PyroPos <- c(
  `Pos1Meth` = "chr8:23023192",
  #`Pos2Meth` = "chr8:23023210",
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
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "grey", linetype = "dotted"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) 


pAnk1Pyro <- input_pyro$ank1 %>% dplyr::select(Age, Group.ID, Pos1Meth, Pos3Meth) %>% 
  reshape2::melt(id = c("Age","Group.ID"), variable.name = "Position", value.name = "methylation") %>% 
  mutate(Age = as.factor(stringr::str_remove(Age,"m"))) %>% 
  mutate(Group.ID = factor(Group.ID, levels = c("WT","TG"))) %>%
  ggplot(., aes(x = Group.ID, y = methylation, fill = Group.ID)) + geom_boxplot(aes(fill = Group.ID), outlier.shape = NA) + 
  facet_grid(~Position, labeller = as_labeller(ank1PyroPos)) +
  #geom_point(aes(fill = Sample.group), size = 2, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values = c(alpha("black",0.2), color_Tg4510_TG),guide="none") +
  scale_colour_manual(values = c("black", color_Tg4510_TG),guide="none") +
  labs(x = "Genotype", y = "Methylation (%)") + theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "grey", linetype = "dotted"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) 


plot_pyro_rrbs_corr(input_pyro$ank1, ank1PyroPosdf, rTg4510_rrbs_beta, phenotype$rTg4510)


# common DMPs in rTg4510 and J20 pathology
intersect(sigRes$J20$Pathology$Position, sigRes$rTg4510$Pathology$Position)

plot_DMP(model="rTg4510", betaMatrix=sigBeta$rTg4510$Pathology, 
         phenotypeFile=phenotype$rTg4510, position = intersect(sigRes$J20$Pathology$Position, sigRes$rTg4510$Pathology$Position), pathology = TRUE)

plot_DMP(model="J20", betaMatrix=sigBeta$J20$Pathology, 
         phenotypeFile=phenotype$J20, position = intersect(sigRes$J20$Pathology$Position, sigRes$rTg4510$Pathology$Position), pathology = TRUE)


sigRes$rTg4510$Pathology[sigRes$rTg4510$Pathology$Position %in% intersect(sigRes$J20$Pathology$Position, sigRes$rTg4510$Pathology$Position),]
sigRes$J20$Pathology[sigRes$J20$Pathology$Position %in% intersect(sigRes$J20$Pathology$Position, sigRes$rTg4510$Pathology$Position),]


# J20 
Nutf2 <- plotGeneTrackDMP(sigRes$J20$Genotype, sigBeta$J20$Genotype, phenotype$J20, "Nutf2", "ENSMUST00000008594.8", colour = "J20", boxplot = TRUE)
Tenm2 <- plotGeneTrackDMP(sigRes$J20$Genotype, sigBeta$J20$Genotype, phenotype$J20, "Tenm2", "ENSMUST00000102801.7", colour = "J20", boxplot = TRUE)
Ncam2 <- plotGeneTrackDMP(sigRes$J20$PathologyCommonInteraction, sigBeta$J20$Pathology, phenotype$J20, "Ncam2", "ENSMUST00000037785.13", colour = "J20", pathology =  TRUE, boxplot = TRUE)
Prmt8 <- plotGeneTrackDMP(sigRes$J20$PathologyCommonInteraction, sigBeta$J20$Pathology, phenotype$J20, "Prmt8", "ENSMUST00000032500.8", colour = "J20", boxplot = TRUE, pathology = TRUE)
Zfp518b <- plotGeneTrackDMP(sigRes$J20$PathologyCommonInteraction, sigBeta$J20$Pathology, phenotype$J20, "Zfp518b", "ENSMUST00000179555.7", colour = "J20", boxplot = TRUE, pathology = TRUE)
Zmiz1 <- plotGeneTrackDMP(sigRes$J20$PathologyCommonInteraction, sigBeta$J20$Pathology, phenotype$J20, "Zmiz1", "ENSMUST00000162645.7", colour = "J20", boxplot = TRUE, pathology = TRUE)

plot_grid(Nutf2, Tenm2, scale = 0.9, labels = c("i","ii"))
plot_grid(Ncam2, Zmiz1, Zfp518b, Prmt8, scale = 0.9, labels = c("i","ii","iii","iv"))

# venn diagram of hippocampus vs entorhinal cortex
HipECXVennrTg4510 <- plot_grid(venn.diagram(
  x = list(rTg4510_array_sig$ECX$Genotype$position, rTg4510_array_sig$ECX$Pathology$position,  
           rTg4510_array_sig$HIP$Genotype$position, rTg4510_array_sig$HIP$Pathology$position),
  category.names = c("ECX_Genotype" , "ECX_Pathology", "HIP_Genotype", "HIP_Pathology"),
  fill = pastelColours,
  filename = NULL
))

HipECXVennJ20 <- plot_grid(venn.diagram(
  x = list(J20_array_sig$ECX$Genotype$position, J20_array_sig$ECX$Pathology$position,  
           J20_array_sig$HIP$Genotype$position, J20_array_sig$HIP$Pathology$position),
  category.names = c("ECX_Genotype" , "ECX_Pathology", "HIP_Genotype", "HIP_Pathology"),
  fill = pastelColours,
  filename = NULL
))
plot_grid(HipECXVennrTg4510,HipECXVennJ20, scale = 0.9, labels = c("i","ii"))

ECXrTg4510Unique2HIP <- setdiff(
  # ECX
  c(rTg4510_array_sig$ECX$Genotype$position,
                  intersect(rTg4510_array_sig$ECX$Interaction$position, rTg4510_array_sig$ECX$Pathology$position)),
  # HIP
  c(rTg4510_array_sig$HIP$Genotype$position,
          intersect(rTg4510_array_sig$HIP$Interaction$position, rTg4510_array_sig$HIP$Pathology$position)))


rbind(
  rTg4510_array_sig$ECX$Genotype[rTg4510_array_sig$ECX$Genotype$position %in% ECXrTg4510Unique2HIP,] %>% 
    arrange(FDR_adj_Genotype) %>% 
    dplyr::select(position, cpg, FDR_adj_Genotype, ChIPseeker_GeneSymbol, distanceToTSS) %>% 
    rename("FDR_adj_Genotype" = "FDR_adj") %>% mutate(test = "genotype"),
  rTg4510_array_sig$ECX$Pathology[rTg4510_array_sig$ECX$Pathology$position %in% ECXrTg4510Unique2HIP,] %>% 
    arrange(FDR_adj_Pathology) %>% 
    dplyr::select(position, cpg, FDR_adj_Pathology, ChIPseeker_GeneSymbol, distanceToTSS) %>%
    rename("FDR_adj_Pathology" = "FDR_adj") %>% mutate(test = "pathology")
)

commonECXHIPplots <- list(
  Dcaf5 = plot_DMP_byTissue(ECXbetaMatrix=sigBeta$rTg4510$Genotype, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                  ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr12:80436248"),
  Satb1  = plot_DMP_byTissue(ECXbetaMatrix=sigBeta$rTg4510$Genotype, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                            ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr17:51746925"),
  Cltc  = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                             ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr11:8670046"),
  Mapt = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                              ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr11:104318231"),
  Ncapg2  = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                            ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr12:116425797"),
  Fgf14 = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                          ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr14:124676565")
)

plot_grid(plotlist = commonECXHIPplots, labels = c("A","B","C","D","E","F"), scale = 0.9)

commonECXHIPrTg4510plots <- list(
  Dennd1a = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                              ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr2:37946161", 
                              pathology = TRUE, gene = "Dennd1a"),
  Rapgefl1 = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                              ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr11:98838683", 
                              pathology = TRUE, gene = "Rapgefl1")
)
plot_grid(plotlist = commonECXHIPrTg4510plots, labels = c("i","ii"), scale = 0.9)

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
plot_grid(plotlist = HIPrTg4510plots, labels = c("i","ii","iii","iv"), scale = 0.9)

HIPJ20plots <- list(
  Mir568  = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                          ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr16:43609394", 
                          gene = "Mir568", model = "J20"),
  Mctp1 = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                            ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr13:76810803", 
                            gene = "Mctp1", model = "J20"),
  Sox4 = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                            ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr13:28949481", 
                            gene = "Sox4", model = "J20", pathology = TRUE),
  Cetn3 = plot_DMP_byTissue(ECXbetaMatrix=rTg4510_array_beta, HIPbetaMatrix=rTg4510_array_HIP_beta, 
                           ECXphenotypeFile=phenotype$rTg4510, HIPphenotypeFile=phenotype$rTg4510_HIP, position ="chr13:81828611", pathology = TRUE, 
                           gene = "Cetn3", model = "J20")
)
plot_grid(plotlist = HIPJ20plots, labels = c("i","ii","iii","iv"), scale = 0.9)


# hippocampus plots
plot_DMP(betaMatrix=rTg4510_array_HIP_beta, phenotypeFile=phenotype$rTg4510_HIP, position = "chr4:109806599", pathology = TRUE)
plot_DMP(betaMatrix=rTg4510_array_HIP_beta, phenotypeFile=phenotype$rTg4510_HIP, position = "chr4:109806599", interaction = TRUE)
plot_DMP(betaMatrix=rTg4510_array_HIP_beta, phenotypeFile=phenotype$rTg4510_HIP, position = "chr2:37946161", pathology = TRUE)
plot_DMP(betaMatrix=rTg4510_array_HIP_beta, phenotypeFile=phenotype$rTg4510_HIP, position = "chr2:37946161", interaction = TRUE)

           humanAllGeneList),
  category.names = c("rTg4510_ECX","rTg4510_HIP", "J20_ECX", "J20_HIP","Human"),
  fill = c(label_colour("rTg4510"), alpha(label_colour("rTg4510"),0.2), label_colour("J20"),alpha(label_colour("J20"),0.2),"yellow"),
  filename = NULL
))


# human comparisons
Prdm16_J20 <- plotGeneTrackDMP(sigRes$J20$Genotype, sigBeta$J20$Genotype, phenotype$J20, "Prdm16", "ENSMUST00000030902.12", colour = "J20", boxplot = TRUE)
Prdm16_rTg4510 <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Prdm16", "ENSMUST00000030902.12", colour = "rTg4510", boxplot = TRUE)
plot_grid(Prdm16_J20, Prdm16_rTg4510, nrow = 1)
sigRes$J20$Genotype[sigRes$J20$Genotype$ChIPseeker_GeneSymbol == "Prdm16",]
sigRes$rTg4510$Genotype[sigRes$rTg4510$Genotype$ChIPseeker_GeneSymbol == "Prdm16",]

plot_DMP(sigBeta$rTg4510$Genotype, phenotype$rTg4510, position = c("chr4:154346846"))

dat <- humanDMPres %>% mutate(humanPosition = paste0("chr",CHR,":", BP)) %>% 
  dplyr::select(humanPosition,Effect_fixed....,P_fixed,UCSC.Nearest.Gene) %>%  
  full_join(sigRes$rTg4510_Human[,c("Position","FDR_adj_Genotype","delta","human")], 
            ., by = c("human" = "UCSC.Nearest.Gene"), relationship = "many-to-many") %>%
  `colnames<-`(c("mousePosition","mouseFDRGenotype","mouseESGenotype","humanGene","humanPosition","humanES","humanPvalue")) %>%
  dplyr::filter(mousePosition != "NA") %>%
  mutate(labels = paste0(mousePosition,",", humanPosition))
length(unique(dat$mousePosition))
length(unique(dat$humanPosition))
ggplot(dat, aes(x = mouseESGenotype, y = humanES, colour = humanGene, label = humanGene)) + 
  geom_point(size = 3) + 
  geom_text_repel(show.legend = FALSE) +
  geom_hline(yintercept=0, linetype = "dotted") +
  geom_vline(xintercept=0, linetype = "dotted") +
  theme_classic() +
  labs(x = "Effect size in rTg4510 TG vs WT", y = "Effect size in human AD studies") +
  scale_colour_discrete(name = "Gene")

dat <- humanDMPres %>% mutate(humanPosition = paste0("chr",CHR,":", BP)) %>% 
  dplyr::select(humanPosition,Effect_fixed....,P_fixed,UCSC.Nearest.Gene) %>%  
  full_join(sigRes$J20_Human[,c("Position","FDR_adj_Genotype","delta","human")], 
            ., by = c("human" = "UCSC.Nearest.Gene"), relationship = "many-to-many") %>%
  `colnames<-`(c("mousePosition","mouseFDRGenotype","mouseESGenotype","humanGene","humanPosition","humanES","humanPvalue")) %>%
  dplyr::filter(mousePosition != "NA") %>%
  mutate(labels = paste0(mousePosition,",", humanPosition))
length(unique(dat$mousePosition))
length(unique(dat$humanPosition))
sigRes$J20$Genotype[sigRes$J20$Genotype$ChIPseeker_GeneSymbol == "Tspan14",]
View(dat)
ggplot(dat, aes(x = mouseESGenotype, y = humanES, colour = humanGene, label = humanGene)) + 
  geom_point(size = 3) + 
  geom_text_repel(show.legend = FALSE) +
  geom_hline(yintercept=0, linetype = "dotted") +
  geom_vline(xintercept=0, linetype = "dotted") +
  theme_classic() +
  labs(x = "Effect size in rTg4510 TG vs WT", y = "Effect size in human AD studies") +
  scale_colour_discrete(name = "Gene")

Tspan14 <- plotGeneTrackDMP(sigResults=sigRes$J20$Genotype, betaMatrix=sigBeta$J20$Genotype, phenotypeFile=phenotype$J20, 
                 gene="Tspan14", transcript="ENSMUST00000047652.5", colour = "J20", boxplot = TRUE,
                 position = "chr14:40966816")

#-------------- Output -------------

plot_grid(pCluster$rTg4510Genotype$gtable, pCluster$rTg4510Pathology$gtable)
plot_grid(Mapt, Prnp, Ncapg2, Fgf14, labels = c("A","B","C","D"), scale = 0.9)
plot_grid(Dcaf5, Arsi, Creb3l4, As3mt, labels = c("i","ii","iii","iv"), scale = 0.9)
plot_grid(Zfp423, Adk, Insyn2b, Ankrd52, labels = c("i","ii","iii","iv"), scale = 0.9)

plot_grid(pPrnPrnPDMP,pPrnPrnpPyro, rel_heights = c(0.6,0.4), labels = c("A","B"))
plot_grid(pAnk1DMP,pAnk1Pyro, rel_heights = c(0.6,0.4), labels = c("A","B"))
plot_grid(
  plot_DMP(model="rTg4510", betaMatrix=sigBeta$rTg4510$Genotype, phenotypeFile=phenotype$rTg4510, sig = sigRes$rTg4510$Genotype),
  plot_DMP(model="rTg4510", betaMatrix=sigBeta$rTg4510$Pathology, phenotypeFile=phenotype$rTg4510, sig=sigRes$rTg4510$Pathology, pathology = TRUE),
  nrow=2
)

plot_grid(
  plot_DMP(model="J20", betaMatrix=sigBeta$J20$Genotype, phenotypeFile=phenotype$J20, sig = sigRes$J20$Genotype), 
  plot_DMP(model="J20", betaMatrix=sigBeta$J20$Pathology, phenotypeFile=phenotype$J20, sig = sigRes$J20$Pathology, pathology = TRUE), 
  nrow=2
)

pdf(paste0(dirnames$paper,"/resFigures/TopResults_GenotypePathology.pdf"), width = 15, height = 8)
plot_DMP(model="rTg4510", betaMatrix=sigBeta$rTg4510$Genotype, phenotypeFile=phenotype$rTg4510, sig = sigRes$rTg4510$Genotype) 
plot_DMP(model="rTg4510", betaMatrix=sigBeta$rTg4510$Pathology, phenotypeFile=phenotype$rTg4510, sig=sigRes$rTg4510$Pathology, pathology = TRUE) 
plot_DMP(model="J20", betaMatrix=sigBeta$J20$Genotype, phenotypeFile=phenotype$J20, sig = sigRes$J20$Genotype) 
plot_DMP(model="J20", betaMatrix=sigBeta$J20$Pathology, phenotypeFile=phenotype$J20, sig = sigRes$J20$Pathology, pathology = TRUE) 
dev.off()