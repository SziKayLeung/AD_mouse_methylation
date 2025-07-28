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

#-------------- Supplementary Figures -------------

## --- Figure 1: plot RRBS CpG sites by annotation

plot_annotate_sites()
plot_grid(pCluster$rTg4510Genotype$gtable, pCluster$rTg4510Pathology$gtable)


## --- Figure 2: correlation of probes for RRBS vs Array

# rTg4510
rTg4510_common_probes <- commonStatsDescription(rTg4510_rrbs_beta, rTg4510_array_beta)
rTg4510_corr_arrayRRBS <- corrPlotCommonProbes(rTg4510_rrbs_beta, rTg4510_array_beta, rTg4510_common_probes)
# J20
J20_common_probes <- commonStatsDescription(J20_rrbs_beta, J20_array_beta)
J20_corr_arrayRRBS <- corrPlotCommonProbes(J20_rrbs_beta, J20_array_beta, J20_common_probes)
plot_grid(rTg4510_corr_arrayRRBS, J20_corr_arrayRRBS)


## --- Figure 3: Top-ranked DMPs in rTg4510 genotype due to transgene

Mapt <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Mapt", "ENSMUST00000100347.10", boxplot = TRUE, colour = "rTg4510")
Prnp <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, c("Prn","Prnp"), "ENSMUST00000091288.12", colour = "rTg4510")
Fgf14 <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Fgf14", "ENSMUST00000095529.9", boxplot = TRUE, colour = "rTg4510")
Ncapg2 <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Ncapg2", "ENSMUST00000084828.4", boxplot = TRUE, colour = "rTg4510")


## --- Figure 4: Pyrosequencing of Prnp

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


## --- Figure 5: Effect size of ECX vs HIP in rTg4510

p5 <- effectSizeComparisons(rTg4510_array_results$Genotype, rTg4510_HIP_array_results$Genotype, "Array", "Genotype", "HIP", animal="rTg4510")
p6 <- effectSizeComparisons(rTg4510_array_results$Pathology, rTg4510_HIP_array_results$Pathology, "Array", "Pathology", "HIP", animal="rTg4510")
p7 <- effectSizeComparisons(rTg4510_array_sig$ECX$Genotype, rTg4510_array_sig$HIP$Genotype, "Array", "Genotype", "HIP", animal = "rTg4510")
p8 <- effectSizeComparisons(rTg4510_array_sig$ECX$Pathology, rTg4510_array_sig$HIP$Pathology, "Array", "Pathology", "HIP", animal = "rTg4510")
pdf(paste0(output,"Figures/rTg4510ECXvsHIP.pdf"),  width = 10, height = 15)
plot_grid(p5,p6,p7,p8,labels = c("A","B","C","D"))
dev.off()


## --- Figure 6: Common DMPs between rTg4510 ECX and HIP

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

## --- Figure 7: Effect size of ECX vs HIP in J20

p9 <- effectSizeComparisons(J20_array_results$Genotype, J20_HIP_array_results$Genotype, "Array", "Genotype", "HIP", animal="J20")
p10 <- effectSizeComparisons(J20_array_results$Pathology, J20_HIP_array_results$Pathology, "Array", "Pathology", "HIP", animal="J20")
p11 <- effectSizeComparisons(J20_array_sig$ECX$Genotype, J20_array_sig$HIP$Genotype, "Array", "Genotype", "HIP", animal = "J20")
p12 <- effectSizeComparisons(J20_array_sig$ECX$Pathology, J20_array_sig$HIP$Pathology, "Array", "Pathology", "HIP", animal = "J20")
pdf(paste0(output,"Figures/J20ECXvsHIP.pdf"),  width = 10, height = 15)
plot_grid(p9,p10,p11,p12,labels = c("A","B","C","D"))
dev.off()


## --- Figure 8: Epigenetic clock

load(file = paste0(zenDir, "/5_epigeneticClock/rTg4510Clock.RData"))
load(file = paste0(zenDir, "/5_epigeneticClock/J20Clock.RData"))

plot_grid(
  plot_clock(rTg4510Clocks$ECX,"ECX",model="rTg4510", boxplot = FALSE),
  plot_clock(rTg4510Clocks$HIP,"HIP",model="rTg4510", boxplot = FALSE),
  plot_clock(J20Clocks$ECX,"ECX",model="J20", boxplot = FALSE),
  plot_clock(J20Clocks$HIP,"HIP",model="J20", boxplot = FALSE),
  labels = c("A","B","C","D")
)

## --- Figure 9: Prdm16

Prdm16_J20 <- plotGeneTrackDMP(sigRes$J20$Genotype, sigBeta$J20$Genotype, phenotype$J20, "Prdm16", "ENSMUST00000030902.12", colour = "J20", boxplot = TRUE)
Prdm16_rTg4510 <- plotGeneTrackDMP(sigRes$rTg4510$Genotype, sigBeta$rTg4510$Genotype, phenotype$rTg4510, "Prdm16", "ENSMUST00000030902.12", colour = "rTg4510", boxplot = TRUE)
plot_grid(Prdm16_rTg4510,Prdm16_J20, labels = c("i","ii"), nrow = 1, scale = 0.95)
plot_DMP(sigBeta$rTg4510$Genotype, phenotype$rTg4510, position = c("chr4:154346846"))