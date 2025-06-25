#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Code to generate main figure (5) - Ank1 (human translational aspect) 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## p1: genotype WT vs TG of Ank1 DMPs - boxplots
## p2: genotype*interaction of Ank1 DMPs - line
## p3: pathology of 9 ANK1 DMPs (Scatter across pathology variable)
## p4: pyrosequencing validation (Rhian experiment)
## p5: track of representative Ank1 transcript & DMPs 
## 
## supplementary figure of Ank1 pyrosequencing vs rrbs: 1_correlate_rrbs_pyro.R output
## Genotype probes: chr8:23023193, chr8:23023211, chr8:23023240, chr8:23023241
## Interaction probes: chr8:23023193, chr8:23023211, chr8:23055131
## Pathology probes: chr8:23058532, chr8:23058569, chr8:23058578, chr8:23058580, chr8:23058583, chr8:23063501, chr8:23063549, chr8:23063550, chr8:23063585


## ------------ packages ------------

suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("gridExtra"))
suppressMessages(library("grid"))


## ------------ input ------------

rootDir = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/"
scriptDir = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config.R"))

## Results after applying BiSeq and smoothing
# load from two mouse models and save into list
RRBS_complete <- list()
load(file = paste0(rootDir, "1_processed/rrbs/rTg4510_RRBSbetasComplete.RData"))
RRBS_complete$rTg4510 <- RRBS_completebetas
load(file = paste0(rootDir, "1_processed/rrbs/J20_RRBSbetasComplete.RData"))
RRBS_complete$J20 <- RRBS_completebetas

# phenotype data
pheno <- list(
  rTg4510 = read.csv(paste0(rootDir, "0_metadata/Tg4510_coldata_RRBS.csv"), stringsAsFactors=FALSE),
  J20 = read.csv(paste0(rootDir, "0_metadata/J20_coldata_RRBS.csv"), stringsAsFactors=FALSE)
)
pheno <- lapply(pheno, function(x) x %>% mutate(Age_months = factor(Age_months), Genotype = factor(Genotype, levels = c("WT","TG"))))

# rTg4510 annotated pathology results 
rTg4510_rrbs_pathology = read.csv(paste0(dirnames$output,"rTg4510_rrbs_pathology.csv"))
Ank1PathologyDMPs <- rTg4510_rrbs_pathology[rTg4510_rrbs_pathology$ChIPseeker_GeneSymbol == "Ank1",c("Position")]

# pyrosequencing data: Ank1
source(paste0(scriptDir, "8_Pyro_assays/pyro.config.R"))

# Ank1 reference gtf 
# subsetted Ank1 annotations from M22 reference GENCODE gtf
Ank1refGtf <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/references/annotation/MouseAnk1_gencodeM22.gtf"
Ank1refGtf <- as.data.frame(rtracklayer::import(Ank1refGtf))

# Table of Ank2 DMPs manually created and uploaded
Ank1DMP <- read.csv(paste0(scriptDir, "0_PaperFigs/otherInput/Ank1DmpPositions.csv")) %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))
Ank1DMP <- distinct(Ank1DMP)

# rTg4510 colour "blue"
color_Tg4510_TG <- "#00AEC9"


## ------------ plot functions ------------

# create_df_plot: generate dataframe for downstream plotting
# datawrangle with annotated rrbs sites and merge with phenotype data
# params:
  # x = list of DMPs <chr8:XXXXXX>
create_df_plot <- function(x){
  
  dat <- RRBS_complete$rTg4510 %>% tibble::rownames_to_column(var = "site") %>% 
    filter(site %in% x) %>% 
    reshape2::melt(variable.name = "Sample_ID", value.name = "methylation", id = "site") %>%
    left_join(., pheno$rTg4510[,c("Sample_ID","Genotype","Age_months","ECX")], by = "Sample_ID") %>% 
    mutate(methylation = methylation * 100, Genotype = factor(Genotype, levels = c("WT","TG")))
  
  return(dat)
}


# ggtranscript_plot: plot track
# params:
  # gexons = df of gtf with exons subsetted 
  # intron = <TRUE/FALSE> to plot line between exons (i.e. FALSE for DMP plotting)
ggtranscript_plot <- function(gexons, intron=TRUE){
  g <- ggplot(gexons, aes(xstart = start, xend = end, y = transcript_id)) +
    geom_range() +
    labs(y = NULL) +  
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.background = element_rect(fill='white'),
          plot.margin = unit(c(0, 0, 0, 0), "cm")) 
  if(isTRUE(intron)){
    g <- g + geom_intron(data = to_intron(gexons, "transcript_id"), aes(strand = strand))
  }
  return(g)
}
  

## ------------ p1: heatmap ------------

source(paste0(scriptDir, "4_HumanGeneListComparisons/4_compGenesSim_txtFiles.R"))
p1 <- compGenes(dirnames$humanAnnot)


## ------------ p2: genotype ------------

p2 <- create_df_plot(c("chr8:23023193","chr8:23023211","chr8:23023240","chr8:23023241")) %>%
  ggplot(., aes(x = Genotype, y = methylation, fill = Genotype)) + geom_boxplot(outlier.shape = NA) + 
  #geom_point(aes(colour = Genotype))  +
  facet_grid(~site) + 
  #scale_fill_manual(values = c(alpha("black",0.2), alpha(color_Tg4510_TG,0.2)),guide="none") +
  scale_fill_manual(values = c(alpha("black",0.2), color_Tg4510_TG),guide="none") +
  scale_colour_manual(values = c("black", color_Tg4510_TG),guide="none") +
  labs(x = "Genotype", y = "Methylation (%)") +
  theme_classic() + 
  theme(panel.border = element_rect(fill = NA, color = "grey", linetype = "dotted"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) 
  

## ------------ p3: interaction - genotype*age ------------

interaction <- create_df_plot(c("chr8:23023193","chr8:23023211","chr8:23055131"))
p3 <- ggplot(interaction, aes(x = Age_months, y = methylation)) + 
  geom_point(aes(colour = Genotype)) + 
  facet_grid(~site) +
  scale_fill_manual(values = c(alpha("black",0.2), color_Tg4510_TG),guide="none") +
  stat_summary(data=interaction, aes(x=Age_months, y=methylation, group=Genotype), fun ="mean", geom="line", linetype = "dotted") +
  labs(y = "Methylation (%)", x = "Age (months)") + theme_classic() + 
  scale_colour_manual(values = c("black", color_Tg4510_TG),guide="none") +
  theme(panel.border = element_rect(fill = NA, color = "grey", linetype = "dotted"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) 


## ------------ p4: pathology ------------

pathology <- create_df_plot(Ank1PathologyDMPs)

## scatter plot of just TG mice and colour coded by site
#p3 <- pathology %>% filter(Genotype == "TG") %>% 
#  ggplot(., aes(x = ECX, y = methylation, colour = site)) + geom_point(size = 2) +
#  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, linetype = "dotted") +
#  #stat_summary(data=dat, aes(x=as.numeric(ECX), y=methylation), fun ="mean", geom="line", linetype = "dotted")+
#  labs(x = "Pathology", y = "Methylation (%)", x = "Age (months)") +
#  theme(panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank()) +
#  theme_classic()

# separate plot for each DMP across pathology
pathologyPlots <-list()
for(n in 1:length(unique(pathology$site))){
  x <- unique(pathology$site)[[n]]
  pathologyPlots[[n]] <- 
    pathology %>% filter(site == x) %>%
    ggplot(., aes(x = ECX, y = methylation, colour = Genotype)) + geom_point(size = 2) +
    geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, linetype = "dotted") +
    labs(x = NULL, y = NULL, subtitle = x) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme_classic() +
    scale_colour_manual(values = c("black", color_Tg4510_TG),guide="none")
  
}
pathologyListPlots <- plot_grid(plotlist = pathologyPlots)
y.grob <- textGrob("Metylation (%)", gp=gpar(fontsize=12), rot=90)
x.grob <- textGrob("Pathology", gp=gpar(fontsize=12))
p4  <- grid.arrange(arrangeGrob(pathologyListPlots, left = y.grob, bottom = x.grob))


## ------------ p5: pyrosequencing validation ------------

ank1PyroPos <- c(
  `Pos1Meth` = "chr8:23023193",
  `Pos2Meth` = "chr8:23023211",
  `Pos3Meth` = "chr8:23023241"
)

p5 <- input_pyro$ank1 %>% select(Age, Sample.group, contains("Pos")) %>% 
  reshape2::melt(id = c("Age","Sample.group"), variable.name = "Position", value.name = "methylation") %>% 
  mutate(Age = as.factor(stringr::str_remove(Age,"m"))) %>% 
  mutate(Sample.group = factor(Sample.group, levels = c("WT","TG"))) %>%
  ggplot(., aes(x = Age, y = methylation, fill = Sample.group)) + geom_boxplot(aes(fill = Sample.group), outlier.shape = NA) + 
  facet_grid(~Position, labeller = as_labeller(ank1PyroPos)) +
  #geom_point(aes(fill = Sample.group), size = 2, shape = 21, position = position_jitterdodge()) +
  scale_fill_manual(values = c(alpha("black",0.2), color_Tg4510_TG),guide="none") +
  scale_colour_manual(values = c("black", color_Tg4510_TG),guide="none") +
  labs(x = "Genotype", y = "Methylation (%)") + theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "grey", linetype = "dotted"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) 


## ------------ p6: Ank1 tracks ------------

# using ggtranscript for visaulilsation
# ENSMUST00000110688.8 representative ANK1 transcript (longest)
gexons <- Ank1refGtf %>% dplyr::filter(type == "exon") %>% select(transcript_id, start, end, seqnames, strand)
gexons <- rbind(gexons,Ank1DMP)

# subset to transcripts and DMP of interest
gexonsT <- gexons %>% filter(transcript_id %in% c("ENSMUST00000110688.8")) 
gexonsD <- gexons %>% filter(transcript_id %in% c("DMP")) 

# plot tracks 
g1 <- ggtranscript_plot(gexonsT)
g2 <- ggtranscript_plot(gexonsD,intron=FALSE)
p6 <- plot_grid(g2,g1,ncol=1,align="hv",rel_heights = c(1, 1))


## ------------ Concatenated plot ------------

# plot A: Isabel heatmap
allP <- plot_grid(
  plot_grid(p1,p2,p3,p4, labels = c("B","C","D","E")),
  p5,
  rel_heights = c(0.8,0.2),nrow=2, labels=c("","F"))


pdf(paste0(dirnames$output,"/MainFigures5_Ank1.pdf"), width = 16, height = 10)
allP
dev.off()