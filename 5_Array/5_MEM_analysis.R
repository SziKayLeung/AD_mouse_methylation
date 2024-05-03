#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Downstream analysis of mixed effects model to plot top-ranked differential methylation differences between ECX and HIP
## Post running 4_MEM_J20_Tissue_sbatch.sh & 4_MEM_rTg4510_Tissue_sbatch.sh
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------
## 
## Motivation/Order
## 1. ADahir:performed array QC and summary 
##    PrepMixedEffectsStats() adapted for output of beta values and phenotype files from rTg4510 and J20 
## 2. EWalker performed mixed effects model regression to identify differential methylation differences between ECX and HIP from array data
##    ~ Tissue + Tissue*Genotype + Genotype + Age_months + (1|MouseID)
##    Model provides two p-values: PrZ.GenotypeTG, PrZ.Tissue.GenotypeTG
## 3. EWalker performed mixed effects model regression to identify differential methylation differences between genotype from array data in ECX
##    ~ Genotype + Age_months + Genotype*Age_months + Chip_ID
##    Model provides two p-values: PrZ.GenotypeTG, PrZ.GenotypeTG.Age
## 4. SLeung to annotate array CpG sites using chipseeker (for consistency with RRBS comparison)
##    FDR correction of two p-values: FDR_adj_GenotypeTG, FDR_adj_TissueGenotype
##    FDR correction of two p-values: FDR_adj_GenotypeTG, FDR_adj_GenotypeAge
## 5. SLeung to identify top ranked differential methylation differences between tissue and plot
##     Tissue difference, but not in Genotype (FDR_adj_TissueGenotype < 0.05)
##     Genotype difference, but consistent in tissue (FDR_adj_GenotypeTG < 0.05)
##     statistically significant: TissueGenotype and Genotype (FDR_adj_TissueGenotype < 0.05 & FDR_adj_GenotypeTG < 0.05)
## 6. SLeung to identify top ranked differential methylation difference between genotype and plot
## 7. SLeung to identify top ranked differential methylation difference between genotype*age (progressive) and plot




#-------------- packages -------------

suppressMessages(library(stringr))
suppressMessages(library(cowplot))
suppressMessages(library(GenomicRanges)) # GenomicRanges_1.38.0
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ChIPseeker))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


#-------------- functions -------------

# chipseeker annotation
source("/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/4_DMP_annotation/0_source_functions.R")
# data wrangle for beta output, phenotype file output, and mixedeffectsmodels results
source("/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/5_Array/5a_MEM_analysisFunctions.R")


#-------------- input -------------

# directory for rTg4510 and J20 RRBS and Genotype 
dirnames <- list(
  AishaArray = "/lustre/projects/Research_Project-191406/Aisha/data/Array/",
  EmmaArray = "/lustre/projects/Research_Project-191406/EmmaW/Array/Results/",
  output = "/lustre/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/0_ZenOutput/3_array/"
)

# input files
input <- list(
  mm10_Manifest = read.csv(paste0(dirnames$AishaArray, "mm10_Manifest.csv"), stringsAsFactors = F, header = T, row.names = 1),
  # array raw values
  J20Array= read.csv(paste0(dirnames$AishaArray,"J20_coldata_VertebrateArray.csv"), header = T, stringsAsFactors = F),
  rTg4510Array = read.csv(paste0(dirnames$AishaArray,"Tg4510_coldata_VertebrateArray.csv"),header = T, stringsAsFactors = F),
  # mixed effects results = tissue
  # results generated from 4a_MixedModelsEffectBetaReg_rTg4510_Tissue_EW.R
  rTg4510_Tissue = read.table(paste0(dirnames$EmmaArray, "rTg4510/MixedModelResults_Tissue_rTg4510.csv"), header = T, sep = ",", row.names = "X"),
  J20_Tissue = read.table(paste0(dirnames$EmmaArray, "J20/MixedModelResults_Tissue_J20.csv"), header = T, sep = ",", row.names = "X"),
  # mixed effects results = genotype ECX
  # results generated from 4b_MixedModelsEffectBetaReg_J20_ECX_EW.R
  rtg4510_Geno = read.table(paste0(dirnames$EmmaArray, "rTg4510/MixedModelResultsECX_rTg4510.csv"), header = T, sep = ",", row.names = "X"),
  J20_Geno = read.table(paste0(dirnames$EmmaArray, "J20/MixedModelResultsECX_J20.csv"), header = T, sep = ",", row.names = "X")
)


#-------------- prepare input files -----------

stats <- list(
  rTg4510_tissue = PrepMixedEffectsStats(input$rTg4510_Tissue, input$mm10_Manifest, intTerm="Tissue.Genotype"),
  J20_tissue = PrepMixedEffectsStats(input$J20_Tissue, input$mm10_Manifest, intTerm="Tissue.Genotype"),
  rTg4510_ECX_Geno = PrepMixedEffectsStats(input$rtg4510_Geno, input$mm10_Manifest, intTerm="Geno.Age"),
  J20_ECX_Geno = PrepMixedEffectsStats(input$J20_Geno, input$mm10_Manifest, intTerm="Geno.Age")
)

# metadata of the combined phenotype and beta values
meta <- list(
  rTg4510_tissue = PrepBetasPhenoFiles(input$rTg4510Array, "rTg4510"),
  rTg4510_ECX_Geno = PrepBetasPhenoFiles(input$rTg4510Array, "rTg4510"),
  J20_tissue = PrepBetasPhenoFiles(input$J20Array, "J20"),
  J20_ECX_Geno = PrepBetasPhenoFiles(input$J20Array, "J20")
)


#-------------- significant top-ranked CpGs ------------

# top-ranked from tissue model 
topRanked_Tissue <- list(
  # Tissue difference, but not in Genotype
  rTg4510TissueOnly = stats$rTg4510_tissue %>% filter(FDR_adj_TissueGenotype < 0.05, FDR_adj_GenotypeTG > 0.05) %>% arrange(FDR_adj_TissueGenotype),
  J20TissueOnly = stats$J20_tissue %>% filter(FDR_adj_TissueGenotype < 0.05, FDR_adj_GenotypeTG > 0.05) %>% arrange(FDR_adj_TissueGenotype),
  
  # Genotype difference, but consistent in tissue
  rTg4510GenoOnly = stats$rTg4510_tissue %>% filter(FDR_adj_TissueGenotype > 0.05, FDR_adj_GenotypeTG < 0.05) %>% arrange(FDR_adj_TissueGenotype),
  J20GenoOnly = stats$J20_tissue %>% filter(FDR_adj_TissueGenotype > 0.05, FDR_adj_GenotypeTG < 0.05) %>% arrange(FDR_adj_TissueGenotype),
  
  # statistically significant: TissueGenotype and Genotype
  rTg4510Geno = stats$rTg4510_tissue %>% filter(FDR_adj_TissueGenotype < 0.05, FDR_adj_GenotypeTG < 0.05) %>% arrange(FDR_adj_GenotypeTG),
  rTg4510TissueGeno = stats$rTg4510_tissue %>% filter(FDR_adj_TissueGenotype < 0.05, FDR_adj_GenotypeTG < 0.05) %>% arrange(FDR_adj_TissueGenotype),
  J20Geno = stats$J20_tissue %>% filter(FDR_adj_TissueGenotype < 0.05, FDR_adj_GenotypeTG < 0.05) %>% arrange(FDR_adj_GenotypeTG),
  J20TissueGeno = stats$J20_tissue %>% filter(FDR_adj_TissueGenotype < 0.05, FDR_adj_GenotypeTG < 0.05) %>% arrange(FDR_adj_TissueGenotype)
)

# top-ranked from genotype model in ECX
topRanked_Genotype <- list(
  rTg4510 = stats$rTg4510_ECX_Geno %>% arrange(FDR_adj_GenotypeTG),
  J20 = stats$J20_ECX_Geno %>% arrange(FDR_adj_GenotypeTG)
)

# top-ranked from genotype and age model in ECX
topRanked_GenotypeAge <- list(
  rTg4510 = stats$rTg4510_ECX_Geno %>% arrange(FDR_adj_GenotypeAge),
  J20 = stats$J20_ECX_Geno %>% arrange(FDR_adj_GenotypeAge)
)


#-------------- plot ----------------

# wrapperplotArrayDMP  
# Aim: to input paramaters for the plotArrayDMP
# Input:
  # cpg = str: cpg value
  # model = str: name of the stats output to plot from
  # column = str: column name of the fdr value to report
  # animal = str: <rTg4510, J20> for colouring purposes
  # feature = str: <Tissue,GenotypeECX> to dictate appearnce of plot (facet to include hippocampus or to filter for ECX data)
# Output:
  # box-plot of the methylation values across age

wrapperplotArrayDMP <- function(cpg, model, column, animal, feature){
  
  p <- plotArrayDMP(cpg = cpg, betaResults = stats[[model]], betas = meta[[model]]$betas, pheno = meta[[model]]$pheno, 
                    colour = label_colour(animal), column=column, feature=feature)
  
  return(p)
}

# plots for tissue model
topRankedPlots_Tissue <- list(
  # Tissue
  rTg4510TiOnly = lapply(topRanked_Tissue$rTg4510TissueOnly$cpg[1:10], 
                         function(x) wrapperplotArrayDMP(cpg=x,model="rTg4510_tissue",column="FDR_adj_TissueGenotype",animal="rTg4510",feature="Tissue")),
  J20TiOnly = lapply(topRanked_Tissue$J20TissueOnly$cpg[1:10], 
                     function(x) wrapperplotArrayDMP(cpg=x,model="J20_tissue",column="FDR_adj_TissueGenotype",animal="J20",feature="Tissue")),
  # Genotype
  rTg4510GeOnly = lapply(topRanked_Tissue$rTg4510GenoOnly$cpg[1:10], 
                         function(x) wrapperplotArrayDMP(cpg=x,model="rTg4510_tissue",column="FDR_adj_GenotypeTG",animal="rTg4510",feature="Tissue")),
  J20GeOnly = lapply(topRanked_Tissue$J20GenoOnly$cpg[1:10], 
                     function(x) wrapperplotArrayDMP(cpg=x,model="J20_tissue",column="FDR_adj_GenotypeTG",animal="J20",feature="Tissue")),
  #Both
  rTg4510Geno = lapply(topRanked_Tissue$rTg4510Geno$cpg[1:10], 
                       function(x) wrapperplotArrayDMP(cpg=x,model="rTg4510_tissue",column="FDR_adj_GenotypeTG",animal="rTg4510",feature="Tissue")),
  rTg4510TissueGeno = lapply(topRanked_Tissue$rTg4510TissueGeno$cpg[1:10], 
                             function(x) wrapperplotArrayDMP(cpg=x,model="rTg4510_tissue",column="FDR_adj_TissueGenotype",animal="rTg4510",feature="Tissue")),
  J20Geno = lapply(topRanked_Tissue$J20Geno$cpg[1:10], 
                   function(x) wrapperplotArrayDMP(cpg=x,model="J20_tissue",column="FDR_adj_GenotypeTG",animal="J20",feature="Tissue")),
  J20TissueGeno = lapply(topRanked_Tissue$J20TissueGeno$cpg[1:10], 
                         function(x) wrapperplotArrayDMP(cpg=x,model="J20_tissue",column="FDR_adj_TissueGenotype",animal="J20",feature="Tissue"))
)

# plots for genotype model
topRankedPlots_Genotype <- list(
  rTg4510 = lapply(topRanked_Genotype$rTg4510$cpg[1:10], function(x) 
    wrapperplotArrayDMP(cpg=x,model="rTg4510_ECX_Geno",column="FDR_adj_GenotypeTG",animal="rTg4510",feature="GenotypeECX")),
  J20 = lapply(topRanked_Genotype$J20$cpg[1:10], function(x) 
    wrapperplotArrayDMP(cpg=x,model="J20_ECX_Geno",column="FDR_adj_GenotypeTG",animal="J20",feature="GenotypeECX"))
)

# plots for genotype and age model (progresssive)
topRankedPlots_GenotypeAge <- list(
  rTg4510 = lapply(stats$rTg4510_ECX_Geno$cpg[1:10], function(x) 
    wrapperplotArrayDMP(cpg=x,model="rTg4510_ECX_Geno",column="FDR_adj_GenotypeAge",animal="rTg4510",feature="GenotypeECX")),
  J20 = lapply(stats$J20_ECX_Geno$cpg[1:10], function(x) 
    wrapperplotArrayDMP(cpg=x,model="J20_ECX_Geno",column="FDR_adj_GenotypeAge",animal="J20",feature="GenotypeECX"))
)


#-------------- output -----------------

# pdf
pdf(paste0(dirnames$output,"/MixedEffectsTopRankedTissue.pdf"), width = 14, height = 8)
plot_grid(plotlist = topRankedPlots_Tissue$rTg4510TiOnly, nrow=2)
plot_grid(plotlist = topRankedPlots_Tissue$rTg4510GeOnly, nrow=2)
plot_grid(plotlist = topRankedPlots_Tissue$rTg4510Geno, nrow=2)
plot_grid(plotlist = topRankedPlots_Tissue$rTg4510TissueGeno, nrow=2)
plot_grid(plotlist = topRankedPlots_Tissue$J20TiOnly, nrow=2)
plot_grid(plotlist = topRankedPlots_Tissue$J20GeOnly, nrow=2)
plot_grid(plotlist = topRankedPlots_Tissue$J20Geno, nrow=2)
plot_grid(plotlist = topRankedPlots_Tissue$J20TissueGeno, nrow=2)
dev.off()

pdf(paste0(dirnames$output,"/MixedEffectsTopRankedGenotypeAge.pdf"), width = 14, height = 8)
plot_grid(plotlist = topRankedPlots_GenotypeAge$rTg4510, nrow=2)
plot_grid(plotlist = topRankedPlots_GenotypeAge$J20, nrow=2)
dev.off()

pdf(paste0(dirnames$output,"/MixedEffectsTopRankedGenotype.pdf"), width = 14, height = 8)
plot_grid(plotlist = topRankedPlots_Genotype$rTg4510, nrow=2)
plot_grid(plotlist = topRankedPlots_Genotype$J20, nrow=2)
dev.off()

pdf(paste0(dirnames$output,"/MixedEffectsMainFigs.pdf"), width = 8, height = 6)
topRankedPlots_Tissue$rTg4510Geno[[1]]
topRankedPlots_Tissue$rTg4510Geno[[2]]
topRankedPlots_Tissue$rTg4510TissueGeno[[1]]
topRankedPlots_Tissue$rTg4510TissueGeno[[2]]
topRankedPlots_Tissue$J20TiOnly[[1]]
topRankedPlots_Tissue$J20Geno[[1]]
dev.off()


# tables
for(i in 1:4){
  write.table(stats[[i]], paste0(dirnames$output,names(stats)[[i]],"_array_mixedeffects.csv"),row.names = F,sep=",")
}


write.table(meta$rTg4510_ECX_Geno$betas, paste0(dirnames$output,"rTg4510_ECX_array_beta.csv"),row.names = T,sep=",")

