#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: identify the dataset (Iso-Seq or ONT) from isoform ID
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------
##



## -------- Load packages & source functions ---------

suppressMessages(
  {library("dplyr"); 
  library("cowplot")}
)

# directory paths 
dirnames <- list(
  source = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/9_methylation_expression_integration",
  rnaseq = "/gpfs/ts0/projects/Research_Project-191406/EmmaW/RNA_expr_comparisons/",
  annotated = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/0_ZenOutput/2_annotated/",
  rrbs = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/",
  figures = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/Figures/"
)


source(paste0(dirnames$source,"/0_simpleCorrelation_functions.R"))

# set colours for plots
color_rTg4510_TG <- "#00AEC9"
color_J20_TG <- "#FF5A62"
WT <- "#000000"

## -------- read in RRBS methylation --------- 

# rTg4510 & J20
load(paste0(dirnames$rrbs, "rTg4510_RRBSbetasComplete.RData")); rTg4510_RRBS <- RRBS_completebetas; rm(RRBS_completebetas)
load(paste0(dirnames$rrbs, "J20_RRBSbetasComplete.RData")); J20_RRBS <- RRBS_completebetas; rm(RRBS_completebetas)


## -------- read in array methylation --------- 

load(paste0(dirnames$arrayRaw,"ArrayBetaRegResultsFull_rTg4510.RData"))


## -------- read in gene normalized counts --------- 

input_RNA <- list(
  rTg4510 = read.csv(paste0(dirnames$rrbs, "RNA-seq/rTg4510_normalized_counts.csv"), stringsAsFactors = F),
  J20 = read.csv(paste0(dirnames$rrbs, "RNA-seq/J20_normalized_counts.csv"), stringsAsFactors = F)
)


## -------- read in phenotype --------- 

input_colData <- list(
  rTg4510 = read.csv(paste0(dirnames$rrbs, "Tg4510_coldata_RRBS.csv"), stringsAsFactors = F),
  J20 = read.csv(paste0(dirnames$rrbs, "J20_coldata_RRBS.csv"), stringsAsFactors = F)
)


## -------- read in results file of signficant differential gene expression --------- 

input_DGE <- list(
  rTg4510 = list(
    genotype = paste0(dirnames$rnaseq, "TableS2_rTg4510_genotype.csv"),
    interaction = paste0(dirnames$rnaseq, "TableS4_rTg4510_interaction.csv"),
    pathology = paste0(dirnames$rnaseq, "TableS5_rTg4510_ECX_pathology.csv")
  ),
  J20 = list(
    genotype = paste0(dirnames$rnaseq, "TableS2_J20_genotype.csv"),
    interaction = paste0(dirnames$rnaseq, "TableS4_J20_interaction.csv"),
    pathology = paste0(dirnames$rnaseq, "TableS5_J20_ECX_pathology.csv")
  )
)

input_DGE$rTg4510 <- lapply(input_DGE$rTg4510, function(x) read.csv(x, stringsAsFactors = F))
input_DGE$rTg4510$pathology <- input_DGE$rTg4510$pathology %>% dplyr::rename("Gene" = "X")

input_DGE$J20 <- lapply(input_DGE$J20, function(x) read.csv(x, stringsAsFactors = F))
input_DGE$J20$pathology <- input_DGE$J20$pathology %>% dplyr::rename("Gene" = "X")


## -------- read in results file of signficant DMPs --------- 

## rrbs
input_DMP <- list(
  rTg4510 = list(
    genotype = paste0(dirnames$annotated, "rTg4510_rrbs_genotype.csv"),
    interaction = paste0(dirnames$annotated, "rTg4510_rrbs_interaction.csv"),
    pathology = paste0(dirnames$annotated, "rTg4510_rrbs_pathology.csv")
  ),
  J20 = list(
    genotype = paste0(dirnames$annotated, "J20_rrbs_genotype.csv"),
    interaction = paste0(dirnames$annotated, "J20_rrbs_interaction.csv"),
    pathology = paste0(dirnames$annotated, "J20_rrbs_pathology.csv")
  )
)
input_DMP$rTg4510 <- lapply(input_DMP$rTg4510, function(x) read.csv(x, stringsAsFactors = F))
input_DMP$J20 <- lapply(input_DMP$J20, function(x) read.csv(x, stringsAsFactors = F))

## array 
input_DMP_array <- list(
  rTg4510 = list(
    genotype = paste0(dirnames$output, "rTg4510_ECX_Geno_array_mixedeffects.csv")
  ),
  J20 = list(
    genotype = paste0(dirnames$output, "J20_ECX_Geno_array_mixedeffects.csv")
  )
)
input_DMP_array$rTg4510 <- lapply(input_DMP_array$rTg4510, function(x) read.csv(x, stringsAsFactors = F) %>% 
                                    dplyr::rename("ChIPseeker_GeneSymbol"="SYMBOL","FDR_adj_genotype"="FDR_adj_GenotypeTG")) 
input_DMP_array$J20 <- lapply(input_DMP_array$J20, function(x) read.csv(x, stringsAsFactors = F) %>% 
                                dplyr::rename("ChIPseeker_GeneSymbol"="SYMBOL","FDR_adj_genotype"="FDR_adj_GenotypeTG"))


## -------- find common genes and plot --------- 

listModelTerm <- c("genotype","interaction")

# rTg4510
rTg4510 <- list()
for(m in listModelTerm){
  rTg4510[[m]] <- findCommonGene(res_DMP = input_DMP[["rTg4510"]][[m]], res_DGE = input_DGE[["rTg4510"]][[m]],
                                 ModelTerm=m, MouseModel = "rTg4510", 
                                 RNA = input_RNA$rTg4510, RRBS = rTg4510_RRBS, colData = input_colData$rTg4510)
}
names(rTg4510) = listModelTerm

# J20
listModelTerm <- c("genotype","interaction","pathology")
J20 <- list()
for(m in listModelTerm){
  J20[[m]] <- findCommonGene(res_DMP = input_DMP[["J20"]][[m]], res_DGE = input_DGE[["J20"]][[m]],
                                 ModelTerm=m, MouseModel = "J20", 
                                 RNA = RNA$J20, RRBS = input_array$rTg4510, colData = colData$J20)
}
names(J20) = listModelTerm

# pathology
length(setdiff(input_DMP$rTg4510$pathology$ChIPseeker_GeneSymbol,input_DGE$rTg4510$pathology$Gene))
topPathology <- input_DMP$rTg4510$pathology %>% arrange(FDR_adj_pathology) %>% .[,c("ChIPseeker_GeneSymbol")]
rTg4510$pathology <- list()
for(gene in c("Cacna1i","Prkar1a","Tgfb1","Ttll11")){
  rTg4510$pathology[[gene]] <- RNAcorMethPlotGenePathology(res = input_DMP$rTg4510$pathology, ModelTerm="pathology",  
                                RNA = input_RNA$rTg4510, RRBS = rTg4510_RRBS, colData = input_colData$rTg4510, gene, "rTg4510")
}

# array

findCommonGene(res_DMP = input_DMP_array[["rTg4510"]][[m]], res_DGE = input_DGE[["rTg4510"]][[m]],
               ModelTerm=m, MouseModel = "rTg4510", 
               RNA = input_RNA$rTg4510, RRBS = rTg4510_RRBS, colData = input_colData$rTg4510)

## -------- plot and save --------- 

# genotype
pdf(paste0(dirnames$figures, "corr_expression_methylation_byprobeGenotype.pdf"),width=20,height=10)
for(i in 1:length(rTg4510$genotype)){print(plot_grid(plotlist = rTg4510$genotype[[i]]))}
dev.off()

# interaction
pdf(paste0(dirnames$annotated, "corr_expression_methylation_byprobeInteraction.pdf"),width=20,height=10)
for(i in 1:length(rTg4510$interaction)){print(plot_grid(plotlist = rTg4510$interaction[[i]]))}
dev.off()

# pathology (subset)
pdf(paste0(dirnames$annotated, "corr_expression_methylation_byprobePathology.pdf"),width=20,height=10)
for(i in 1:length(rTg4510$pathology)){print(plot_grid(plotlist = rTg4510$pathology[[i]]))}
dev.off()



res <- res %>% mutate(ucsc = paste0(word(Position,c(1),sep=fixed(":"))," ",word(Position,c(2),sep=fixed(":")), " ", word(Position,c(2),sep=fixed(":"))))
res[res$ChIPseeker_GeneSymbol %in% inBoth,"ucsc"]
write.csv(res[res$ChIPseeker_GeneSymbol %in% inBoth,"ucsc"],
          "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/Figures/rTg4510_integration_probes.csv",
          row.names = F, quote = F, col.names = F)