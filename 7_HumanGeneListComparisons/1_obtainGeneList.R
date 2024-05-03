#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Generate gene lists associated to DMP from rTg4510 and J20 RRBS and genotype
##          Generate gene list from RSmith meta-analysis and GShireby FANS dataset
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------
## 
## 1. EWalker obtained meta-analysis results from RSmith and converted to human orthologue
## 2. EWalker ran array analysis for rTg4510 and J20 (note missing interaction for J20 array)
## 3. SLeung obtained gene list from rTg4510 and J20 RRBS analysis for genes with differentially methylated probes (using chipseeker)
 

#-------------- packages -------------

setwd("/gpfs/ts0/projects/Research_Project-191406/EmmaW/HumanGeneListComparisons/")
source("function_makeGeneList_RRBS.r")
source("function_cbind.fill.r")


#-------------- Bex's meta-analysis results -------------

# read in gene conversion file
mouse2human <- read.csv("mousehumangeneconversion.csv")

#read in results file from meta analysis
resFile <- read.csv("/gpfs/ts0/projects/Research_Project-191406/EmmaW/HumanGeneListComparisons/DMPs_meta2021_fixed_effect.csv", stringsAsFactors = F)

## reduce gene annotation to unique genes
resFile$UCSC_REFGENE_NAME<-unlist(lapply(resFile$UCSC_REFGENE_NAME, uniqueAnno))
## exclude sites with no gene anno
resFile<-resFile[-c(which(resFile$UCSC_REFGENE_NAME == "")),]
# create list of genes
human <-unique(unlist(strsplit(resFile$UCSC_REFGENE_NAME, "\\;")))

# convert to mouse
converted <- left_join(as.data.frame(human), mouse2human)
converted <- as.character(converted$mouse[!is.na(converted$mouse)])


#-------------- RRBS -------------

rrbsOutputDir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/Figures/"
rrbs_names <- list.files(rrbsOutputDir,pattern="rrbs",full.names = T)
rrbs <- lapply(rrbs_names, function(x) read.csv(x))
names(rrbs) <- word(list.files(rrbsOutputDir,pattern="rrbs"),c(1),sep=fixed("."))

rrbsGeneList <- list(
  rTg4510_rrbs_genotype = rrbs$rTg4510_rrbs_genotype %>% filter(FDR_adj_genotype < 0.05) %>% select(ChIPseeker_GeneSymbol),
  rTg4510_rrbs_interaction = rrbs$rTg4510_rrbs_interaction %>% filter(FDR_adj_interaction < 0.05) %>% select(ChIPseeker_GeneSymbol),
  rTg4510_rrbs_pathology = rrbs$rTg4510_rrbs_pathology %>% filter(FDR_adj_pathology < 0.05) %>% select(ChIPseeker_GeneSymbol),
  J20_rrbs_genotype = rrbs$J20_rrbs_genotype %>% filter(FDR_adj_genotype < 0.05) %>% select(ChIPseeker_GeneSymbol),
  J20_rrbs_interaction = rrbs$J20_rrbs_interaction %>% filter(FDR_adj_interaction < 0.05) %>% select(ChIPseeker_GeneSymbol),
  J20_rrbs_pathology = rrbs$J20_rrbs_pathology %>% filter(FDR_adj_pathology < 0.05) %>% select(ChIPseeker_GeneSymbol)
)


#-------------- Array -------------

arrayOutputDir = "/gpfs/ts0/projects/Research_Project-191406/EmmaW/Array/Results/rTg4510/SigResults/1500bpTSSAnnotations/"
array_names <- list.files(arrayOutputDir,pattern="ECX",full.names = T)[c(2,4,5)]
rTg4510_array <- lapply(array_names, function(x) read.csv(x))
names(rTg4510_array) <- c("rTg4510_array_interaction","rTg4510_array_genotype","rTg4510_array_pathology")

arrayOutputDir = "/gpfs/ts0/projects/Research_Project-191406/EmmaW/Array/Results/J20/SigResults/1500bpTSSAnnotations/"
array_names <- list.files(arrayOutputDir,pattern="ECX",full.names = T)[c(2,4)]
J20_array <- lapply(array_names, function(x) read.csv(x))
names(J20_array) <- c("J20_array_genotype","J20_array_pathology")

arrayGeneList <- lapply(rTg4510_array, function(x) x %>% filter(FDR.Adj.Pvalues < 0.05) %>% select(Gene_Symbol))
names(arrayGeneList) <- c("rTg4510_array_interaction","rTg4510_array_genotype","rTg4510_array_pathology")

arrayGeneList[4:5] <- lapply(J20_array, function(x) x %>% filter(FDR.Adj.Pvalues < 0.05) %>% select(Gene_Symbol))
names(arrayGeneList)[4:5] <-c("J20_array_genotype","J20_array_pathology")


#-------------- write out gene list -------------

for(i in 1:length(rrbsGeneList)){
  name = names(rrbsGeneList[i])
  write.table(rrbsGeneList[[i]], paste0(rrbsOutputDir,name,"_geneList.txt"),row.names=F,col.names=F,sep="\t",quote=F)
}

for(i in 1:length(arrayGeneList)){
  name = names(arrayGeneList[i])
  write.table(arrayGeneList[[i]], paste0(rrbsOutputDir,name,"_geneList.txt"),row.names=F,col.names=F,sep="\t",quote=F)
}
write.table(converted,paste0(rrbsOutputDir,"RSmith_human2mouse_geneList.txt"),row.names=F,col.names=F,sep="\t",quote=F)
