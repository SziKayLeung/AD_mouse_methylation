#!/usr/bin/env Rscript
## ----------Script-----------------  
##
## Purpose: Re-annotate differentially-methylated positions using chipseeker and alignment of gtf
##          Barplot of chipseeker annotations in rTg4510 and J20 mouse model (location)
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------
## 
## Motivation
## 1. E.Walker used chipseeker for annotation of DMP
## 2. I.Castanho noticed misannotation of top DMP (in exon) to Mir6388 rather than Gm30948.
## 3. S.Leung to validate and correct annotations (see ChipSeeker_wrongannotations.html)
## 4. S.Leung noticed another separate method for annotation, by annotating the EnsemblID in the the “annotation” column from ChipSeeker output with mouse reference gtf
##    Such approach, while limited to DMPs located in exon and intron, further annotated DMP to Gm30948     
## 5. After discussion with J.Mill, E.Hannon..., to include column for Chipseeker annotations and annotation from merging with reference gtf
## 6. Use Chipseeker column to manually output annotation locations in order to include results across all and significant sites, and in both mouse models
## 
## Function
## ChipAnnotatePeaks - original function called for chipseeker annotation used by E.Walker
## 

##--------------- packages 

suppressMessages(library(stringr))
suppressMessages(library(GenomicRanges)) # GenomicRanges_1.38.0
suppressMessages(library(dplyr))
suppressMessages(library(ChIPseeker))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
suppressMessages(library(dplyr))
suppressMessages(library(naniar))
suppressMessages(library(kableExtra))
suppressMessages(library(ggplot2))
suppressMessages(library(wesanderson))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# directory input
rootDir = "/lustre/projects/Research_Project-191406"
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/"
  
# source functions to run script
source(paste0(scriptDir, "4_DMP_annotation/0_source_functions.R"))

##-------------- input -------

# directory for rTg4510 and J20 RRBS and Genotype 
dirnames <- list(
  rTg4510_rrbs = paste0(rootDir, "/isabel/RRBS_new/DMPs/rTg4510"),
  J20_rrbs = paste0(rootDir, "/isabel/RRBS_new/DMPs/J20"),
  mouse_array = paste0(rootDir, "/Aisha/data/Array"),
  EW_results = paste0(rootDir, "/EmmaW/RRBSAnnotatedResults/"),
  ref_dir = "/lustre/projects/Research_Project-MRC148213/sl693/reference/annotation",
  output = paste0(scriptDir, "0_ZenOutput/2_annotated/")
)


# results after applying BiSeq and smoothing
# load from two mouse models and save into list
RRBS_complete <- list()
load(file = paste0(rootDir, "/isabel/RRBS_new/rTg4510_RRBSbetasComplete.RData"))
RRBS_complete$rTg4510 <- RRBS_completebetas
load(file = paste0(rootDir, "/isabel/RRBS_new/J20_RRBSbetasComplete.RData"))
RRBS_complete$J20 <- RRBS_completebetas
RRBS_complete <- lapply(RRBS_complete, function(x) x %>% tibble::rownames_to_column(., var = "position"))

# read input files
input_files <- list(
  # genotype 
  rTg4510_rrbs_genotype = paste0(dirnames$rTg4510_rrbs, "/DMPsInteractionModel_rTg4510_sig_genotype.csv"),
  J20_rrbs_genotype = paste0(dirnames$J20_rrbs, "/DMPsInteractionModel_J20_sig_genotype.csv"),
  # interaction
  rTg4510_rrbs_interaction = paste0(dirnames$rTg4510_rrbs, "/DMPsInteractionModel_rTg4510_sig_interaction.csv"),
  J20_rrbs_interaction = paste0(dirnames$J20_rrbs, "/DMPsInteractionModel_J20_sig_interaction.csv"),
  # pathology 
  rTg4510_rrbs_pathology = paste0(dirnames$rTg4510_rrbs, "/DMPsPathology_rTg4510_sig_pathology.csv"),
  J20_rrbs_pathology = paste0(dirnames$J20_rrbs, "/DMPsPathology_J20_sig_pathology.csv"),
  # reference
  mouse_gtf = paste0(dirnames$ref_dir, "/gencode.vM22.annotation.gtf")
)
input_reg <- lapply(input_files[1:6], read.csv)
input_reg <- lapply(input_reg, function(x) x %>% mutate(position = X) %>% dplyr::select(-X, -betaResultsClean....Position..))
input_gtf <- rtracklayer::import(input_files$mouse_gtf) # creates a GRanges object
input_gtf <- as.data.frame(input_gtf)

input_reg$rTg4510_all <- RRBS_complete$rTg4510
input_reg$J20_all <- RRBS_complete$J20


##-------------- annotate peaks using ChipSeeker

anno_chipseeker <- lapply(input_reg, function(x) ChipAnnotatePeaks(x))


##-------------- classify and plot annotations

# classify peaks into simplified categories 
class_anno_chipseeker <- lapply(anno_chipseeker, function(x) ChipClassifyAnno_df(x))

# plot annotations of all and significant sites using chipseeker annotation column
pAnno <- ChipPlotAnno(class_anno_chipseeker)


##-------------- further annotate peaks in intron and exon (defined from ChipSeeker) using reference gtf

anno_chipseeker_mergedref <- lapply(anno_chipseeker, function(x) ChipMergeRefGtf(x, input_gtf))


##-------------- finalise annotations from ChipSeeker and merging with reference gtf

anno_final <- list()
for(i in 1:6){anno_final[[i]] <- ChipFinaliseAnno(anno_chipseeker[[i]], anno_chipseeker_mergedref[[i]], input_reg[[i]])}
names(anno_final) <- names(anno_chipseeker)[1:6]


##-------------- QC: crossreference ChipSeeker results with Emma's latest results (April 2022)

# read in Emma Walker's rTg4510 genotype results 
EW_Chipseeker <- list()
EW_Chipseeker$rTg4510_genotype = paste0(dirnames$EW_results, "rTg4510/DMPs/DMPsInteractionModel_rTg4510_sig_genotype_1500bptssAnno.csv")

# create a column for comparison: location_gene = <position>_<chipseeker_output_gene>
EW_Chipseeker <- lapply(EW_Chipseeker, function(x) read.csv(x) %>% mutate(location_gene = paste0(X,"_",ENSEMBL)))
anno_final$rTg4510_rrbs_genotype <- anno_final$rTg4510_rrbs_genotype %>% mutate(location_gene = paste0(Position,"_",ChIPseeker_GeneEnsembl))

# check no annotations that are misaligned between latest results (Jan 2023) and Emma's results (April 2022)
if(length(setdiff(EW_Chipseeker$rTg4510_genotype$location_gene, anno_final$rTg4510_rrbs_genotype$location_gene)) > 0){print("Missing annotations in current vs 2022 results")}
if(length(setdiff(anno_final$rTg4510_rrbs_genotype$location_gene, EW_Chipseeker$rTg4510_genotype$location_gene)) > 0){print("Additional annotations in current vs 2022 results")}                             


##-------------- Write output to xlsx file 

for(i in 1:6){
  write.csv(anno_final[[i]], paste0(dirnames$output, names(anno_final)[[i]],".csv"))
}

pdf(paste0(dirnames$output,"/ChipseekerAnnotation.pdf"), width = 10, height = 8)
pAnno
dev.off()