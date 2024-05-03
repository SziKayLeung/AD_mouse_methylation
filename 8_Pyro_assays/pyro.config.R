#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: correlate rrbs and pyro effect size for prnp
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------


##--------------- packages 

suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(grid))


##-------------- input -------

dirnames <- list(
  biseq = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/",
  rTg4510_root =  "/lustre/projects/Research_Project-MRC148213/sl693/rTg4510/D_Methylation/",
  rTg4510_rrbs = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMPs/rTg4510/",
  metaroot = "/lustre/projects/Research_Project-MRC148213/sl693/rTg4510/0_metadata/D_methylation/",
  pyroRaw = "/lustre/projects/Research_Project-MRC148213/sl693/rTg4510/1_raw/G_pyro/",
  output = "/lustre/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/Figures/"
)
load(file = paste0(dirnames$biseq, "/rTg4510_RRBSbetasComplete.RData"))

# raw data before smoothing
load(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/biseq_rTg4510.RData")
# Filter raw methylation data by coverage
dim(rrbs) # 4,824,627 sites
rrbs.reduced <- filterBySharedRegions(object=rrbs, groups=colData(rrbs)$Genotype, perc.samples=0.5, minCov=10)
rrbs.rel <- rawToRel(rrbs.reduced) # convert to methylation values
RRBS_rawbetas <- as.data.frame(methLevel(rrbs.rel))
coordinates <- as.data.frame(rrbs.reduced@rowRanges)
RRBS_rawbetas$position <- paste(coordinates$seqnames, ":", coordinates$start, sep="")


# datawrangle smoothed data
rrbs <- RRBS_completebetas %>% tibble::rownames_to_column(., var = "position") 
rrbs$chrom <- word(rrbs$position,1,sep=fixed(":"))
rrbs$coordinate <-word(rrbs$position,2,sep=fixed(":"))
rrb <- rrbs %>% reshape2::melt() %>% select(position, chrom, coordinate, variable, value)
colnames(rrbs) <- c("position","chr","coordiate","sample","rrbs_perc")

# read input files
input_files <- list(
  
  # rT4510 full model: DMP data files
  rTg4510_rrbs_stats = paste0(dirnames$rTg4510_rrbs, "/DMPsInteractionModel_stats_table_rTg4510.csv"),
  
  # phenotype
  RRBS_phenotype = paste0(dirnames$metaroot, "Tg4510_phenotype_RRBS.csv")

)
input_files <- lapply(input_files, function(x) read.csv(x))

# input pyro raw data
input_pyro <- list(
  bin1 = read.csv(paste0(dirnames$pyroRaw, "Toni_bin1_failed.csv")),
  prnp = read.csv(paste0(dirnames$pyroRaw, "prnp_dmr.csv")),
  ank1 = read.csv(paste0(dirnames$pyroRaw, "ank1_dmp.csv"))
)
