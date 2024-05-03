# Emma Walker
# E.M.Walker@exeter.ac.uk
# 07/02/2022

# This script annotates the rTg4510 RRBS data using the "Function_annotateRRBS.r" function
# It includes the "InGeneBodyOr1500bpTSS" column to specify sites which are either genic or 1.5kbp from the TSS
# And also a "InGeneBodyOr1500bpTSS_SYMBOL" which contains the gene symbol if InGeneBodyOr1500bpTSS == T


# title: "RRBS - BiSeq rTg4510 - 'Full' Annotations"
# author: "Emma Walker - adapted from Isabel Castanho (I.S.Castanho@exeter.ac.uk)"

# load packages
suppressMessages(library(stringr))
suppressMessages(library(GenomicRanges)) # GenomicRanges_1.38.0
suppressMessages(library(dplyr))
suppressMessages(library(ChIPseeker))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
suppressMessages(library(dplyr))
suppressMessages(library(naniar))


setwd("/lustre/projects/Research_Project-191406/EmmaW")
source("Function_annotateRRBS.r")


########################## DMRs #################################

# read in file locations
temp <- list.files("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMRs/rTg4510", full.names = T, pattern = ".csv")
temp2 <- gsub('.{4}$', '', list.files("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMRs/rTg4510/", pattern = ".csv"))

#remove empty interaction file as it causes issues in the loop
temp <- temp[-2]
temp2 <- temp2[-2]


#loop over all files to annotate them using the function and write to csv
for(i in 1:length(temp)){
  
  # Get file
  File <- read.csv(temp[i])
  
  # add columns to csv files for chromosome, location, min location and max location
  File$seqnames <- gsub("M", "MT", File$seqnames)
  File$location <-paste0(File$seqnames, ":", File$start, "-", File$end)
  
  
  # annotate file
  out <- annotateRRBS(File)
  
  #join original file and chipseeker annotations together
  out$location <- paste0(out$seqnames, ":", out$start, "-", out$end)
  out <- left_join(File, out, by="location")
  
  
  #write csv file
  outpath <- paste0("/gpfs/ts0/projects/Research_Project-191406/EmmaW/RRBSAnnotatedResults/rTg4510/DMRs/", temp2[i], "_1500bptssAnno.csv")
  write.csv(out, file = outpath)
}




############################### DMPs ############################


## Get significant locations and annotate them
temp <- list.files("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMPs/rTg4510/", full.names = T)
temp2 <- gsub('.{4}$', '', list.files("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/DMPs/rTg4510/", pattern = ".csv"))


#loop over all files to annotate them using the function and write to csv
for(i in 1:length(temp)){
  
  # Get file
  File <- read.csv(temp[i])
  File$location <- File$X
  
  # annotate file
  out <- annotateRRBS(File)
  
  #join original file and chipseeker annotations together
  out$location <- paste0(out$seqnames, ":", out$start)
  out <- left_join(File, out, by="location")
  

  outpath <- paste0("/gpfs/ts0/projects/Research_Project-191406/EmmaW/RRBSAnnotatedResults/rTg4510/DMPs/", temp2[i], "_1500bptssAnno.csv")
  write.csv(out, file = outpath)
  
}


