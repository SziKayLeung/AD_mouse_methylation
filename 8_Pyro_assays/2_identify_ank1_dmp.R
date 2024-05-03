#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Select differentially methylated positions associated with Ank1 for downstream validation with pyrosequencing
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Notes -----------------
## 
## Filtering Ank1-associated DMPs after reannotation with ChipSeeker and merging with gtf (1_reannotate_peaks.R)
## Filtering DMPs identified with genotype, interaction and pathology in rTg4510 model (RRBS)


##--------------- packages 

suppressMessages(library(dplyr))
suppressMessages(library(stringr))


##-------------- input -------

dirnames <- list(output = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/Figures/")

# read in files generated from 1_reannotate_peaks.R 
# genotype, interaction, pathology.csv
rrbs_annotated <- list.files(path = dirnames$output, pattern = "rrbs", full.names = T)
rrbs_annotated <- lapply(rrbs_annotated, function(x) read.csv(x))
names(rrbs_annotated) <- list.files(path = dirnames$output, pattern = "rrbs")


##-------------- selection of Ank1 probes -------

# filter for Ank1 in ChIPseeker_GeneSymbol across all .csv
# merge to one dataset
ank1_dmp <- lapply(rrbs_annotated, function(x) x %>% filter(ChIPseeker_GeneSymbol == "Ank1" | RefGtf_GeneSymbol == "Ank1") %>% 
                     select(Position, contains("Meth"), contains("FDR"), ChIPseeker_Annotation, ChIPseeker_GeneSymbol, RefGtf_GeneSymbol))


all_ank1_dmp <- dplyr::bind_rows(ank1_dmp, .id = "file")

# dissect postiion of probes to ease input to genome browser
all_ank1_dmp %>% mutate(chrom = word(Position,c(1),sep=fixed(":")),
                        probe = word(Position,c(2),sep=fixed(":")),
                        probe2 = word(Position,c(2),sep=fixed(":"))) 

# write output
write.csv(all_ank1_dmp, paste0(dirnames$output,"Ank1_DMP_positions.csv"), quote = F)