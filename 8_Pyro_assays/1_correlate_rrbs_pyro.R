#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: correlate pyrosequencing data and RRBS data (complete and smoothed)
## 1 Bin1 DMP 
## 4 Prnp DMPs (in DMR)
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------------------------


# source input variables
library(ggrepel)
library(stringr)
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/8_Pyro_assays/pyro.config.R")
color_Tg4510_TG <- "#00AEC9"

# subset rrbs data by position (referring to validated DMP using pyrosequencing)
valRRBS <- list(
  bin1 = rrbs %>% filter(position == "chr18:32431597"),
  prnp1 = rrbs %>% filter(position == "chr2:131910163"),
  prnp2 = rrbs %>% filter(position == "chr2:131910165"),
  prnp3 = rrbs %>% filter(position == "chr2:131910181"),
  prnp4 = rrbs %>% filter(position == "chr2:131910202"),
  ank1 = rrbs %>% filter(position == "chr8:23023193"),
  ank2 = rrbs %>% filter(position == "chr8:23023211"),
  ank3 = rrbs %>% filter(position == "chr8:23023241")
)

valRawRRBS <- list(
  ank1 = RRBS_rawbetas %>% filter(position == "chr8:23023193"),
  ank2 = RRBS_rawbetas %>% filter(position == "chr8:23023211"),
  ank3 = RRBS_rawbetas %>% filter(position == "chr8:23023241")
)
valRawRRBS <- lapply(valRawRRBS, function(x) x %>% reshape2::melt(variable.name = "sample", value.name = "rrbs_perc", id = c("position")))


plot_rrbs_across_age <- function(input.position){
  
  subset.rrbs <- rrbs %>% filter(position == input.position) %>% filter(!sample %in% c("K18","L24","M22","N20"))
  merged.dat <- merge(input_files$RRBS_phenotype[,c("Sample_ID","Genotype","Age_months")], subset.rrbs, by.x="Sample_ID",by.y="sample",all.y=T) 
  
  p <- merged.dat %>%
    mutate(run = ifelse(Sample_ID %in% samples2run,TRUE,FALSE)) %>%
    mutate(Genotype = factor(Genotype,levels=c("WT","TG"))) %>%
    ggplot(., aes(x = Genotype,y=rrbs_perc)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(colour=run)) +
    labs(x = "Genotype", y = "Methylation (%)", subtitle = input.position) + theme_classic() + 
    facet_grid(~Age_months) +
    geom_label_repel(aes(label = Sample_ID, fill = run),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50')
  
  return(p)
  
}


## ------------------- corr_rrbs_pyro()

# Aim: correlate the pyrosequencing and rrbs values 
# Input:
  # rrbs = df of subsetted position of DMP; cols: <position> <chr> <coordinate> <sample> <rrbs_perc>
  # pyroInput = df of pyrosequencing results; cols : <SAMPLE> <Pos1Meth> <Pos1Quality> <Pos2Meth> 
  # phenotype = df of phenotype data; cols <Sample_ID> <Age_months> <Genotype>
  # positionname = string to extract columns in pyroInput (i.e Pos1, Pos2)
  # noquality = TRUE/FALSE; whether to do filtering on the quality of pyrosequencing data
# Output:
  # scatter plot of the methylation percentage using RRBS (y-axis) and pyrosequencing (x-axis), and correlation

corr_rrbs_pyro <- function(rrbs, pyroInput, phenotype, positioname, noquality){
  # select position name of interest
  pyro_pos <- pyroInput %>% select(SAMPLE, contains({{ positioname }}))
  print(head(pyro_pos))
  
  # remove samples that failed
  if(noquality == "FALSE"){
    pyro_pos <- pyro_pos %>% 
      mutate(quality = ifelse(.[[3]] == "Passed","Passed","Failed")) %>% 
      filter(quality == "Passed")
  }
  
  # merge with rrbs values
  merged_pyro_pos <- merge(rrbs, pyro_pos, by.x = "sample", by.y = "SAMPLE") %>% mutate(rrbs_perc = rrbs_perc * 100) 
  
  merged_pyro_pos <- merge(merged_pyro_pos, phenotype[,c("Sample_ID","Age_months","Genotype")], 
                            by.x = "sample", by.y = "Sample_ID") 
  
  # correlation
  x.var <- rlang::sym(paste0(quo_name(enquo(positioname)),"Meth"))
  merged_pyro_pos[[x.var]] <- as.numeric(as.character(merged_pyro_pos[[x.var]]))

  print(cor.test(merged_pyro_pos[[x.var]],merged_pyro_pos[["rrbs_perc"]], use = "pairwise.complete.obs"))
  corr.value <- cor(merged_pyro_pos[[x.var]],merged_pyro_pos[["rrbs_perc"]], use = "pairwise.complete.obs")
  corr <- grobTree(textGrob(paste("r = ", round(corr.value, 2)), 
                            x = 0.05, y = 0.90, hjust = 0, 
                            gp = gpar(col = "black", fontsize = 14, fontface = "italic",family="CM Roman")))
  
  
  # convert age to factor for shape plotting
  merged_pyro_pos <- merged_pyro_pos %>% mutate(Age_months = as.factor(Age_months))
  
  # plot 
  p <- merged_pyro_pos %>% 
    ggplot(., aes(x = !! x.var, y = rrbs_perc, colour = Genotype, shape = Age_months)) + geom_point(size = 3) +
    geom_label_repel(aes(label = sample), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') + 
    labs(x = "Pyrosequencing - methylation (%)", y = "RRBS - methylation (%)", title = unique(merged_pyro_pos$position)) +
    annotation_custom(corr) + theme_bw()
  
  return(p)
  
}


## ------------------- output

# Bin1 
corr_rrbs_pyro(valRRBS$bin1, input_pyro$bin1, input_files$RRBS_phenotype, "Pos1", noquality="FALSE")

# Ank1
Ank1p1 <- corr_rrbs_pyro(valRRBS$ank1, input_pyro$ank1, input_files$RRBS_phenotype, "Pos1", noquality="TRUE")
Ank1p2 <- corr_rrbs_pyro(valRRBS$ank2, input_pyro$ank1, input_files$RRBS_phenotype, "Pos2", noquality="TRUE")
Ank1p3 <- corr_rrbs_pyro(valRRBS$ank3, input_pyro$ank1, input_files$RRBS_phenotype, "Pos3", noquality="TRUE")

# rawBeta - Ank1
Ank1p4 <- corr_rrbs_pyro(valRawRRBS$ank1, input_pyro$ank1, input_files$RRBS_phenotype, "Pos1", noquality="TRUE")
Ank1p5 <- corr_rrbs_pyro(valRawRRBS$ank2, input_pyro$ank1, input_files$RRBS_phenotype, "Pos2", noquality="TRUE")
Ank1p6 <- corr_rrbs_pyro(valRawRRBS$ank3, input_pyro$ank1, input_files$RRBS_phenotype, "Pos3", noquality="TRUE")


# Pnrp
p_corrPrnp <- list(
  corr_rrbs_pyro(valRRBS$prnp1, input_pyro$prnp, input_files$RRBS_phenotype, "Pos1", noquality="TRUE"),
  corr_rrbs_pyro(valRRBS$prnp2, input_pyro$prnp, input_files$RRBS_phenotype, "Pos2", noquality="TRUE"),
  corr_rrbs_pyro(valRRBS$prnp3, input_pyro$prnp, input_files$RRBS_phenotype, "Pos3", noquality="TRUE"),
  corr_rrbs_pyro(valRRBS$prnp4, input_pyro$prnp, input_files$RRBS_phenotype, "Pos4", noquality="TRUE")
)

plot_grid(plotlist = p_corrPrnp)
pdf(paste0(dirnames$output,"/SuppFigures_Ank1PyroCorr.pdf"), width = 6, height = 6)
Ank1p1
Ank1p2
Ank1p3
Ank1p4
Ank1p5
Ank1p6
dev.off()