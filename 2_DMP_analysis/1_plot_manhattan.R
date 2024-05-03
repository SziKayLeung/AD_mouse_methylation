#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Re-plot manhattan plots
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
## 
## Function
## ChipAnnotatePeaks - original function called for chipseeker annotation used by E.Walker
## 


#-------------- packages -------------
library("dplyr")
library("CMplot")
library("stringr")
library("wesanderson")
library("cowplot")

source("/lustre/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/2_DMP_analysis/0_source_functions.R")
source("/lustre/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/methylation.config.R")

##-------------- input -------

### Results after applying BiSeq and smoothing
# load from two mouse models and save into list
RRBS_complete <- list()
load(file = paste0(dirnames$biseq, "/rTg4510_RRBSbetasComplete.RData"))
RRBS_complete$rTg4510 <- RRBS_completebetas
load(file = paste0(dirnames$biseq, "/J20_RRBSbetasComplete.RData"))
RRBS_complete$J20 <- RRBS_completebetas

# number of sites before FDR correction
cat("Number of sites in rTg4510:", nrow(RRBS_complete$rTg4510),"\n")
cat("Number of sites in J20:", nrow(RRBS_complete$J20),"\n")


#-------------- processed results ----------

# manhattan plots 
# remmove chromosome X, Y and M for plotting
prep_mhat <- list(
  array = lapply(anno$array, function(x) prepare_manhattan(x, "PrZ.GenotypeTG","FDR_correct")),
  rrbs = lapply(anno$rrbs, function(x) prepare_manhattan(x, "FDR_adj_genotype","FDR_present"))
)

# remove top ranked rTg4510 ECX point to avoid skewing graph
prep_mhat$rrbs$rTg4510_ECX_geno <- prep_mhat$rrbs$rTg4510_ECX_geno %>% filter(Location != "chr12:115283603")
prep_mhat$array$rTg4510_ECX_geno <- prep_mhat$array$rTg4510_ECX_geno %>% filter(!Location %in% c("chr11:104318231","chr12:80436248","chr14:124676565","chr12:116425797"))

prepped_mhat <- list(
  rTg4510_ECX_geno = merge_manhattan(prep_mhat$rrbs$rTg4510_ECX_geno,prep_mhat$array$rTg4510_ECX_geno,"mFDR_RRBS","mFDR_Array"),
  J20_ECX_geno = merge_manhattan(prep_mhat$rrbs$J20_ECX_geno,prep_mhat$array$J20_ECX_geno,"mFDR_RRBS","mFDR_Array"),
  rTg4510_HIP_geno = prep_mhat$array$rTg4510_HIP_geno %>% mutate(gene = SYMBOL),
  J20_HIP_geno = prep_mhat$array$J20_HIP_geno %>% mutate(gene = SYMBOL)
)


cumu_ls <- lapply(prepped_mhat, function(x) cumulative(x))

final_mhat_p <- list(
  rTg4510_ECX_full = plot_manhattan_3(prepped_mhat$rTg4510_ECX_geno,"ECX","Tg4510","no"),
  rTg4510_HIP_full = plot_manhattan_3(cumu_ls$rTg4510_HIP_geno,"HIP","Tg4510","no"),
  rTg4510_ECX = plot_manhattan_3(prepped_mhat$rTg4510_ECX_geno,"ECX","Tg4510","yes"),
  rTg4510_HIP = plot_manhattan_3(cumu_ls$rTg4510_HIP_geno,"HIP","Tg4510","yes"),
  J20_ECX_full = plot_manhattan_3(prepped_mhat$J20_ECX_geno,"ECX","J20","no"),
  J20_HIP_full = plot_manhattan_3(cumu_ls$J20_HIP_geno,"HIP","J20","no"),
  J20_ECX = plot_manhattan_3(prepped_mhat$J20_ECX_geno,"ECX","J20","yes"),
  J20_HIP = plot_manhattan_3(cumu_ls$J20_HIP_geno,"HIP","J20","yes") 
)



pdf(paste0(output_dir,"/manhattan.pdf"), width = 10, height = 8)
# rTg4510
plot_manhattan_2(prep_mhat$rrbs$rTg4510_ECX_geno,"Tg4510","single","rTg4510 ECX RRBS - Genotype")
plot_manhattan_2(prep_mhat$array$rTg4510_ECX_geno,"Tg4510","single","rTg4510 ECX Array - Genotype")
plot_manhattan_2(prepped_mhat$rTg4510_HIP_geno,"Tg4510","single","rTg4510 HIP Array - Genotype")
# J20
plot_manhattan_2(prep_mhat$rrbs$J20_ECX_geno,"J20","single","J20 ECX RRBS - Genotype")
plot_manhattan_2(prep_mhat$array$J20_ECX_geno,"J20","single","J20 ECX Array - Genotype")
plot_manhattan_2(prepped_mhat$J20_HIP_geno,"J20","single","J20 HIP Array - Genotype")
##
# p7: rTg4510 ECX not filtered by pvalue (p<e-50), rTg4510 HIP not filtered by p-value (p<e-50)
# p8: rTg4510 ECX filtered by pvalue, rTg4510 HIP filtered by pvalue
# p9: J20 ECX not filtered by pvalue, J20 HIP not filtered by p-value
# p10: J20 ECX filtered by pvalue, J20 HIP filtered by p-value
plot_grid(final_mhat_p$rTg4510_ECX_full,final_mhat_p$rTg4510_HIP_full,ncol=1,nrow=2)
plot_grid(final_mhat_p$rTg4510_ECX,final_mhat_p$rTg4510_HIP,ncol=1,nrow=2)
plot_grid(final_mhat_p$J20_ECX_full,final_mhat_p$J20_HIP_full,ncol=1,nrow=2)
plot_grid(final_mhat_p$J20_ECX,final_mhat_p$J20_HIP,ncol=1,nrow=2)
dev.off()

# output list
final_plots <- list(
  # rTg4510: top = ECX, bottom = HIP
  rTg4510 = plot_grid(plot_manhattan_3(prepped_mhat$rTg4510_ECX_geno,"ECX_upwards","Tg4510","yes") +
              scale_y_continuous(lim=c(0,8), breaks = seq(0, 8, 2)),
            plot_manhattan_3(cumu_ls$rTg4510_HIP_geno %>% mutate(mFDR),"HIP","Tg4510","yes") + 
              scale_y_reverse(lim=c(50,0), breaks = seq(50, 0, by = -10)), ncol = 1),
  
  # J20: top = ECX, bottom = HIP
  J20 = plot_grid(plot_manhattan_3(prepped_mhat$J20_ECX_geno,"ECX_upwards","J20","yes") +
              scale_y_continuous(lim=c(0,8), breaks = seq(0, 8, 2)),
            plot_manhattan_3(cumu_ls$J20_HIP_geno%>% mutate(mFDR),"HIP","J20","yes") + 
              scale_y_reverse(lim=c(8,0), breaks = seq(8, 0, by = -2)), ncol = 1)
  
)

pdf(paste0(output_dir,"/final_manhattan.pdf"), width = 10, height = 8)
final_plots$rTg4510
final_plots$J20
# ECX genotype
plot_grid(
  plot_manhattan_3(cumu_ls$rTg4510_ECX_geno %>% mutate(mFDR = mFDR_RRBS) %>% filter(Location != "chr12:115283603"), "upwards","Tg4510","no") + 
    scale_y_continuous(lim=c(0,8), breaks = seq(0, 16, 2)),
  plot_manhattan_3(cumu_ls$rTg4510_ECX_geno %>% mutate(mFDR = mFDR_RRBS) %>% filter(Location != "chr12:115283603"), "upwards","Tg4510","no") + scale_y_continuous(lim=c(0,8), breaks = seq(0, 16, 2)),
  plot_manhattan_3(cumu_ls$J20_ECX_geno %>% mutate(mFDR = mFDR_RRBS),"HIP","J20","no") + scale_y_reverse(lim=c(8,0), breaks = seq(8, 0, by = -2)),
  ncol = 1)

dev.off()

names(final_plots) <- c("rTg4510manhattan","J20manhattan")
saveRDS(final_plots, paste0(dirnames$output,"PlotManhattan.rds"))
