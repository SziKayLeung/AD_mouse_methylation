oldRRBS <- read.csv(paste0(dirnames$annotated, "/rTg4510_rrbs_genotype.csv"))
oldRRBS[oldRRBS$Position == "chr19:46730273",]
oldRRBS[oldRRBS$Position == "chr8:23023241",c("Position","FDR_adj_genotype","ChIPseeker_GeneSymbol","meth.group1.WT","meth.group2.TG")]
oldRRBS[oldRRBS$Position == "chr8:23023240",c("Position","FDR_adj_genotype","ChIPseeker_GeneSymbol","meth.group1.WT","meth.group2.TG")]


rTg4510_DMP$Genotype[rTg4510_DMP$Genotype$Position == "chr19:46730273",]

load(paste0(dirnames$processed, "/rTg4510_RRBS_SmoothBetas.RData"))
RRBS_smoothbetas["chr19:46730273",]
RRBS_smoothbetas["chr8:23023240",]
RRBS_smoothbetas["chr8:23023241",]

load(file = "/lustre/projects/Research_Project-191406/isabel/RRBS_new/rTg4510_RRBSbetasComplete.RData")
merge(RRBS_smoothbetas["chr8:23023241",] %>% reshape2::melt(variable.name = "sample", value.name = "methylation"), 
      phenotype$rTg4510, by.x = "sample", by.y = 0) %>% 
  ggplot(., aes(x = Genotype, y = methylation)) + geom_boxplot()


merge(RRBS_smoothbetas_SK["chr8:23023240",] %>% reshape2::melt(variable.name = "sample", value.name = "methylation"), 
      phenotype$rTg4510, by.x = "sample", by.y = 0) %>% 
  ggplot(., aes(x = Genotype, y = methylation)) + geom_boxplot()


as.data.frame(rowRanges(rrbs.clust.unlim)) %>% filter(start == 46730273)
as.data.frame(rowRanges(rrbs)) %>% filter(start == 46730273)

# raw coverage data 
load(paste0(dirnames$processed,"/biseq_raw_rTg4510.RData"))
RRBS_rawbetas["chr8:23023240",]


rrbs.reduced <- filterBySharedRegions(object=rrbs, groups=colData(rrbs)$Genotype, perc.samples=0.1, minCov=1)
rrbs.rel <- rawToRel(rrbs.reduced) # convert to methylation values
RRBS_rawbetas <- as.data.frame(methLevel(rrbs.rel))
coordinates <- as.data.frame(rrbs.reduced@rowRanges)
rownames(RRBS_rawbetas) <- paste0(coordinates$seqnames, ":", coordinates$start)
RRBS_rawbetas_SKL <- RRBS_rawbetas

load(file = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/biseq_rTg4510.RData")
rrbs_IC <- rrbs
dim(rrbs) # 4,793,943 sites
rrbs.reduced <- filterBySharedRegions(object=rrbs_IC, groups=colData(rrbs_IC)$Genotype, perc.samples=0.1, minCov=1)
rrbs.rel <- rawToRel(rrbs.reduced) # convert to methylation values
RRBS_rawbetas <- as.data.frame(methLevel(rrbs.rel))
nrow(RRBS_rawbetas) # 1,559,379 sites
coordinates <- as.data.frame(rrbs.reduced@rowRanges)
rownames(RRBS_rawbetas) <- paste0(coordinates$seqnames, ":", coordinates$start)
RRBS_rawbetas_IC <- RRBS_rawbetas

p1 <- merge(RRBS_rawbetas_IC["chr8:23023240",] %>% reshape2::melt(variable.name = "sample", value.name = "methylation"), 
      phenotype$rTg4510, by.x = "sample", by.y = 0) %>% 
  mutate(methylation = as.numeric(ifelse(methylation == "NaN", 0, methylation))) %>% 
  ggplot(., aes(x = Genotype, y = methylation)) + geom_boxplot() + labs(title = "chr8:23023240")

p2 <- merge(RRBS_rawbetas_IC["chr8:23023241",] %>% reshape2::melt(variable.name = "sample", value.name = "methylation"), 
      phenotype$rTg4510, by.x = "sample", by.y = 0) %>% 
  mutate(methylation = as.numeric(ifelse(methylation == "NaN", 0, methylation))) %>% 
  ggplot(., aes(x = Genotype, y = methylation)) + geom_boxplot() + labs(title = "chr8:23023241")


p3 <- merge(RRBS_rawbetas_SKL["chr8:23023240",] %>% reshape2::melt(variable.name = "sample", value.name = "methylation"), 
      phenotype$rTg4510, by.x = "sample", by.y = 0) %>% 
  mutate(methylation = as.numeric(ifelse(methylation == "NaN", 0, methylation))) %>% 
  ggplot(., aes(x = Genotype, y = methylation)) + geom_boxplot() + 
  labs(title = "chr8:23023240 -23023241 - Bismark mergeCpG")


p4 <- merge(RRBS_rawbetas_SKL[c("chr8:23023240","chr8:23023241"),] %>% reshape2::melt(variable.name = "sample", value.name = "methylation"), 
      phenotype$rTg4510, by.x = "sample", by.y = 0) %>% 
  mutate(methylation = as.numeric(ifelse(methylation == "NaN", 0, methylation))) %>% 
  ggplot(., aes(x = Genotype, y = methylation)) + 
  geom_boxplot() + 
  labs(title = "chr8:23023240,23023241 - manual merge")

plot_grid(p1,p2)
plot_grid(p3,p4)

RRBS_rawbetas_IC["chr8:23023241",]

RRBS_rawbetas_SKL <- RRBS_rawbetas
RRBS_rawbetas_SKL["chr8:23023240",]

plot_raw <- function(rawbetas, coordinate){
  
  p <- merge(rawbetas[coordinate,] %>% reshape2::melt(variable.name = "sample", value.name = "methylation"), 
        phenotype$rTg4510, by.x = "sample", by.y = 0) %>% 
    #mutate(methylation = as.numeric(ifelse(methylation == "NaN", 0, methylation))) %>% 
    ggplot(., aes(x = Genotype, y = methylation)) + geom_boxplot() + labs(title = coordinate)
  
  return(p)
}

p1 <- plot_raw(RRBS_rawbetas_IC, "chr19:46730291")
p2 <- plot_raw(RRBS_rawbetas_IC, "chr19:46730292")
p3 <- plot_raw(RRBS_rawbetas_SKL, "chr19:46730291")
plot_raw(RRBS_rawbetas_SKL, "chr8:23023240")

plot_grid(p1,p2)

compareRawBeta <- function(coordinate1, coordinate2){
  dat <- cbind(RRBS_rawbetas_IC[coordinate1,] %>% reshape2::melt(variable.name = "sample", value.name = paste0("methylation_",coordinate1)),
               RRBS_rawbetas_IC[coordinate2,] %>% reshape2::melt(variable.name = "sample2", value.name = paste0("methylation_",coordinate2)))
  final <- merge(dat, 
        RRBS_rawbetas_SKL[coordinate1,] %>% reshape2::melt(variable.name = "sample", value.name = paste0("methylationSL_", coordinate1)), by = "sample")
  
  return(final)
}

compareRawBeta("chr8:23023241","chr8:23023242")
compareRawBeta("chr19:46730291","chr19:46730292")

oldRRBS[oldRRBS$Position == "chr19:46730291",c("Position","FDR_adj_genotype","ChIPseeker_GeneSymbol","meth.group1.WT","meth.group2.TG")]
oldRRBS[oldRRBS$Position == "chr19:46730292",c("Position","FDR_adj_genotype","ChIPseeker_GeneSymbol","meth.group1.WT","meth.group2.TG")]


