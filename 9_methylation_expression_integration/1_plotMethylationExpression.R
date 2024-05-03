library("dplyr")
library("ggplot2")
library("cowplot")

dirnames <- list(
  script = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/9_methylation_expression_integration/",
  biseq = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/",
  isabel = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/C_RNASeq/1_RNASeq_Isabel/",
  sigResults = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/Figures/",
  output = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/C_RNASeq/3_integration/"
)

source(paste0(dirnames$script,"0_MethReg_functions.R"))
source(paste0(dirnames$script,"0_MethReg_interactionModel.R"))
source(paste0(dirnames$script,"0_MethReg_plotInteractionModel.R"))
source(paste0(dirnames$script,"0_MethReg_stratifiedModel.R"))


# dna methylation
message("Read in DNA methylation")
load(file = paste0(dirnames$biseq, "rTg4510_RRBSbetasComplete.RData"))

# gene expression
message("Read in gene expression")
gene.exp <- read.csv(paste0(dirnames$isabel, "rTg4510_DESeq2normalizedcounts.csv"))

# phenotype 
coldata <- read.csv("/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/Tg4510_coldata_RRBS.csv", row.names=1, stringsAsFactors=FALSE)
coldata$Age_months <- as.factor(coldata$Age_months) # originally ran as.factor but we later changed to numeric
coldata$Genotype <- as.factor(coldata$Genotype)
coldata$Genotype <- relevel(coldata$Genotype, "WT")
levels(coldata$Genotype)

plot_integration_raw <- function(targetGene,tranFactor,probe){
  expG <- gene.exp[gene.exp$X == targetGene,] %>% reshape2::melt()
  expTF <- gene.exp[gene.exp$X == tranFactor,] %>% reshape2::melt()
  dnaMG <- RRBS_completebetas[rownames(RRBS_completebetas) == probe,] %>% reshape2::melt()
  
  DiffExp <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/C_RNASeq/1_RNASeq_Isabel/Isabel_Supp2_Tg4510GenotypeDEG.csv")
  DiffExp[DiffExp$Gene == targetGene,]
  DiffExp[DiffExp$Gene == tranFactor,]
  
  
  merge(expG, coldata, by.x = "variable", by.y = 0) %>% ggplot(., aes(x = Genotype, y = value)) + geom_boxplot()
  p1 <- merge(expTF, coldata, by.x = "variable", by.y = 0) %>% ggplot(., aes(x = Genotype, y = value, fill = Genotype)) + geom_boxplot() +
    labs(x = "Genotype", y= "Normalised counts", title = paste0(targetGene," gene expression")) + theme_classic() +
    scale_fill_manual(values = c("grey","red")) + theme(legend.position = "None")
  
  p2 <- merge(dnaMG, coldata, by.x = "variable", by.y = 0) %>% ggplot(., aes(x = Genotype, y = value, fill = Genotype)) + geom_boxplot() +
    labs(x = "Genotype", y= "Methylation", title = paste0(probe,"")) + theme_classic() +
    scale_fill_manual(values = c("grey","red")) + theme(legend.position = "None")

  p3 <- merge(merge(expG, coldata, by.x = "variable", by.y = 0) %>% 
                dplyr::rename(exp = value), dnaMG, by = "variable") %>% dplyr::rename(meth = value) %>% 
    ggplot(., aes(x = meth, y = exp, colour = Genotype)) + geom_point() + theme_classic() + 
    labs(x = paste0(substitute(probe), " methylation"),
         y = paste0("Target ", substitute(targetGene), " expression")) +
      scale_colour_manual(values = c("grey","red")) + theme(legend.position = "None")
  
  
  p4 <- merge(merge(expG %>% rename(targetGene = value), expTF %>% rename(tranFactor = value), by = "variable"), coldata, by.x = "variable", by.y = 0) %>% 
    ggplot(., aes(x = tranFactor, y = targetGene, colour = Genotype)) + geom_point() + theme_classic() + 
    labs(x = paste0("TF ", substitute(tranFactor), " expression"),
         y = paste0("Target ", substitute(targetGene), " expression")) +
    scale_colour_manual(values = c("grey","red")) + theme(legend.position = "None")
  

  return(plot_grid(p1,p2,p3,p4))
}

plot_integration_raw("Ninj1","Gata3","chr13:49341599")
plot_integration_raw("Cul4b","Foxj3","chrX:37666823")
plot_integration_raw("Piwil1","Stat2","chr5:128780460")
