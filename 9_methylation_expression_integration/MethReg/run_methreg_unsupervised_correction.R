
## ---------- packages -----------------

suppressMessages({
  library("tidyr")
  library("dplyr")
  library("stringr")
  library("GenomicRanges")
  library("JASPAR2022")
  library("TFBSTools")
  library("motifmatchr")
  library("parallel")
  library("doParallel")
  library("BSgenome.Mmusculus.UCSC.mm10")
  library("stageR")
  library("MASS")
  library("pscl")
  library("tibble")
  library("ggplot2")
  library("ggpubr")
  library("MethReg")
  library("optparse")
})


## ---------- arguments -----------------
option_list <- list( 
  make_option(c("-c", "--chr"), type="character", default=NULL, 
              help="chromosome", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


## ---------- directory and scripts -----------------

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
gene.exp$ens <- map_symbol_to_ensg(gene.exp$X , genome=genome)
gene.exp <- gene.exp[gene.exp$ens != "NA",]
rownames(gene.exp) <- gene.exp$ens
gene.exp <- gene.exp %>% dplyr::select(-c(X,ens))


## ---------- functions -----------------

perform_correction <- function(chr){
  
  message("Input: /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/C_RNASeq/3_integration/ret_",chr,".csv")
  inputRet <- read.csv(paste0("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/C_RNASeq/3_integration/ret_",chr,".csv"))
  
  message("Performing Stage wise correction for triplets")
  ret <- calculate_stage_wise_adjustment(inputRet)
  
  message("Filtering results to have interaction, TF or DNAm significant")
  sig.threshold = 0.05
  ret <- ret %>% filter_at(vars(
    intersect(
      contains("triplet_stage_wise_adj_pvalue", ignore.case = TRUE),
      contains("RLM", ignore.case = TRUE)
    )
  ), 
  any_vars(. < sig.threshold)
  )
  
  message("Filtering results to remove the significant in the wilcoxon test TF Q1 vs Q4")  
  ret <- ret %>% dplyr::filter(.data$TF_DNAm_high_vs_TF_DNAm_low_wilcoxon_pvalue > sig.threshold)
  
  message("Writing output")
  write.csv(ret,paste0(dirnames$output,"retCorrected_",chr,".csv"),quote=F,row.names = F)
}

plot_interaction <- function(chr, gene.exp){
  ret <- read.csv(paste0("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/C_RNASeq/3_integration/retCorrected_", chr,".csv"))
  
  message("Subsetting methylation data to: ",chr)
  dnam <- RRBS_completebetas[which(grepl(opt$chr, rownames(RRBS_completebetas))),]
  message("Number of sites: ", nrow(dnam))
  rownames(dnam) <- sub("(.*):(\\d+)", "\\1:\\2-\\2", rownames(dnam))
  
  common.samples <- intersect(colnames(dnam),colnames(gene.exp))
  dnam <- dnam[common.samples]
  exp <- gene.exp[common.samples]
  if(!all(colnames(dnam) == colnames(exp))){
    print("Missing samples in gene expression and dna methylation")
    system.exit()
  }
  
  retSigTarget <- ret %>% filter(Target_gene_DNAm_high_vs_Target_gene_DNAm_low_wilcoxon_pvalue < 0.1) %>% 
    arrange(Target_gene_DNAm_high_vs_Target_gene_DNAm_low_wilcoxon_pvalue)
  
  pdf(paste0(dirnames$output,"MethRegOutput_",chr,".pdf"), width = 15, height = 10)
  for(i in retSigTarget$regionID){
    print(plot_interaction_model(triplet.results = ret[ret$regionID == i,], dnam = as.matrix(dnam), exp = as.matrix(exp)))
  }
  dev.off()
}


#perform_correction(opt$chr)
plot_interaction(opt$chr, gene.exp)