## ---------- Script -----------------
##
## Script name: run_methreg_unsupervised.R
##
## Purpose of script: run methreg analysis (promoter) through rTg4510 methylation and gene expression datasets
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------


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
  make_option(c("-c", "--chr"), type="character", default=NULL, help="chromosome", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL,help="output name", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, help="choice of platform: rrbs, array", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$name)){
  stop("need name argument, --name")
}

if(!opt$type %in% c("rrbs","array")){
  stop("type argument needs to be either rrbs or array, --type <array,rrbs>")
}

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

message("Output files deposited to:", dirnames$output)

## ---------- Data input ---------------------------

# dna methylation
if(opt$type == "rrbs"){
  message("Read in DNA methylation from RRBS data")
  load(file = paste0(dirnames$biseq, "rTg4510_RRBSbetasComplete.RData"))
  #sigResults <- list(rTg4510.genotype = read.csv(paste0(dirnames$sigResults,"rTg4510_rrbs_genotype.csv")) %>% filter(FDR_adj_genotype < 0.05))
  #dnam <- RRBS_completebetas[rownames(RRBS_completebetas) %in% sigResults$rTg4510.genotype$Position, ]
  
  message("Subsetting methylation data to: ",opt$chr)
  dnam <- RRBS_completebetas[which(grepl(opt$chr, rownames(RRBS_completebetas))),]
  message("Number of sites: ", nrow(dnam))
  rownames(dnam) <- sub("(.*):(\\d+)", "\\1:\\2-\\2", rownames(dnam))
}else{
  message("Read in DNA methylation from array data")
  dnam = read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rrbs-ad-mice/Figures/rTg4510_ECX_array_beta.csv")
  manifest = read.csv("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/mm10_Manifest.csv", stringsAsFactors = F, header = T, row.names = 1)
  manifest$position = paste0(manifest$seqnames,":",manifest$start,"-",manifest$end)
  dnam <- merge(dnam, manifest[,c("cpg","position")], by = 0)
  dnam <- dnam[!duplicated(dnam$position),]
  
  message("Subsetting methylation data to: ",opt$chr)
  dnam$chr <- word(dnam$position,c(1), sep = fixed(":"))
  dnam <- dnam[dnam$chr == opt$chr,]
  message("Number of sites: ", nrow(dnam))
  rownames(dnam) <- dnam$position
}


# gene expression
message("Read in gene expression")
gene.exp <- read.csv(paste0(dirnames$isabel, "rTg4510_DESeq2normalizedcounts.csv"))
gene.exp$ens <- map_symbol_to_ensg(gene.exp$X , genome=genome)
gene.exp <- gene.exp[gene.exp$ens != "NA",]
rownames(gene.exp) <- gene.exp$ens
gene.exp <- gene.exp %>% dplyr::select(-c(X,ens))

# match samples
common.samples <- intersect(colnames(dnam),colnames(gene.exp))
dnam <- dnam[common.samples]
exp <- gene.exp[common.samples]
if(!all(colnames(dnam) == colnames(exp))){
  print("Missing samples in gene expression and dna methylation")
  system.exit()
}


## ---------- Mapping regions ----------------------

regions <- rownames(dnam)
regions.gr <- make_granges_from_names(regions)
regions.names <- make_names_from_granges(regions.gr)

# paramaters
target.method = "window"
target.window.size = 500 * 10^3
target.num.flanking.genes = 5
target.promoter.upstream.dist.tss = 2000
target.promoter.downstream.dist.tss = 2000
target.rm.promoter.regions.from.distal.linking = TRUE
motif.search.window.size = 50
motif.search.p.cutoff = 10^-3
TF.peaks.gr = NULL
max.distance.region.target =  10^6
genome = "mm10"
window.size = 500 * 10^3
cores = 1

# gene target genes
# use the window size rather than the promoter overlap as not as specific and missing target genes
message("Get target region from methylation")
if(target.method == "promoter_overlap"){
  message("Using promoter overlap")
  region.target = get_region_target_gene_by_promoter_overlap(regions.gr)
}else{
  message("Using window size 500kbp")
  region.target = get_region_target_gene_window(regions.gr, genome, window.size)
}

# get transcript factor
message("Get transcription factor")
region.tf <- get_tf_in_region(region = regions.gr,
                              genome = genome,
                              window.size = motif.search.window.size,
                              p.cutoff = motif.search.p.cutoff,
                              cores = cores,
                              TF.peaks.gr = TF.peaks.gr
)

write.table(region.tf,paste0(dirnames$output,opt$name,"_region.tf"),quote=F)


triplet <- dplyr::inner_join(region.target, region.tf,by = join_by(regionID),relationship = "many-to-many")
triplet <- triplet %>% dplyr::filter(!is.na(.data$TF))

if(target.method == "promoter_overlap"){
  message("Removing regions and target genes from different chromosomes")
  triplet <- triplet %>% dplyr::filter(!is.na(.data$distance_region_target_tss))
  
  message("Removing regions and target genes with ditance higher than ", max.distance.region.target, " bp")
  triplet <- triplet %>% dplyr::filter(abs(.data$distance_region_target_tss) < max.distance.region.target)
  
  triplet <- triplet %>%
    dplyr::relocate(.data$distance_region_target_tss, .after = dplyr::last_col()) %>%
    dplyr::relocate(contains("pos"), .after = dplyr::last_col())
}


write.csv(triplet,paste0(dirnames$output,"triplet_",opt$name,"_",opt$chr,".csv"),quote=F,row.names = F)
#triplet <- read.csv(paste0(dirnames$script,"triplet.csv"))


## ---------- Run interaction model ----------------

tf.activity.es = NULL
dnam.group.threshold = 0.25
cores=1
sig.threshold = 0.05

gene.exp <- filter_genes_zero_expression_all_samples(gene.exp)
regions.keep <- (rowSums(is.na(dnam)) < (ncol(dnam) * 0.75)) %>% which %>% names
dnam <- dnam[regions.keep,,drop = FALSE]
triplet <- triplet %>% dplyr::filter(
  .data$target %in% rownames(gene.exp) &
    .data$regionID %in% rownames(dnam)
)

triplet <- triplet %>% dplyr::filter(
  .data$TF %in% rownames(gene.exp)
)

# Remove cases where target is also the TF if it exists
triplet <- triplet %>% dplyr::filter(
  as.character(.data$TF) != as.character(.data$target)
)

if(!"TF_symbol" %in% colnames(triplet))
  triplet$TF_symbol <- map_ensg_to_symbol(triplet$TF)

if(!"target_symbol" %in% colnames(triplet))
  triplet$target_symbol <- map_ensg_to_symbol(triplet$target)

if(!"target_region" %in% colnames(triplet))
  triplet$target_region <- map_ensg_to_region(triplet$target)

if(!"distance_region_target_tss" %in% colnames(triplet)){
  triplet$distance_region_target_tss <- get_target_tss_to_region_distance(triplet$regionID,triplet$target)
}


parallel <- register_cores(cores)


ret <- plyr::adply(
  .data = triplet,
  .margins = 1,
  .fun = function(
  row.triplet
  ){
    
    data <- get_triplet_data(
      exp = exp,
      dnam = dnam,
      row.triplet = row.triplet,
      tf.es = tf.activity.es
    )
    
    
    upper.cutoff <-  quantile(data$met,na.rm = TRUE,  1 - dnam.group.threshold)
    low.cutoff <-  quantile(data$met,na.rm = TRUE,  dnam.group.threshold)
    
    quant.diff <- data.frame("met.IQR" = upper.cutoff - low.cutoff)
    
    data.high.low <- data %>% filter(.data$met <= low.cutoff | .data$met >= upper.cutoff)
    data.high.low$metGrp <- ifelse(data.high.low$met <= low.cutoff, 0, 1)
    
    # pct.zeros.in.samples <- sum(data$rna.target == 0, na.rm = TRUE) / nrow(data)
    
    if(nrow(data.high.low %>% filter(.data$metGrp == 1)) != 0){
      
      suppressWarnings({
        # Add information to filter TF if differenly expressed between DNAm high and DNAm low groups
        wilcoxon.tf.q4met.vs.q1met <- wilcox.test(
          data.high.low %>% filter(.data$metGrp == 1) %>% pull(.data$rna.tf),
          data.high.low %>% filter(.data$metGrp == 0) %>% pull(.data$rna.tf),
          exact = FALSE
        )$p.value
      })
      
      suppressWarnings({
        # Add information to filter Target if differently expressed between DNAm high and DNAm low groups
        wilcoxon.target.q4met.vs.q1met <- wilcox.test(
          data.high.low %>% filter(.data$metGrp == 1) %>% pull(.data$rna.target),
          data.high.low %>% filter(.data$metGrp == 0) %>% pull(.data$rna.target),
          exact = FALSE
        )$p.value
      })
      
      pct.zeros.in.quant.samples <- sum(
        data.high.low$rna.target == 0,
        na.rm = TRUE) / nrow(data.high.low)
      
      if (pct.zeros.in.quant.samples > 0.25) {
        itx.quant <- interaction_model_quant_zeroinfl(data.high.low)
      } else {
        itx.quant <- interaction_model_quant_rlm(data.high.low)
      }
      
      # Create output
      interaction_model_output(
        quant.diff,
        itx.quant,
        pct.zeros.in.quant.samples,
        wilcoxon.target.q4met.vs.q1met = wilcoxon.target.q4met.vs.q1met,
        wilcoxon.tf.q4met.vs.q1met = wilcoxon.tf.q4met.vs.q1met
      )
    }
   
  },
  .progress = "time",
  .parallel = parallel,
  .inform = TRUE,
  .paropts = list(.errorhandling = 'pass')
)

write.csv(ret,paste0(dirnames$output,"ret_",opt$name,"_",opt$chr,".csv"),quote=F,row.names = F)
#ret <- read.csv(paste0(dirnames$script,"ret.csv"))


## ---------- Stage wise correction ----------------

perform_correction <- function(name, chr){
  
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
  write.csv(ret,paste0(dirnames$output,"retCorrected_", name,"_", chr,".csv"),quote=F,row.names = F)
}

perform_correction(opt$name, opt$chr)

## ---------- Output ----------------

plot_interaction <- function(name, chr, gene.exp){
  ret <- read.csv(paste0("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/C_RNASeq/3_integration/retCorrected_",name,"_", chr,".csv"))
  
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
  
  pdf(paste0(dirnames$output,"MethRegOutput_",name,"_", chr,".pdf"), width = 15, height = 10)
  for(i in retSigTarget$regionID){
    print(plot_interaction_model(triplet.results = ret[ret$regionID == i,], dnam = as.matrix(dnam), exp = as.matrix(exp)))
  }
  dev.off()
}

plot_interaction(opt$name, opt$chr, gene.exp)