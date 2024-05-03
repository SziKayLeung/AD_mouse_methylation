# Reference: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.annotation.gtf.gz

suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# chipseeker annotation
source("/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/4_DMP_annotation/0_source_functions.R")
# data wrangle for beta output, phenotype file output, and mixedeffectsmodels results
source("/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rrbs-ad-mice/5_Array/5a_MEM_analysisFunctions.R")

output_dir = "/lustre/projects/Research_Project-MRC148213/lsl693/rrbs_ad_mice/0_ZenOutput/2_annotated/"
dirnames <- list(
  biseq = "/gpfs/ts0/projects/Research_Project-191406/isabel/RRBS_new/",
  AishaArray = "/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/",
  array_anno = "/lustre/projects/Research_Project-191406/EmmaW/Array/Results",
  rrbs_anno = "/lustre/projects/Research_Project-191406/EmmaW/RRBSAnnotatedResults"
)

## RRBS
# results after applying BiSeq and smoothing
# load from two mouse models and save into list
RRBS_complete <- list()
load(file = paste0(dirnames$biseq, "rTg4510_RRBSbetasComplete.RData"))
RRBS_complete$rTg4510 <- RRBS_completebetas
load(file = paste0(dirnames$biseq, "J20_RRBSbetasComplete.RData"))
RRBS_complete$J20 <- RRBS_completebetas
RRBS_complete <- lapply(RRBS_complete, function(x) x %>% tibble::rownames_to_column(., var = "position"))


# array files
mm10_Manifest = read.csv(paste0(dirnames$AishaArray, "mm10_Manifest.csv"), stringsAsFactors = F, header = T, row.names = 1)

raw <- list(
  array = list(
    J20Array= read.csv(paste0(dirnames$AishaArray,"J20_coldata_VertebrateArray.csv"), header = T, stringsAsFactors = F),
    rTg4510Array = read.csv(paste0(dirnames$AishaArray,"Tg4510_coldata_VertebrateArray.csv"),header = T, stringsAsFactors = F)
  )
)

interaction <- list(
  array = list(
    rTg4510_ECX_geno = paste0(dirnames$array_anno,"/rTg4510/MixedModelResultsECX_rTg4510.csv"),
    rTg4510_HIP_geno = paste0(dirnames$array_anno,"/rTg4510/MixedModelResultsHIP_rTg4510.csv"),
    J20_ECX_geno = paste0(dirnames$array_anno,"/J20/MixedModelResultsECX_J20.csv"),
    J20_HIP_geno = paste0(dirnames$array_anno,"/J20/MixedModelResultsHIP_J20.csv"))
)
interaction$array <- lapply(interaction$array, function(x) read.table(x,header = T, sep = ",", row.names = "X"))

anno <- list(
  array = list(
    rTg4510_ECX_geno = PrepMixedEffectsStats(interaction$array$rTg4510_ECX_geno, mm10_Manifest, intTerm="Geno.Age"),
    J20_ECX_geno = PrepMixedEffectsStats(interaction$array$J20_ECX_geno, mm10_Manifest, intTerm="Geno.Age"),
    rTg4510_HIP_geno = PrepMixedEffectsStats(interaction$array$rTg4510_HIP_geno, mm10_Manifest, intTerm="Geno.Age"),
    J20_HIP_geno = PrepMixedEffectsStats(interaction$array$J20_HIP_geno, mm10_Manifest, intTerm="Geno.Age"),
    rTg4510_ECX_interaction = PrepMixedEffectsStats(interaction$array$rTg4510_ECX_geno, mm10_Manifest, intTerm="Geno.Age")
  ),
  rrbs = list(
    rTg4510_ECX_geno = paste0(output_dir,"rTg4510_rrbs_genotype.csv"),
    J20_ECX_geno = paste0(output_dir,"J20_rrbs_genotype.csv"),
    rTg4510_ECX_interaction = paste0(output_dir,"rTg4510_rrbs_interaction.csv"),
    rTg4510_ECX_pathology = paste0(output_dir,"rTg4510_rrbs_pathology.csv"),
    J20_ECX_interaction = paste0(output_dir,"J20_rrbs_interaction.csv"),
    J20_ECX_pathology = paste0(output_dir,"J20_rrbs_pathology.csv"))
)
anno$rrbs <- lapply(anno$rrbs, function(x) read.csv(x, header = T))

# phenotype data
phenotype <- list(
  rTg4510 = read.csv("/lustre/projects/Research_Project-191406/isabel/RRBS_new/Tg4510_coldata_RRBS.csv", stringsAsFactors=FALSE),
  J20 = read.csv("/lustre/projects/Research_Project-191406/isabel/RRBS_new/J20_coldata_RRBS.csv", stringsAsFactors=FALSE)
)
phenotype <- lapply(phenotype, function(x) x %>% mutate(Age_months = factor(Age_months), Genotype = factor(Genotype, levels = c("WT","TG"))))
