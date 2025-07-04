library("dplyr")


dirnames <- list(
  script = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/AD_mouse_methylation/",
  metadata = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/0_metadata",
  processed = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/1_processed",
  differential = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/2_differentialAnalysis",
  annotated = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/3_annotated",
  clock = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/5_epigeneticClock",
  
  utils = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/0_utils/",
  
  # bismark
  rTg4510bismark = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/RRBS/rTg4510/3_bismark/",
  J20bismark = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/RRBS/J20/3_bismark/",
  
  # pyrosequencing output
  pyroRaw = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/0_pyro/",
  
  # human annotations 
  humanAnnot = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/4_humanAnnotations/",
  
  # GO annotations list 
  GO = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/7_GO/",
  
  # Paper output
  paper = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/AD_Mouse_Model/rTg4510_mice_methylation_paper/0_ZenOutput/PaperOutput"
)

datawrangle_phenotype <- function(path_pheno){
  
  # Phenotypic data
  coldata <- read.csv(path_pheno, row.names=1, stringsAsFactors=FALSE)
  coldata$Age_months <- as.numeric(coldata$Age_months)
  coldata$Genotype <- as.factor(coldata$Genotype)
  coldata$Genotype <- relevel(coldata$Genotype, "WT")
  
  
  return(coldata)
}

phenotype <- lapply(list.files(path = dirnames$metadata, pattern = "coldata_RRBS", full.names = T), function(x) datawrangle_phenotype(x))
phenotype <- lapply(phenotype, function(x) x %>% dplyr::select(Genotype, Age_months, Histology_no, HIP, ECX))
names(phenotype) <- c("J20", "rTg4510")

# full phenotype with pathology data
phenotype_path <- list(
  J20 = read.csv(paste0(dirnames$metadata, "/J20_coldata_VertebrateArray.csv"), header = T, stringsAsFactors = F),
  rTg4510 = read.csv(paste0(dirnames$metadata, "/Tg4510_coldata_VertebrateArray.csv"), header = T, stringsAsFactors = F)
)

phenotype$J20_HIP <- phenotype_path$J20 %>% tibble::column_to_rownames(., var = "Sample_ID_HIP") %>% 
  dplyr::select("Genotype", "Age_months","Histology_no","Pathology_HIP", "Pathology_ECX") %>%  `colnames<-`(colnames(phenotype$J20))
phenotype$rTg4510_HIP <- phenotype_path$rTg4510 %>% tibble::column_to_rownames(., var = "Sample_ID_HIP") %>% 
  dplyr::select("Genotype", "Age_months","Histology_no","Pathology_HIP", "Pathology_ECX") %>%  `colnames<-`(colnames(phenotype$rTg4510)) %>% 
  mutate(Genotype = factor(Genotype, levels = c("WT","TG")))

phenotype$rTg4510_all <- rbind(phenotype$rTg4510, phenotype$rTg4510_HIP)
phenotype$J20_all <- rbind(phenotype$J20, phenotype$J20_HIP)

## mm10 mannifest file (created from below by A.Dahir)
#mm10_Manifest <- read.csv(paste0(dirnames$utils, "mm10_Manifest.csv"), stringsAsFactors = F, header = T, row.names = 1)
#mm10_Manifest$position <- paste0(mm10_Manifest$seqnames,":",mm10_Manifest$start)
#mm10_Manifest <- mm10_Manifest[,c("position","Gene_Symbol","annotation")]
#write.csv(mm10_Manifest, paste0(dirnames$utils, "mm10_Final_Manifest.csv"))
mm10_Manifest <- read.csv(paste0(dirnames$utils, "mm10_Final_Manifest.csv"), stringsAsFactors = F, header = T, row.names = 1)

## Horvath mammalian chip coordinates
HorvathMammalChip <- data.table::fread(paste0(dirnames$utils, "HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.1.csv"), data.table = F)
##filter for mouse samples
#species_mappability_name = "MusMusculus"
#species_probes <- HorvathMammalChip[, c('probeID','MusMusculus')] 
#species_probes <- species_probes[species_probes$MusMusculus != "",]
#nrow(species_probes)
#species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
#species_probes$position <- paste0("chr", species_probes$MusMusculus)
#head(species_probes)
#write.csv(species_probes, paste0(dirnames$utils,"HorvathMammalChip_mm10GenomeCoordsUnique_v1.1.csv"), row.names = F, quote = F)

# Horvath epigenetic Clock
OutputMouseClocks <- read.csv(paste0(dirnames$metadata, "/HorvathArray/OutputMouseClocks.csv"), header = T,stringsAsFactors = F)

## pyrosequencing
# pyrosequencing data: Ank1
input_pyro <- list(
  bin1 = read.csv(paste0(dirnames$pyroRaw, "Toni_bin1_failed.csv")),
  prnp = read.csv(paste0(dirnames$pyroRaw, "prnp_dmr.csv")),
  ank1 = read.csv(paste0(dirnames$pyroRaw, "ank1_dmp.csv"))
)

## human annotations
# read in gene conversion file
mouse2human <- read.csv(paste0(dirnames$script, "4_HumanGeneListComparisons/mousehumangeneconversion.csv"))

# read in gene lists from G.Shireby and R.Smith lists (EWAS)
inputHumanListDir <- paste0(dirnames$script, "4_HumanGeneListComparisons/inputLists")
inputHumanLists <- list.files(path = inputHumanListDir, pattern = "txt", full.names = T)
inputHumanLists <- lapply(inputHumanLists, function(x) read.table(x,col.names=c("human")))
names(inputHumanLists) <- list.files(path = inputHumanListDir, pattern = "txt")

# number of raw reads RRBS
rawRRBSReads <- list(
  rTg4510 = read.csv(paste0(dirnames$processed,"/rrbs/Tg4510NumberRawReads.csv")),
  J20 = read.csv(paste0(dirnames$processed,"/rrbs/J20NumberRawReads.csv"))
)
