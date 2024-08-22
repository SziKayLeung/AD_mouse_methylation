
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "2_Array/functions/finalise_betaRegressionArray.R"))

# save as Rdata
ResultDir <- "/gpfs/ts0/projects/Research_Project-191406/EmmaW/Array/Results/J20/"

J20_array_ECX_DMP <- list(
  Genotype = read.csv(paste0(ResultDir, "MixedModelResultsECX_J20.csv")),
  Interaction = read.csv(paste0(ResultDir, "MixedModelResultsECX_J20.csv")),
  Pathology = read.csv(paste0(ResultDir, "MixedModelResultsECX_Pathology_J20.csv"))
)

J20_array_HIP_DMP <- list(
  Genotype = read.csv(paste0(ResultDir, "MixedModelResultsHIP_J20.csv")),
  Interaction = read.csv(paste0(ResultDir, "MixedModelResultsHIP_J20.csv")),
  Pathology = read.csv(paste0(ResultDir, "MixedModelResultsHIP_Pathology_J20.csv"))
)

J20_array_tissue_DMP <- list(
  Genotype = read.csv(paste0(ResultDir, "MixedModelResults_Tissue_J20.csv"))
)

save(J20_array_ECX_DMP, file = paste0(dirnames$differential, "/array/J20_array_ECX_allResultsDMPs.RData"))
save(J20_array_HIP_DMP, file = paste0(dirnames$differential, "/array/J20_array_HIP_allResultsDMPs.RData"))
save(J20_array_tissue_DMP, file = paste0(dirnames$differential, "/array/J20_array_tissue_allResultsDMPs.RData"))

### significance

J20_array_ECX_DMP <- list(
  Genotype = PrepMixedEffectsStats(mixedEffectsResults=J20_array_ECX_DMP$Genotype, manifest=mm10_Manifest, test="Geno"),
  Interaction = PrepMixedEffectsStats(mixedEffectsResults=J20_array_ECX_DMP$Genotype, manifest=mm10_Manifest, test="Geno.Age"),
  Pathology = PrepMixedEffectsStats(mixedEffectsResults=J20_array_ECX_DMP$Pathology, manifest=mm10_Manifest, test="Pathology")
)

J20_array_HIP_DMP <- list(
  Genotype = PrepMixedEffectsStats(mixedEffectsResults=J20_array_HIP_DMP$Genotype, manifest=mm10_Manifest, test="Geno"),
  Interaction = PrepMixedEffectsStats(mixedEffectsResults=J20_array_HIP_DMP$Genotype, manifest=mm10_Manifest, test="Geno.Age"),
  Pathology = PrepMixedEffectsStats(mixedEffectsResults=J20_array_HIP_DMP$Pathology, manifest=mm10_Manifest, test="Pathology")
)

J20_array_tissue_DMP <- list(
  Genotype = PrepMixedEffectsStats(mixedEffectsResults=J20_array_tissue_DMP$Genotype, manifest=mm10_Manifest, test="Tissue.Genotype")
)

save(J20_array_ECX_DMP, file = paste0(dirnames$differential, "/array/J20_array_ECX_sigResultsDMPs.RData"))
save(J20_array_HIP_DMP, file = paste0(dirnames$differential, "/array/J20_array_HIP_sigResultsDMPs.RData"))
save(J20_array_tissue_DMP, file = paste0(dirnames$differential, "/array/J20_array_tissue_sigResultsDMPs.RData"))

