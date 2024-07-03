
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "2_Array/2_rTg4510_DMP/import.config"))
source(paste0(scriptDir, "2_Array/functions/finalise_betaRegressionArray.R"))

# save as Rdata
ResultDir <- "/gpfs/ts0/projects/Research_Project-191406/EmmaW/Array/Results/rTg4510/"

rTg4510_array_ECX_DMP <- list(
  Genotype = read.csv(paste0(ResultDir, "MixedModelResultsECX_rTg4510.csv")),
  Interaction = read.csv(paste0(ResultDir, "MixedModelResultsECX_rTg4510.csv")),
  Pathology = read.csv(paste0(ResultDir, "MixedModelResultsECX_Patholgy_rTg4510.csv"))
)

rTg4510_array_HIP_DMP <- list(
  Genotype = read.csv(paste0(ResultDir, "MixedModelResultsHIP_rTg4510.csv")),
  Interaction = read.csv(paste0(ResultDir, "MixedModelResultsHIP_rTg4510.csv")),
  Pathology = read.csv(paste0(ResultDir, "MixedModelResultsHIP_Pathology_rTg4510.csv"))
)

rTg4510_array_tissue_DMP <- list(
  Genotype = read.csv(paste0(ResultDir, "MixedModelResults_Tissue_rTg4510.csv"))
)

save(rTg4510_array_ECX_DMP, file = paste0(dirnames$differential, "/array/rTg4510_array_ECX_allResultsDMPs.RData"))
save(rTg4510_array_HIP_DMP, file = paste0(dirnames$differential, "/array/rTg4510_array_HIP_allResultsDMPs.RData"))
save(rTg4510_array_tissue_DMP, file = paste0(dirnames$differential, "/array/rTg4510_array_tissue_allResultsDMPs.RData"))

### significance

rTg4510_array_ECX_DMP <- list(
  Genotype = PrepMixedEffectsStats(mixedEffectsResults=rTg4510_array_ECX_DMP$Genotype, manifest=mm10_Manifest, test="Geno"),
  Interaction = PrepMixedEffectsStats(mixedEffectsResults=rTg4510_array_ECX_DMP$Genotype, manifest=mm10_Manifest, test="Geno.Age"),
  Pathology = PrepMixedEffectsStats(mixedEffectsResults=rTg4510_array_ECX_DMP$Pathology, manifest=mm10_Manifest, test="Pathology")
)

rTg4510_array_HIP_DMP <- list(
  Genotype = PrepMixedEffectsStats(mixedEffectsResults=rTg4510_array_HIP_DMP$Genotype, manifest=mm10_Manifest, test="Geno"),
  Interaction = PrepMixedEffectsStats(mixedEffectsResults=rTg4510_array_HIP_DMP$Genotype, manifest=mm10_Manifest, test="Geno.Age"),
  Pathology = PrepMixedEffectsStats(mixedEffectsResults=rTg4510_array_HIP_DMP$Pathology, manifest=mm10_Manifest, test="Pathology")
)

rTg4510_array_tissue_DMP <- list(
  Genotype = PrepMixedEffectsStats(mixedEffectsResults=rTg4510_array_tissue_DMP$Genotype, manifest=mm10_Manifest, test="Tissue.Genotype")
)

save(rTg4510_array_ECX_DMP, file = paste0(dirnames$differential, "/array/rTg4510_array_ECX_sigResultsDMPs.RData"))
save(rTg4510_array_HIP_DMP, file = paste0(dirnames$differential, "/array/rTg4510_array_HIP_sigResultsDMPs.RData"))
save(rTg4510_array_tissue_DMP, file = paste0(dirnames$differential, "/array/rTg4510_array_tissue_sigResultsDMPs.RData"))

