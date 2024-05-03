LOGEN_ROOT="/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen"
source(paste0(LOGEN_ROOT, "/aesthetics_basics_plots/draw_density.R"))
source(paste0(LOGEN_ROOT, "/aesthetics_basics_plots/pthemes.R"))


probe_mapping_file <- data.table::fread("/gpfs/ts0/projects/Research_Project-191406/Horvath/MammalianArrayNormalizationTools/mappabilityData/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.1.csv", data.table = F)

#filter for mouse samples
species_mappability_name = "MusMusculus"
species_probes <- probe_mapping_file[, c('probeID','MusMusculus')] 
species_probes <- species_probes[species_probes$MusMusculus != "",]
nrow(species_probes)
species_probes$Chr <- gsub(":.*","",species_probes$MusMusculus)
species_probes$position <- paste0("chr", species_probes$MusMusculus)
head(species_probes)

length(RRBS_complete$rTg4510$position)
length(species_probes$position)
length(intersect(RRBS_complete$rTg4510$position, species_probes$position))

length(intersect(RRBS_complete$J20$position, species_probes$position))

load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame.rdat")


pheno <- QCmetrics
betas <- Normalised_Sesame_Betas

pheno_ecx <- pheno[which(pheno$Tissue %in% "CortexEntorhinalis"),]
betas_ecx <- betas[, colnames(betas) %in% pheno_ecx$Basename]
pheno_ecx_J20 <- pheno_ecx[which(pheno_ecx$AD_model %in% "J20"),] 
betas_ecx_J20 <- betas_ecx[, colnames(betas_ecx) %in% pheno_ecx_J20$Basename] %>% tibble::rownames_to_column(., var = "probeID")
pheno_ecx_Tg4510 <- pheno_ecx[which(pheno_ecx$AD_model %in% "rTg4510"),] 
betas_ecx_Tg4510 <- betas_ecx[, colnames(betas_ecx) %in% pheno_ecx_Tg4510$Basename] %>% tibble::rownames_to_column(., var = "probeID")

length(betas_ecx_Tg4510$position)
betas_ecx_Tg4510 <- merge(betas_ecx_Tg4510, species_probes, by = "probeID", all.x = T)
length(intersect(RRBS_complete$rTg4510$position, betas_ecx_Tg4510$position))
betas_ecx_J20 <- merge(betas_ecx_J20, species_probes, by = "probeID", all.x = T)
length(betas_ecx_J20$position)
length(intersect(RRBS_complete$J20$position, betas_ecx_J20$position))

# J20 
commonJ20 <- intersect(RRBS_complete$J20$position, betas_ecx_J20$position)
length(commonJ20)
length(intersect(RRBS_complete$J20$position, betas_ecx_J20$position))/length(betas_ecx_J20$position) * 100
length(intersect(RRBS_complete$J20$position, betas_ecx_J20$position))/length(RRBS_complete$J20$position) * 100

# rTg4510 
commonrTg4510 <- intersect(RRBS_complete$rTg4510$position, betas_ecx_Tg4510$position)
length(commonrTg4510)
length(intersect(RRBS_complete$rTg4510$position, betas_ecx_Tg4510$position))/length(betas_ecx_Tg4510$position) * 100
length(intersect(RRBS_complete$rTg4510$position, betas_ecx_Tg4510$position))/length(RRBS_complete$rTg4510$position) * 100

RRBS_complete$J20 <- replace(RRBS_complete$J20, is.na(RRBS_complete$J20), 0)
commonJ20RRBS <- RRBS_complete$J20[RRBS_complete$J20$position %in% commonJ20,] 
rownames(commonJ20RRBS) <- commonJ20RRBS$position

# note there are two probes for array: chr15:82917534, chr15:8444274 with same position (2 different positions/examples)
# took the representative probe with higher methylation
commonJ20Array <- betas_ecx_J20[betas_ecx_J20$position %in% commonJ20,]
commonJ20Array <- commonJ20Array[!commonJ20Array$probeID %in% c("cg12544505.2","cg01079397.1"),]
rownames(commonJ20Array) <- commonJ20Array$position

commonJ20RRBSdf <- commonJ20RRBS %>% dplyr::select(-position) %>% apply(., 1, mean) %>% reshape2::melt(value.name = "RRBS")
commonJ20Arraydf <- commonJ20Array %>% dplyr::select(-position, -probeID, -Chr, -MusMusculus) %>% apply(., 1, mean) %>% reshape2::melt(value.name = "Array")


density_plot(merge(commonJ20Arraydf, commonJ20RRBSdf, by = 0),x.var = "RRBS", y.var = "Array", x_lab = "RRBS", y_lab = "Array", title = "") + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10")

# rTg4510
RRBS_complete$rTg4510 <- replace(RRBS_complete$rTg4510, is.na(RRBS_complete$rTg4510), 0)
commonTg4510RRBS <- RRBS_complete$rTg4510[RRBS_complete$rTg4510$position %in% commonrTg4510,] 
rownames(commonTg4510RRBS) <- commonTg4510RRBS$position

# note there are two probes for array: chr15:8444274, chr15:82917534 with same position (2 different positions/examples)
# took the representative probe with higher methylation
commonTg4510Array <- betas_ecx_Tg4510[betas_ecx_Tg4510$position %in% commonrTg4510,]
commonTg4510Array <- commonTg4510Array[!commonTg4510Array$probeID %in% c("cg01079397.2","cg12544505.1"),]
rownames(commonTg4510Array) <- commonTg4510Array$position

commonTg4510RRBSdf <- commonTg4510RRBS %>% dplyr::select(-position) %>% apply(., 1, mean) %>% reshape2::melt(value.name = "RRBS")
commonTg4510Arraydf <- commonTg4510Array %>% dplyr::select(-position, -probeID, -Chr, -MusMusculus) %>% apply(., 1, mean) %>% reshape2::melt(value.name = "Array")

density_plot(merge(commonTg4510Arraydf, commonTg4510RRBSdf, by = 0),x.var = "RRBS", y.var = "Array", x_lab = "RRBS", y_lab = "Array", title = "") + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10")


