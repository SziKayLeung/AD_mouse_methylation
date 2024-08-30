LOGEN_ROOT="/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen"
source(paste0(LOGEN_ROOT, "/aesthetics_basics_plots/draw_density.R"))
source(paste0(LOGEN_ROOT, "/aesthetics_basics_plots/pthemes.R"))

suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))

# common probes between beta dataset and number of probes from array
commonStatsDescription <- function(rrbsBeta, arrayBeta){
  message("Number of probes in RRBS after smoothing: ", length(row.names(rrbsBeta)))
  message("Number of probes in array detected (Horvath mouse array): ", length(row.names(arrayBeta)))
  common_probes <- intersect(row.names(rrbsBeta), row.names(arrayBeta))
  message("Number of common probes: ", length(common_probes))
  message("% common probes in RRBS: ", round(length(common_probes)/length(row.names(rrbsBeta)) * 100, 2))
  message("% common probes in RRBS: ", round(length(common_probes)/length(row.names(arrayBeta)) * 100, 2))
  
  return(common_probes)
}

# correlation plot of the common probes between RRBS and Array
corrPlotCommonProbes <- function(rrbsBeta, arrayBeta, commonProbes){
  common <- list(
    rrbs = rrbsBeta[commonProbes,] %>% apply(., 1, mean, na.rm=TRUE) %>% reshape2::melt(value.name = "RRBS"),
    array = arrayBeta[commonProbes,] %>% apply(., 1, mean, na.rm=TRUE) %>% reshape2::melt(value.name = "Array")
  )
  
  pCommonCorr <- density_plot(merge(common$rrbs, common$array, by = 0),
                              x.var = "RRBS", y.var = "Array", x_lab = "RRBS - mean methylation", y_lab = "Array - mean methylation", title = "") + 
    scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10")
  
  return(pCommonCorr)
}


# calculate absolute difference between TG and WT
extract_delta <- function(betaMatrix, sigMatrix = FALSE, phenotypeFile, position = NULL){
  if(!isFALSE(sigMatrix)){
    dat <- betaMatrix %>% filter(row.names(betaMatrix) %in% sigMatrix[["position"]])
  }else{
    dat <- betaMatrix %>% filter(row.names(betaMatrix) %in% position)
  }
  
  
  if(!is.null(sigMatrix)){
    colnames(dat) <- unlist(lapply(colnames(betaMatrix), function(x) paste0(x,"_", phenotypeFile[row.names(phenotypeFile) == x,"Genotype"])))
    dat <- dat %>% tibble::rownames_to_column(., var = "position") %>% reshape2::melt(variable.name = "sample",value.name = "methylation", id = "position")
    dat <- dat %>% dplyr::mutate(Genotype = word(as.character(sample), c(2), sep = stringr::fixed("_")))
    dat <- dat %>% group_by(position, Genotype) %>% dplyr::summarise(mean = mean(methylation, na.rm = T), .groups = 'keep') 
    dat <- tidyr::spread(dat, Genotype, mean) %>% mutate(delta = TG - WT, 
                                                         directionTG = ifelse(delta < 0, "down", "up"))
    dat <- dat %>% dplyr::rename(meanTG = TG, meanWT = WT)
  }else{
    dat <- data.frame()
  }
  
  return(dat)
}


MergeArrayRRBSSigResults <- function(rrbsBeta, arrayBeta, rrbsSigResults, arraySigResults, phenotypeInput){
  
  # delta
  # calculate absolute difference between TG and WT --> delta = TG - WT
  # therefore negative delta, decrease methylation in TG
  message("Calculating delta difference between WT and TG")
  delta <- list(
    Genotype = list(
      array = extract_delta(betaMatrix = arrayBeta, sigMatrix = arraySigResults$ECX$Genotype, phenotypeFile = phenotypeInput),
      rrbs = extract_delta(betaMatrix = rrbsBeta, sigMatrix = rrbsSigResults$Genotype, phenotypeFile = phenotypeInput)
    ),
    Interaction = list(
      array = extract_delta(betaMatrix = arrayBeta, sigMatrix = arraySigResults$ECX$Interaction, phenotypeFile = phenotypeInput),
      rrbs = extract_delta(betaMatrix = rrbsBeta, sigMatrix = rrbsSigResults$Interaction, phenotypeFile = phenotypeInput)
    ),
    Pathology = list(
      array = extract_delta(betaMatrix = arrayBeta, sigMatrix = arraySigResults$ECX$Pathology, phenotypeFile = phenotypeInput),
      rrbs = extract_delta(betaMatrix = rrbsBeta, sigMatrix = rrbsSigResults$Pathology, phenotypeFile = phenotypeInput)
    )
  )
  
  ## merge array and RRBS results for genotype: ECX
  # arrange by P.value rather than FDR as merging array and RRBS results with different significant values
  # i.e. array would be more sensitive with higher FDR than RRBS given fewer sites (reduced multiple testing)
  cols2Keep <- c("meanWT", "meanTG", "delta", "directionTG", "ChIPseeker_GeneEnsembl","ChIPseeker_GeneSymbol", "ChIPseeker_GeneName","ChIPseeker_Annotation","ChIPseeker_TransEnsembl","distanceToTSS")
  
  # Process Genotype data
  message("Merge ECX Genotype")
  sig_ECX_Genotype <- rbind(
    mergeSigResults("Genotype", arraySigResults$ECX$Genotype, "Array", delta$Genotype$array, "PrZ.GenotypeTG", "FDR_adj_Genotype", cols2Keep),
    mergeSigResults("Genotype",rrbsSigResults$Genotype, "RRBS", delta$Genotype$rrbs, "p.val.Genotype", "FDR_adj_genotype", cols2Keep)
  ) %>% arrange(P.value_Genotype) 
  
  # Process GenotypeAge data
  message("Merge ECX GenotypeAge")
  sig_ECX_GenotypeAge <- rbind(
    mergeSigResults("GenotypeAge", arraySigResults$ECX$Interaction, "Array", delta$Interaction$array, "PrZ.GenotypeTG.Age_months", "FDR_adj_GenotypeAge", cols2Keep),
    mergeSigResults("GenotypeAge", rrbsSigResults$Interaction, "RRBS", delta$Interaction$rrbs, "p.val.Interaction", "FDR_adj_interaction", cols2Keep)
  ) %>% arrange(P.value_GenotypeAge)
  
  # Process Pathology data 
  message("Merge ECX Pathology")
  sig_ECX_Pathology <- rbind(
    mergeSigResults("Pathology", arraySigResults$ECX$Pathology, "Array", delta$Pathology$array, "PrZ.Pathology", "FDR_adj_Pathology", cols2Keep),
    mergeSigResults("Pathology", rrbsSigResults$Pathology, "RRBS", delta$Pathology$rrbs, "p.val.Pathology", "FDR_adj_pathology", cols2Keep)
  ) %>% arrange(P.value_Pathology)
  
  
  # overlap of probes that are genotype only vs interaction
  probe_acrossModels <- list(
    GenotypeInteraction = intersect(sig_ECX_GenotypeAge$Position, sig_ECX_Genotype$Position),
    GenotypeOnly = setdiff(sig_ECX_Genotype$Position, sig_ECX_GenotypeAge$Position),
    InteractionOnly = setdiff(sig_ECX_GenotypeAge$Position, sig_ECX_Genotype$Position),
    GenotypePathology = intersect(sig_ECX_Genotype$Position, sig_ECX_Pathology$Position)
  )
  
  sig_ECX <- list(
    Genotype = sig_ECX_Genotype %>% mutate(InteractionEffect = ifelse(Position %in% probe_acrossModels$GenotypeInteraction,TRUE,FALSE)),
    GenotypeAge = sig_ECX_GenotypeAge %>% mutate(GenotypeEffect = ifelse(Position %in% probe_acrossModels$GenotypeInteraction,TRUE,FALSE)),
    Pathology = sig_ECX_Pathology %>% mutate(GenotypeEffect = ifelse(Position %in% probe_acrossModels$GenotypePathology,TRUE,FALSE))
  )
  
  
  # merge beta dataset
  common_samples <- intersect(colnames(rrbsBeta), colnames(arrayBeta))
  sig_beta <- list(
    Genotype = rbind(arrayBeta %>% filter(row.names(.) %in% sig_ECX$Genotype$Position) %>% dplyr::select(all_of(common_samples)),
                     rrbsBeta %>% filter(row.names(.) %in% sig_ECX$Genotype$Position) %>% dplyr::select(all_of(common_samples))),
    Pathology = rbind(arrayBeta %>% filter(row.names(.) %in% sig_ECX$Pathology$Position) %>% dplyr::select(all_of(common_samples)),
                      rrbsBeta %>% filter(row.names(.) %in% sig_ECX$Pathology$Position) %>% dplyr::select(all_of(common_samples)))
  )
  
  output <- list(delta, sig_ECX_Genotype, sig_ECX_GenotypeAge, sig_ECX_Pathology, sig_ECX, sig_beta)
  names(output) <- c("delta","sig_ECX_Genotype", "sig_ECX_GenotypeAge", "sig_ECX_Pathology", "sig_ECX_Full", "sig_betaECX_Full")
  return(output)
}


sigArrayResults <- function(arrayBeta, arraySigResults, phenotypeInput){
  cols2Keep <- c("meanWT", "meanTG", "delta", "directionTG", "ChIPseeker_GeneEnsembl","ChIPseeker_GeneSymbol", "ChIPseeker_GeneName","ChIPseeker_Annotation","ChIPseeker_TransEnsembl","distanceToTSS")
  
  delta <- list(
    Genotype = extract_delta(betaMatrix = arrayBeta, sigMatrix = arraySigResults[["Genotype"]], phenotypeFile = phenotypeInput),
    GenotypeAge = extract_delta(betaMatrix = arrayBeta, sigMatrix = arraySigResults[["Interaction"]], phenotypeFile = phenotypeInput),
    Pathology = extract_delta(betaMatrix = arrayBeta, sigMatrix = arraySigResults[["Pathology"]], phenotypeFile = phenotypeInput)
  )
  
  sigRes <- list(
    Genotype = mergeSigResults("Genotype", arraySigResults[["Genotype"]], "Array", delta$Genotype, "PrZ.GenotypeTG", "FDR_adj_Genotype", cols2Keep),
    GenotypeAge = mergeSigResults("GenotypeAge", arraySigResults[["Interaction"]], "Array", delta$GenotypeAge, "PrZ.GenotypeTG.Age_months", "FDR_adj_GenotypeAge", cols2Keep),
    Pathology = mergeSigResults("Pathology", arraySigResults[["Pathology"]], "Array", delta$Pathology, "PrZ.Pathology", "FDR_adj_Pathology", cols2Keep)
  )
  
  
  sigRes$Pathology <- sigRes$Pathology %>% dplyr::mutate(
    GenotypeEffect = ifelse(Position %in% sigRes$Genotype$Position,TRUE,FALSE),
    InteractionEffect = ifelse(Position %in% sigRes$GenotypeAge$Position,TRUE,FALSE)
  )
  
  # remove duplicated array probes, taking the most significant probe as representative
  sigRes$Genotype <- sigRes$Genotype  %>% group_by(Position) %>% dplyr::slice(which.min(FDR_adj_Genotype))
  sigRes$GenotypeAge <- sigRes$GenotypeAge %>% group_by(Position) %>% dplyr::slice(which.min(FDR_adj_GenotypeAge))
  sigRes$Pathology <- sigRes$Pathology %>% group_by(Position) %>% dplyr::slice(which.min(FDR_adj_Pathology))
  
  sigRes <- lapply(sigRes, function(x) as.data.frame(x))
  return(sigRes)
}


merge_beta_phenotype <- function(betaMatrix, phenotypeFile, position){
  
  dat <- betaMatrix[position,]
  dat <- dat %>% tibble::rownames_to_column(., var = "position") %>% reshape2::melt(variable.name = "sample",value.name = "methylation", id = "position")
  dat <- merge(dat, phenotypeFile, by.y = 0, by.x = "sample")

  return(dat)  
}

plot_DMP <- function(betaMatrix, phenotypeFile, position, interaction = FALSE, pathology = FALSE, model = "rTg4510", dat = NULL){
  
  colourPoints <- label_colour(model)
  
  if(is.null(dat)){
    dat <- merge_beta_phenotype(betaMatrix, phenotypeFile, position)
  }else{
    dat <- dat
  }
  
  if(isTRUE(interaction)){
    p <- ggplot(dat, aes(x = Age_months, y = methylation)) + 
      geom_point(aes(colour = Genotype)) + 
      scale_fill_manual(values = c(alpha("black",0.2), colourPoints),guide="none") +
      stat_summary(data=dat, aes(x=Age_months, y=methylation, group=Genotype, colour=Genotype), fun ="mean", geom="line", linetype = "dotted") +
      labs(y = "Methylation", x = "Age (months)") + theme_classic() + 
      scale_colour_manual(values = c("black", colourPoints),guide="none") +
      theme(panel.border = element_rect(fill = NA, color = "grey", linetype = "dotted"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank()) 
  
  }else if(isTRUE(pathology)){
     p <- ggplot(dat, aes(x = ECX, y = methylation)) + 
      geom_point(aes(colour = Genotype)) + 
      scale_colour_manual(values = c(alpha("black",0.2), colourPoints),guide="none") +
      labs(x = "Pathology", y = "Methylation") +
      geom_smooth(method=lm, formula = y~poly(x,3),colour="black",fill = alpha("gray",0.8))

  }else{
    p <- ggplot(dat, aes(x = Genotype, y = methylation, fill = Genotype)) + geom_boxplot(outlier.shape = NA) +
      scale_fill_manual(values = c(alpha("black",0.2), colourPoints),guide="none") +
      scale_colour_manual(values = c("black", colourPoints),guide="none") +
      labs(x = "Genotype", y = "Methylation") 
  }
  
  if(length(position > 1)){
    p <- p + facet_grid(~position)
  }
  
  p <- p +
    theme_classic() + 
    theme(panel.border = element_rect(fill = NA, color = "grey", linetype = "dotted"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()) 

  
  return(p)
}

plot_DMP_byTissue <- function(ECXbetaMatrix, HIPbetaMatrix, ECXphenotypeFile, HIPphenotypeFile, position, 
                              interaction = FALSE, pathology = FALSE, model = "rTg4510", gene = NULL){
  
  ECX_dat <- merge_beta_phenotype(ECXbetaMatrix, ECXphenotypeFile, position) %>% mutate(tissue = "ECX")
  HIP_dat <- merge_beta_phenotype(HIPbetaMatrix, HIPphenotypeFile, position)%>% mutate(tissue = "HIP")
  dat <- rbind(ECX_dat, HIP_dat)

  p <- plot_DMP(betaMatrix=NULL,phenotypeFile=NULL,position=NULL, interaction=interaction,pathology=pathology,model=model,dat=dat) + facet_grid(~ tissue) 
  
  if(is.null(gene)){
    p <- p + labs(subtitle = position)
  }else{
    p <- p + labs(subtitle = paste0(gene, " (", position, ")"))
  }
  
  return(p)
}


mergeSigResults <- function(test, sigResults, platform_name, logFC_data, value_col, fdr_col, cols_to_keep) {
  
  if(is.null(sigResults)){
    dat <- data.frame(matrix(ncol = 5 + length(cols_to_keep), nrow = 0))
    colnames(dat) <- c("Platform", "Position", "cpg", paste0("P.value_",test), paste0("FDR_adj_",test), cols_to_keep)
    
  }else{
    if(platform_name == "Array"){
      sigResults <- sigResults %>% mutate(platform = platform_name)
    }else{
      sigResults <- sigResults %>% mutate(platform = platform_name, cpg = NA)
    }
    dat <- sigResults %>%
      mutate(position = Position) %>%
      merge(., logFC_data, by = "position", all = TRUE) %>%
      dplyr::select(platform, position, cpg, !!value_col, !!fdr_col, all_of(cols_to_keep)) %>%
      setNames(c("Platform", "Position", "cpg", paste0("P.value_",test), paste0("FDR_adj_",test), cols_to_keep))
  }
  
  
  return(dat)
}


plot_platform_commonProbes <- function(rrbsBeta, arrayBeta, SigResults, phenotypeInput, commonProbes){
  
  p <- rbind(
    extract_delta(betaMatrix=rrbsBeta,
                  phenotypeFile=phenotypeInput,
                  position=intersect(SigResults$Position, commonProbes)) %>% 
      mutate(platform = "RRBS"),
    extract_delta(betaMatrix=arrayBeta,
                  phenotypeFile=phenotypeInput,
                  position=intersect(SigResults$Position, commonProbes)) %>% 
      mutate(platform = "Array")
  ) %>% dplyr::select(-meanTG, -meanWT, -directionTG) %>% tidyr::spread(., platform, delta) %>% 
    ggplot(., aes(x = Array, y = RRBS)) + geom_point() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    mytheme + labs(x = "Array - delta", y = "RRBS - delta")
  
  return(p)
}


datawrangle_tissue_array <- function(res){
  res <- res %>% mutate(Platform = "Array") %>% 
    dplyr::select(Platform, position, cpg, PrZ.Tissue, PrZ.Tissue.GenotypeTG, FDR_adj_TissueGenotype,
                  ChIPseeker_GeneEnsembl, ChIPseeker_GeneSymbol, ChIPseeker_GeneName, ChIPseeker_Annotation, ChIPseeker_TransEnsembl, distanceToTSS) %>% 
    dplyr::rename("Position" = "position", "P.value_Tissue" = "PrZ.Tissue", "P.value_TissueGenotype" = "PrZ.Tissue.GenotypeTG") %>% 
    arrange(FDR_adj_TissueGenotype)
  
  return(res)
}
