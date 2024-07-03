suppressMessages(library(ChIPseeker))

label_colour <- function(var){
  if(var %in% c("rTg4510","Tg4510")){colour = "#00AEC9"}else{
    if(var == "J20"){colour = "#FF5A62"}else{
    }}
  return(colour)
}


##-------------- PrepMixedEffectsStats

# Aim: merge array statistics results with manifest, calculate FDR and order the FDR interactive term
# Input:
  # mixedEffectsResults = df from running the statistics model
  # mm10_Mannifest = read in manifest file
  # inTerm = str: <Tissue.Genotype> <Geno.Age> = the additional p.value to calculate FDR and further ordered
# Output: 
  # mixedEffectsResults = merged df of probe annotation and additional FDR columns

PrepMixedEffectsStats <- function(mixedEffectsResults,manifest,test){
  
  mixedEffectsResults <- merge(manifest, mixedEffectsResults, by.x = "cpg", by.y = "X")
  mixedEffectsResults <-mixedEffectsResults %>% mutate(Position = paste0(seqnames,":",start))
  
  
  if(test == "Tissue.Genotype"){
    mixedEffectsResults$FDR_adj_TissueGenotype <- p.adjust(mixedEffectsResults[,"PrZ.Tissue.GenotypeTG"], method = "fdr")
    mixedEffectsResults <-mixedEffectsResults[order(mixedEffectsResults$PrZ.Tissue.GenotypeTG),]
    
    mixedEffectsResults <- mixedEffectsResults[which(mixedEffectsResults$FDR_adj_TissueGenotype < 0.05),]
    cat("Number of significant sites in Tissue*Genotype", nrow(mixedEffectsResults), "\n") 
    
  }else if(test == "Geno.Age"){
    mixedEffectsResults$FDR_adj_GenotypeAge <- p.adjust(mixedEffectsResults[,"PrZ.GenotypeTG.Age_months"], method = "fdr")
    mixedEffectsResults$FDR_adj_Genotype <- p.adjust(mixedEffectsResults[,"PrZ.GenotypeTG"], method = "fdr")
    
    mixedEffectsResults <- mixedEffectsResults[which(mixedEffectsResults$FDR_adj_GenotypeAge < 0.05),]
    cat("Number of significant sites in Age*Genotype", nrow(mixedEffectsResults), "\n") 
  
  }else if(test == "Geno"){
    mixedEffectsResults$FDR_adj_GenotypeAge <- p.adjust(mixedEffectsResults[,"PrZ.GenotypeTG.Age_months"], method = "fdr")
    mixedEffectsResults$FDR_adj_Genotype <- p.adjust(mixedEffectsResults[,"PrZ.GenotypeTG"], method = "fdr")
    
    mixedEffectsResults <- mixedEffectsResults[which(mixedEffectsResults$FDR_adj_Genotype < 0.05),]
    cat("Number of significant sites in Genotype", nrow(mixedEffectsResults), "\n") 
    
  }else if(test == "Pathology"){
    mixedEffectsResults$FDR_adj_Pathology <- p.adjust(mixedEffectsResults[,"PrZ.Pathology"], method = "fdr")
    mixedEffectsResults <- mixedEffectsResults[which(mixedEffectsResults$FDR_adj_Pathology < 0.05),]
    cat("Number of significant sites in Pathology", nrow(mixedEffectsResults), "\n") 
    
    
  }else{
    print("Interaction term required either: <Tissue.Genotype> <Geno.Age> <Geno> <Pathology>")
  }
  
  return(mixedEffectsResults)
}


##-------------- PrepBetasPhenoFiles

# Aim: prepare the beta values from the normalised array, QC 
# code taken from E.Walker 
# Input:
  # model_path = df of phenotype array data 
  # mouseModel = str: <rTg4510, J20>
# Output: 
# combined df of phneotype and beta values

PrepBetasPhenoFiles <- function(model_path, mouseModel){
  
  load("/gpfs/ts0/projects/Research_Project-191406/Aisha/data/Array/Normalised_Data_Sesame_2.rdat")
  QCmetrics$SampleID <- gsub(".*_", "", QCmetrics$ExternalSampleID)
  ## ! The same samples do not have the same sample ID, these are matched on the pathology file. If we get the 
  ## same order of the sample ID as the pathology file then they are matched correctly
  ### ECX
  pheno_model_ECX <- QCmetrics[which(QCmetrics$Tissue == "CortexEntorhinalis" & QCmetrics$AD_model == mouseModel),]
  betas_model_ECX <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno_model_ECX$Basename ]
  #identical(pheno_model_ECX$Basename, colnames(betas_model_ECX))
  colnames(betas_model_ECX) <- pheno_model_ECX$SampleID
  rownames(pheno_model_ECX) <- pheno_model_ECX$SampleID 
  pheno_model_ECX <- pheno_model_ECX[model_path$Sample_ID_ECX, ] #order based on pathology samples id (NA rows will be generated)
  pheno_model_ECX <- cbind(pheno_model_ECX, model_path)
  ### Hip
  pheno_model_HIP <- QCmetrics[which(QCmetrics$Tissue == "Hippocampus" & QCmetrics$AD_model == mouseModel),]
  betas_model_HIP <- Normalised_Sesame_Betas[,colnames(Normalised_Sesame_Betas) %in% pheno_model_HIP$Basename ]
  #identical(pheno_model_HIP$Basename, colnames(betas_model_HIP))
  colnames(betas_model_HIP) <- pheno_model_HIP$SampleID
  rownames(pheno_model_HIP) <- pheno_model_HIP$SampleID 
  pheno_model_HIP <- pheno_model_HIP[model_path$Sample_ID_HIP, ] #order based on pathology samples id (NA rows will be generated)
  pheno_model_HIP <- cbind(pheno_model_HIP, model_path)
  
  # Match the same mouse with 2 brain regions
  samples <- as.data.frame(cbind(pheno_model_ECX$SampleID, pheno_model_HIP$SampleID))
  colnames(samples) <- c("ECX", "HIP")
  # samples <- samples[!is.na(samples$ECX),]
  # samples <- samples[!is.na(samples$HIP),]
  samples$ECX <- as.character(samples$ECX)
  samples$HIP <- as.character(samples$HIP)
  for(i in 1:nrow(samples)){
    samples[i, "MouseID"] <- paste("M", i, sep="")
  }
  ## 36 matched samples!!
  
  # Now filter these matching samples on both tissues
  ### ecx
  #pheno_model_ECX <- pheno_model_ECX[samples$ECX,]
  #identical(pheno_model_ECX$SampleID, samples$ECX)
  pheno_model_ECX$MouseID <- samples$MouseID[match(samples$ECX, pheno_model_ECX$SampleID)] #Now the orders are the same - add in the mouseIDs
  #betas_model_ECX <- betas_model_ECX[,pheno_model_ECX$SampleID]
  #identical(colnames(betas_model_ECX) , pheno_model_ECX$SampleID)
  pheno_model_ECX <- pheno_model_ECX[colnames(betas_model_ECX), ]
  #identical(colnames(betas_model_ECX) , pheno_model_ECX$SampleID)
  
  
  ### hip
  #pheno_model_HIP <- pheno_model_HIP[samples$HIP,]
  #identical(pheno_model_HIP$SampleID, samples$HIP)
  pheno_model_HIP$MouseID <- samples$MouseID[match(samples$HIP, pheno_model_HIP$SampleID)] #Now the orders are the same - add in the mouseIDs
  #betas_model_HIP <- betas_model_HIP[,pheno_model_HIP$SampleID]
  #identical(colnames(betas_model_HIP) , pheno_model_HIP$SampleID)
  pheno_model_HIP <- pheno_model_HIP[colnames(betas_model_HIP), ]
  #identical(colnames(betas_model_HIP) , pheno_model_HIP$SampleID)
  
  
  pheno <- rbind(pheno_model_ECX, pheno_model_HIP)
  betas <- cbind(betas_model_ECX, betas_model_HIP)
  if(!identical(pheno$SampleID, colnames(betas))){
    print("betas and pheno not in same order")
  }
  
  pheno$Genotype <- factor(pheno$Genotype, levels = c("WT","TG"))
  pheno$Tissue <- factor(pheno$Tissue, levels = c("CortexEntorhinalis","Hippocampus"))
  pheno$Chip_ID <- as.factor(pheno$Chip_ID)
  pheno$MouseID <- as.factor(pheno$MouseID)
  
  pheno <- pheno[,c("Tissue","Genotype", "Age_months", "Chip_ID", "PlateNumber", "MouseID","Pathology_ECX", "Pathology_HIP")]
  betas <- as.matrix(betas)
  
  #combine pathology into one column
  pheno$Pathology <- rep(NA, nrow(pheno))
  for(i in 1:nrow(pheno)){
    pheno[i,"Pathology"] <- ifelse(pheno[i,"Tissue"] == 'CortexEntorhinalis', pheno[i,"Pathology_ECX"], 
                                   pheno[i,"Pathology_HIP"])
  }
  pheno <- pheno[,c("Tissue","Genotype", "Age_months", "Chip_ID", "PlateNumber", "MouseID", "Pathology")]
  
  output <- list(pheno, betas)
  names(output) <- c("pheno","betas")
  
  return(output)
  
}


##-------------- plotArrayDMP

# Aim: plot the DMP from array data
# Input:
  # cpg = str: cpg site
  # betaResults = output from linear model 
  # betas = beta values 
  # pheno = phenotype values
  # colour = rTg4510 or J20
  # column = column with the fdr p value to report
  # feature = <Tissue,GenotypeECX> to dictate appearnce of plot (facet to include hippocampus or to filter for ECX data)
# Output: 
  # box-plot of the methylation values across age

plotArrayDMP <- function(cpg, betaResults, betas, pheno, colour, column, feature){
  data <- betas
  coldata <- pheno
  index <- cpg
  ages <- pheno$Age_months
  data_site <- reshape2::melt(data[cpg,]) %>% tibble::rownames_to_column(., var = "sample") 
  data_plot <- merge(coldata[,c("Age_months","Genotype", "Tissue")], data_site, by.x = 0, by.y = "sample")
  colnames(data_plot) <- c("sample","Age_months","Genotype", "Tissue", "site")
  data_plot <- data_plot[complete.cases(data_plot), ]
  data_betaResults <- betaResults[betaResults$cpg == cpg,]
  
  if(feature=="Tissue"){
    data_plot$Tissue <- factor(data_plot$Tissue, levels = c("CortexEntorhinalis", "Hippocampus"),
                               labels = c("Entorhinal Cortex", "Hippocampus"))
  }else if(feature=="GenotypeECX"){
    data_plot <- data_plot %>% filter(Tissue == "CortexEntorhinalis")
  }else{
    print("feature required: <Tissue> <GenotypeECX>")
  }

  p <- ggplot(data_plot, aes(x=factor(Age_months), y=site, color=Genotype)) + 
    geom_boxplot() + 
    geom_point(aes(fill = Genotype), size = 1.5, shape = 21, position = position_jitterdodge()) +
    facet_wrap(~Tissue) +
    scale_color_manual(values=c("black", colour)) +
    scale_fill_manual(values=c("black", colour)) +
    labs(y = "Methylation", x = "Age (months)", title = paste(data_betaResults$Gene_Symbol, betaResults$cpg, sep = " - "),
         subtitle =paste("FDR Pvalue: ", signif(data_betaResults[[column]],3), sep = ""))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          strip.background = element_rect(colour="white",fill="white"),
          strip.text.x = element_text(size = 10),
          legend.position = "None")
  
  if(feature=="Tissue"){
    p <- p +  facet_wrap(~Tissue)
  }
  
  return(p)
}
