RNAcorMethPlot <- function(res, ModelTerm, RNA, RRBS, colData, colours, probe, MouseModel, stats = FALSE){
  
  # res is significant DMP results file from regression
  # RNA is normalised counts for the mouse model
  # RRBS is RRBS methyaltion matrix for the mouse model
  # colData is pheno data for the mouse model
  # colours is the colours for TG/WT
  # gene_number refers to signficance level. e.g. 1 is top most significant, 2 is 2nd most significant etc N.B only includes sites annotated to a          gene also in the RNA data
  
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  
  colnames(RNA)[1] <- "Gene" 
  
  pvalcol <- paste0("FDR_adj_", ModelTerm)
  
  
  #susbet DMP results to those in gene of interest
  res <- res[order(res[, pvalcol]),]# order results file by pval
  res_sub <- res[which(res$Position == probe),] #remove rows where DMPs weren't annotated to a gene in the RNA data
  gene <- res_sub$ChIPseeker_GeneSymbol
  #res_sub <- res[which(res$ChIPseeker_GeneSymbol == unique(res$Genes_3000bp_TSS)[gene_number]),] # subset to sites annotated to gene_number
  
  #Subset RRBS data to sites in res_sub
  RRBS_sub <- RRBS[rownames(RRBS) %in% res_sub$Position,]
  RRBS_sub <- as.data.frame(t(RRBS_sub))
  RRBS_sub$Sample_ID <- rownames(RRBS_sub)
  
  # subset RNA data to gene in res_sub
  RNA_sub <- RNA[which(RNA$Gene == gene),]
  RNA_sub <- as.data.frame(t(RNA_sub))
  colnames(RNA_sub) <- "RNAexpr"
  RNA_sub$Sample_ID <- rownames(RNA_sub)
  RNA_sub <- RNA_sub[-c(1),]
  RNA_sub$Sample_ID <- rownames(RNA_sub)
  RNA_sub$RNAexpr <- as.numeric(as.character(RNA_sub$RNAexpr))
  
  
  ########################## plot ====================
  
  # create dataframe containing all data for plot
  plotdf <- left_join(colData %>% dplyr::select(Sample_ID, Genotype, Age_months), RNA_sub, by = "Sample_ID") %>%
    left_join(., RRBS_sub, by = "Sample_ID")
  plotdf <- melt(plotdf, id.vars = c("Sample_ID", "Genotype", "Age_months","RNAexpr")) %>% dplyr::rename("Methylation" = "value")
  plotdf <- plotdf[complete.cases(plotdf),]
  
  if(isFALSE(stats)){
    # create plot
    p <- ggplot(plotdf, aes(x = Methylation, y = RNAexpr))+
      geom_point(aes(color = Genotype, shape = as.factor(Age_months)), size = 2) +
      geom_smooth(method=lm, aes(color=Genotype), formula = y ~ x)+
      #geom_smooth(method=lm, se=FALSE, fullrange=FALSE, linetype = "dashed", formula = y ~ x)+
      scale_color_manual(values=colours)+
      labs(x = "Methylation", y = "RNA-Seq normalized counts", shape = "Age (months)") +
      ggtitle(paste0(gene,", ", probe)) + theme_classic()
    
    return(p)
  }else{
    return(plotdf)
  }
}

RNAcorMethPlotGene <- function(res, ModelTerm, RNA, RRBS, colData, gene, MouseModel, stats = FALSE){
  print(gene)
  sub_res <- res %>% filter(ChIPseeker_GeneSymbol == gene)
  output<- list()
  for(i in 1:length(sub_res$Position)){
    output[[i]] <- RNAcorMethPlot(res = sub_res, ModelTerm = ModelTerm, RNA = RNA, 
                                  RRBS = RRBS, colData = colData, colours = c(color_rTg4510_TG, WT), 
                                  probe = sub_res$Position[i], MouseModel = MouseModel, stats = stats)
  }
  return(output)
}

findCommonGene <- function(res_DMP, res_DGE, ModelTerm, RNA, RRBS, colData, colours, MouseModel){
  #find which genes are sig in both RNA and methylation data
  inBoth <- intersect(res_DMP$ChIPseeker_GeneSymbol, res_DGE$Gene)
  resDMPCommon <- res_DMP %>% filter(res_DMP$ChIPseeker_GeneSymbol %in% inBoth)
  print(nrow(resDMPCommon))
  if(ModelTerm == "pathology"){
    resDMPCommon <- resDMPCommon %>% filter(Promoter == "TRUE")
    print(nrow(resDMPCommon))
  }
  
  message("Total number of common genes", length(inBoth))
  message("Common genes: ", paste0(inBoth,","))
  if(ModelTerm == "genotype"){
    Output <- lapply(inBoth, function(gene) RNAcorMethPlotGene(res_DMP, ModelTerm, RNA, RRBS, colData, gene, MouseModel)) 
  }else{
    Output <- lapply(inBoth, function(gene) RNAcorMethPlotGeneInteraction(res_DMP, ModelTerm, RNA, RRBS, colData, gene, MouseModel)) 
  }
  return(Output)
}

RNAcorMethPlotGeneInteraction <- function(res, ModelTerm, RNA, RRBS, colData, gene, MouseModel,stats=TRUE){
  dat <- RNAcorMethPlotGene(res, ModelTerm, RNA,RRBS,colData,gene,MouseModel,stats)
  dat <- lapply(dat, function(x) x %>% mutate(rep = paste0(Genotype,Age_months)))
  
  output <- list()
  for(i in 1:length(dat)){
    probe <- dat[[i]][["variable"]]
    output[[i]] <- ggplot(dat[[i]], aes(x = Methylation, y = RNAexpr))+
      #geom_point(data = df2, colour = "grey70") +
      geom_point(aes(color = Genotype), size = 2) + 
      facet_grid(~Age_months) +
      geom_smooth(aes(group=Genotype,colour=Genotype), linetype="dashed", method="lm",formula = y ~ x) +
      labs(x = "Methylation", y = "RNA-Seq normalized counts", shape = "Age (months)") +
      ggtitle(paste0(gene,", ",probe)) +
      theme_classic()
  }
  
  return(output)
}

RNAcorMethPlotGenePathology <- function(res, ModelTerm, RNA, RRBS, colData, gene, MouseModel,stats=TRUE){
  dat <- RNAcorMethPlotGene(res, ModelTerm, RNA,RRBS,colData,gene,MouseModel,stats)
  dat <- lapply(dat, function(x) merge(x, colData[,c("Sample_ID","ECX")], by = "Sample_ID"))
  dat <- lapply(dat, function(x) x %>% mutate(rep = paste0(Genotype,Age_months)))
  dat <- lapply(dat, function(x) x %>% filter(Genotype != "WT"))
  
  output <- list()
  for(i in 1:length(dat)){
    probe <- dat[[i]][["variable"]]
    output[[i]] <- ggplot(dat[[i]], aes(ECX, as.numeric(RNAexpr))) + 
      geom_point(aes(colour=Methylation),size=3) +
      geom_smooth(method=lm, formula = y ~ x, linetype = "dashed", size = 0.5) +
      #facet_grid(~as.factor(Genotype)) +
      scale_colour_viridis() +
      theme_classic() +
      ggtitle(paste0(gene,", ",probe)) +
      labs(x = "Pathology", y = "RNA-Seq normalized counts")
  }
  return(output)
}

