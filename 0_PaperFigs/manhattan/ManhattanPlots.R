library("stringr")
library("ggrepel")

plot_manhattan <- function(rrbs, array, mode){
  
  if(mode == "Genotype"){
    manhattanPlot <- rrbs %>% dplyr::select(Position, p.val.Genotype) 
    manhattanPlotArray <- array %>% dplyr::select(X,PrZ.GenotypeTG) %>% 
      merge(., mm10_Manifest, by.x = "X", by.y = 0)
  }else{
    manhattanPlot <- rrbs %>% dplyr::select(Position, p.val.Pathology) 
    manhattanPlotArray <- array %>% dplyr::select(X,PrZ.Pathology) %>% 
      merge(., mm10_Manifest, by.x = "X", by.y = 0)
  }
  manhattanPlot$CHR <- stringr::word(manhattanPlot$Position,c(1),sep=stringr::fixed(":"))
  manhattanPlot$BP <- stringr::word(manhattanPlot$Position,c(2),sep=stringr::fixed(":"))
  manhattanPlot[,"CHR"][which(manhattanPlot[,"CHR"] == "chrX")]<-"chr23"
  manhattanPlot[,"CHR"][which(manhattanPlot[,"CHR"] == "chrY")]<-"chr24"
  manhattanPlot$CHR <- as.numeric(str_remove(manhattanPlot$CHR,"chr"))
  manhattanPlot$BP <- as.numeric(manhattanPlot$BP)
  manhattanPlot <- manhattanPlot %>% mutate(Platform = "RRBS")
  #manhattanPlot <- merge(manhattanPlot, rTg4510_rrbs_anno$Genotype[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position") 
  
  manhattanPlotArray$CHR <- stringr::word(manhattanPlotArray$position,c(1),sep=stringr::fixed(":"))
  manhattanPlotArray$BP <- stringr::word(manhattanPlotArray$position,c(2),sep=stringr::fixed(":"))
  manhattanPlotArray[,"CHR"][which(manhattanPlotArray[,"CHR"] == "chrX")]<-"chr23"
  manhattanPlotArray[,"CHR"][which(manhattanPlotArray[,"CHR"] == "chrY")]<-"chr24"
  manhattanPlotArray$CHR <- as.numeric(str_remove(manhattanPlotArray$CHR,"chr"))
  manhattanPlotArray$BP <- as.numeric(manhattanPlotArray$BP)
  manhattanPlotArray <- manhattanPlotArray  %>% mutate(Platform = "Array")
  
  if(mode == "Genotype"){
    manhattanPlotArray <- manhattanPlotArray %>% dplyr::select(position, PrZ.GenotypeTG, CHR, BP, Platform)
  }else{
    manhattanPlotArray <- manhattanPlotArray %>% dplyr::select(position, PrZ.Pathology, CHR, BP, Platform)
  }
  colnames(manhattanPlotArray) <- colnames(manhattanPlot)
  
  mergedManhattanPlot <- rbind(manhattanPlotArray,manhattanPlot)
  
  don <- mergedManhattanPlot %>% 
    
    # Compute chromosome size
    dplyr::group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(mergedManhattanPlot, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot)
  
  axisdf = don %>%
    group_by(CHR) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) 
  axisdf <- axisdf %>% filter(CHR != "NA")
  
  donRRBS <- merge(don[don$Platform == "RRBS",], 
                   rTg4510_rrbs_anno$Genotype[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position", all = T)
  donArray <- merge(don[don$Platform == "Array",], 
                    rTg4510_array_anno$ECX$Genotype[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position", all = T)
  don <- rbind(donRRBS, donArray)
  if(mode == "Genotype"){
    don <- don %>% arrange(p.val.Genotype)
  }else{
    don <- don %>% arrange(p.val.Pathology)
  }
  top <- unique(don$ChIPseeker_GeneSymbol)[1:200]
  don <- don %>% mutate(label = ifelse(ChIPseeker_GeneSymbol %in% top, ChIPseeker_GeneSymbol, NA))
  don <- don %>% mutate(label_unique = if_else(duplicated(label), NA_character_, label))
  
  don <- don %>% filter(!is.na(CHR))
  don$Platform <- as.factor(don$Platform)
  if(mode == "Genotype"){
    p <- ggplot(don, aes(x=BPcum, y=-log10(p.val.Genotype))) 
  }else{
    p <- ggplot(don, aes(x=BPcum, y=-log10(p.val.Pathology))) 
  }
  p <- p +
    
    # Show all points
    geom_point(aes(color = Platform, fill = as.factor(Platform)), 
               size = 1.3, shape = 21) +  # Adjust alpha to 0.6 for better visibility
    scale_color_manual(values = rep(c(wes_palette("Rushmore1")[4], alpha(wes_palette("Rushmore1")[3],0.2)), 22)) +
    scale_fill_manual(values = rep(c(wes_palette("Rushmore1")[4], alpha(wes_palette("Rushmore1")[3],0.2)), 22)) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    #geom_text_repel(data = don[!is.na(don$label_unique), ], 
    #                aes(label = label_unique), 
    #                size = 3, 
    #                box.padding = 0.3, 
    #                point.padding = 0.5) +
    
    # Custom the theme:
    mytheme +
    theme(legend.position="none") +
    labs(x = "Chromosome", y = expression(-log[10](italic(p))))
  
  return(p)
  
}


# export plots as png (width = 1200, height = 260)
pManhattan <- list()
pManhattan$rTg4510_Genotype <- plot_manhattan(rTg4510_rrbs_results$Genotype, rTg4510_array_results$Genotype, mode = "Genotype")
pManhattan$rTg4510_Genotype
pManhattan$rTg4510_Pathology <- plot_manhattan(rTg4510_rrbs_results$Pathology, rTg4510_array_results$Pathology, mode = "Pathology")
pManhattan$rTg4510_Pathology
pManhattan$J20_Genotype <- plot_manhattan(J20_rrbs_results$Genotype, J20_array_results$Genotype, mode = "Genotype")
pManhattan$J20_Genotype 
pManhattan$J20_Pathology <- plot_manhattan(J20_rrbs_results$Pathology, J20_array_results$Pathology, mode = "Pathology")

