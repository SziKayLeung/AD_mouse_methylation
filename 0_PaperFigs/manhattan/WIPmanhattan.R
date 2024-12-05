library("stringr")
library("ggrepel")

manhattanPlot <- rTg4510_rrbs_results$Genotype %>% dplyr::select(Position, p.val.Genotype) 
manhattanPlot$CHR <- stringr::word(manhattanPlot$Position,c(1),sep=stringr::fixed(":"))
manhattanPlot$BP <- stringr::word(manhattanPlot$Position,c(2),sep=stringr::fixed(":"))
manhattanPlot[,"CHR"][which(manhattanPlot[,"CHR"] == "chrX")]<-"chr23"
manhattanPlot[,"CHR"][which(manhattanPlot[,"CHR"] == "chrY")]<-"chr24"
manhattanPlot$CHR <- as.numeric(str_remove(manhattanPlot$CHR,"chr"))
manhattanPlot$BP <- as.numeric(manhattanPlot$BP)
manhattanPlot <- manhattanPlot %>% mutate(Platform = "RRBS")
#manhattanPlot <- merge(manhattanPlot, rTg4510_rrbs_anno$Genotype[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position") 

rTg4510_array_results <- get(load(paste0(dirnames$differential, "/array/rTg4510_array_ECX_allResultsDMPs.RData")))
manhattanPlotArray <- rTg4510_array_results$Genotype %>% dplyr::select(X,PrZ.GenotypeTG) %>% 
  merge(., mm10_Manifest, by.x = "X", by.y = 0)
manhattanPlotArray$CHR <- stringr::word(manhattanPlotArray$position,c(1),sep=stringr::fixed(":"))
manhattanPlotArray$BP <- stringr::word(manhattanPlotArray$position,c(2),sep=stringr::fixed(":"))
manhattanPlotArray[,"CHR"][which(manhattanPlotArray[,"CHR"] == "chrX")]<-"chr23"
manhattanPlotArray[,"CHR"][which(manhattanPlotArray[,"CHR"] == "chrY")]<-"chr24"
manhattanPlotArray$CHR <- as.numeric(str_remove(manhattanPlotArray$CHR,"chr"))
manhattanPlotArray$BP <- as.numeric(manhattanPlotArray$BP)
manhattanPlotArray <- manhattanPlotArray  %>% mutate(Platform = "Array")
manhattanPlotArray <- manhattanPlotArray %>% dplyr::select(position, PrZ.GenotypeTG, CHR, BP, Platform)
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

donRRBS <- merge(don[don$Platform == "RRBS",], 
                 rTg4510_rrbs_anno$Genotype[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position", all = T)
donArray <- merge(don[don$Platform == "Array",], 
                  rTg4510_array_anno$ECX$Genotype[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position", all = T)
don <- rbind(donRRBS, donArray)
don <- don %>% arrange(p.val.Genotype)
top <- unique(don$ChIPseeker_GeneSymbol)[1:200]
don <- don %>% mutate(label = ifelse(ChIPseeker_GeneSymbol %in% top, ChIPseeker_GeneSymbol, NA))
don <- don %>% mutate(label_unique = if_else(duplicated(label), NA_character_, label))

don <- don %>% filter(!is.na(CHR))
don$Platform <- as.factor(don$Platform)
p <- ggplot(don, aes(x=BPcum, y=-log10(p.val.Genotype))) +
  
  # Show all points
  geom_point(aes(color = Platform, fill = as.factor(Platform)), 
             size = 1.3, shape = 21) +  # Adjust alpha to 0.6 for better visibility
  scale_color_manual(values = rep(c(wes_palette("Rushmore1")[4], alpha(wes_palette("Rushmore1")[3],0.2)), 22)) +
  scale_fill_manual(values = rep(c(wes_palette("Rushmore1")[4], alpha(wes_palette("Rushmore1")[3],0.2)), 22)) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  geom_text_repel(data = don[!is.na(don$label_unique), ], 
                  aes(label = label_unique), 
                  size = 3, 
                  box.padding = 0.3, 
                  point.padding = 0.5) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="bottom",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Chromosome", y = expression(-log[10](italic(p))))

pdf(paste0(dirnames$paper, "/ManhattanPlot_rTg4510genotype.pdf"), width = 6, height = 6)
p
dev.off()



manhattanPlot <- datawrangle_manhattan(rTg4510_rrbs_results$Genotype, "p.val.Genotype", "RRBS")
manhattanPlotArray <- datawrangle_manhattan(rTg4510_array_results$Genotype, "PrZ.GenotypeTG", "Array")
colnames(manhattanPlotArray) <- colnames(manhattanPlot)
mergedManhattanPlot <- rbind(manhattanPlotArray,manhattanPlot)
p <- plot_merged_manhattan(mergedManhattanPlot, rTg4510_rrbs_anno$Genotype,rTg4510_array_anno$ECX$Genotype)


