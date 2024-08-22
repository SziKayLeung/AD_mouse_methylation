
rootDir = "/lustre/projects/Research_Project-191406"
scriptDir = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_mouse_methylation/"
source(paste0(scriptDir, "import.config"))
source(paste0(scriptDir, "0_PaperFigs/Functions.R"))
source(paste0(scriptDir, "0_PaperFigs/paper_import.config.R"))

library("stringr")
library("ggrepel")
library("dplyr")

manhattanPlot <- rTg4510_rrbs_results$Genotype %>% dplyr::select(Position, p.val.Genotype) 
rTg4510_array_results <- get(load(paste0(dirnames$differential, "/array/rTg4510_array_ECX_allResultsDMPs.RData")))
manhattanPlotArray <- rTg4510_array_results$Genotype 

prepare_manhattan <- function(rrbs, array){
  
  prep_coordinates <- function(manhattanPlot){
    manhattanPlot$CHR <- stringr::word(manhattanPlot$Position,c(1),sep=stringr::fixed(":"))
    manhattanPlot$BP <- stringr::word(manhattanPlot$Position,c(2),sep=stringr::fixed(":"))
    manhattanPlot[,"CHR"][which(manhattanPlot[,"CHR"] == "chrX")]<-"chr23"
    manhattanPlot[,"CHR"][which(manhattanPlot[,"CHR"] == "chrY")]<-"chr24"
    manhattanPlot$CHR <- as.numeric(str_remove(manhattanPlot$CHR,"chr"))
    manhattanPlot$BP <- as.numeric(manhattanPlot$BP)
    return(manhattanPlot)
  }
  
  # rrbs
  rrbs <- rrbs %>% dplyr::select(Position, p.val.Genotype) 
  rrbsPrepared <- prep_coordinates(rrbs) %>% mutate(Platform = "RRBS")
  array <- array %>% 
    dplyr::select(X,PrZ.GenotypeTG) %>% 
    merge(., mm10_Manifest, by.x = "X", by.y = 0) %>% 
    mutate(Position = position)
  arrayPrepared <- prep_coordinates(array) %>% mutate(Platform = "Array") %>% dplyr::select(position, PrZ.GenotypeTG, CHR, BP, Platform)
  colnames(arrayPrepared) <- colnames(rrbsPrepared)
  
  mergedPrepared <- rbind(rrbsPrepared,arrayPrepared)
  mergedPrepared <- mergedPrepared %>% mutate(BP = as.numeric(BP))
  
  don <- mergedPrepared %>% 
    
    # Compute chromosome size
    dplyr::group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    dplyr::mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(mergedPrepared, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot)
  
  return(don)
}


plot_manhattan <- function(don, sigResRRBS, sigResArray){
  axisdf = don %>%
    group_by(CHR) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  # merge with significant values to label
  donRRBS <- merge(don[don$Platform == "RRBS",], sigResRRBS[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position", all = T)
  donArray <- merge(don[don$Platform == "Array",], sigResArray[,c("position","ChIPseeker_GeneSymbol")], by.x = "Position", by.y = "position", all = T)
  don <- rbind(donRRBS, donArray)
  don <- don %>% arrange(p.val.Genotype)
  top <- unique(don$ChIPseeker_GeneSymbol)[1:200]
  don <- don %>% mutate(label = ifelse(ChIPseeker_GeneSymbol %in% top, ChIPseeker_GeneSymbol, NA))
  don <- don %>% mutate(label_unique = if_else(duplicated(label), NA_character_, label))
  
  # remove random chromosomes
  don <- don %>% filter(!is.na(CHR))
  p <- ggplot(don, aes(x=BPcum, y=-log10(p.val.Genotype))) +
    
    # Show all points
    geom_point( aes(color=as.factor(Platform)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("red", "skyblue"), 22 )) +
    
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
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(x = "Chromosome", y = expression(-log[10](italic(p))))
  
  return(p)
  
}


preppedManhattan <- list(
  rTg4510Genotype = prepare_manhattan(rTg4510_rrbs_results$Genotype,rTg4510_array_results$Genotype),
)
preppedManhattan$rTg4510Pathology = prepare_manhattan(rTg4510_rrbs_results$,rTg4510_array_results$Genotype)
pManhattanTg4510Geno <- plot_manhattan(preppedManhattan$rTg4510Genotype, rTg4510_rrbs_anno$Genotype, rTg4510_array_anno$ECX$Genotype)

